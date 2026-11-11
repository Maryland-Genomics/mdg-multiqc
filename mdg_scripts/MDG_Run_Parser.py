import csv
import json
import logging
import re
import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
RUNINFO_COLUMNS_FILE = SCRIPT_DIR / "runinfo_columns.yaml"
ASSETS_DIR = SCRIPT_DIR / "assets"

Runinfo = dict[str, Any]
FieldValue = str | int | float

# Tables built from runinfo.json are tab-separated: MultiQC's custom content parser doesn't
# understand CSV quoting, so commas inside LIMS values (contact lists, comments) would split columns.
# get_field() strips tabs and line breaks from values, so no quoting is needed.
TSV_ARGS: dict[str, Any] = {"sep": "\t", "quoting": csv.QUOTE_NONE}


@dataclass(frozen=True)
class PlotImage:
    """One plot image and the report section that shows it"""

    section_id: str
    title: str
    source: Path
    # Name of the copy in the data folder's images/ folder
    file_name: str


@dataclass(frozen=True)
class Platform:
    """
    Everything that differs between sequencing platforms. A "folder" is an Unaligned folder
    (Illumina) or a cell folder (PacBio); each gets its own {outdir}/{folder}_data.
    """

    display_name: str
    # wrapper_multiqc.py option that picks a single folder
    folder_option: str
    find_folders: Callable[[Path, Runinfo, str, str | None], list[str]]
    # Lanes or cells holding the project, and the runinfo_columns.yaml field naming them
    get_units: Callable[[Runinfo, str], list[str]]
    units_key: str
    unit_field: str
    # MultiQC input folders read in place, besides the generated data folders
    module_search_dirs: Callable[[Path, str, str], list[Path]]
    tophits_dirs: Callable[[Path, str, str], list[Path]]
    tophits_row_id: Callable[[Path, str, str], str]
    plot_images: Callable[[Path, Runinfo, str, list[str]], list[PlotImage]]
    write_read_stats: Callable[[Path, Path, str], None]
    # Files covering the whole run, and files for one folder
    copy_run_files: Callable[[Path, Runinfo, Path, str, list[str]], None]
    copy_folder_files: Callable[[Path, Runinfo, Path, str, str], None]
    software_versions: dict[str, str]


def parse(run_path: Path, project: str, platform_name: str, folders: list[str], outdir: Path) -> None:
    """
    Generate the parsed data files MultiQC reads. Files covering the whole run go to
    {outdir}/run_data, and files specific to one folder go to {outdir}/{folder}_data.
    """
    platform = PLATFORMS[platform_name]
    runinfo = load_runinfo(run_path)
    columns = load_column_config()

    copy_assets(outdir)

    data_dir = get_run_data_dir(outdir)
    reset_dir(data_dir)
    write_work_order_info(runinfo, project, data_dir, columns["work_order"])
    write_sample_info(runinfo, project, platform, data_dir, columns["sample_info"][platform_name])
    write_sequencing_info(runinfo, platform, data_dir, columns["sequencing_info"][platform_name])
    platform.write_read_stats(run_path, data_dir, project)
    platform.copy_run_files(run_path, runinfo, data_dir, project, folders)
    tidy_tables(data_dir)

    # Tophits numbering continues across folders so report section ids stay unique
    tophits_count = 0
    for folder in folders:
        data_dir = get_folder_data_dir(outdir, folder)
        reset_dir(data_dir)
        tophits_count = write_megablast_qc(run_path, runinfo, platform, data_dir, project, folder, tophits_count)
        platform.copy_folder_files(run_path, runinfo, data_dir, project, folder)
        tidy_tables(data_dir)


def get_run_data_dir(outdir: Path) -> Path:
    return outdir / "run_data"


def get_folder_data_dir(outdir: Path, folder: str) -> Path:
    return outdir / f"{folder}_data"


def load_runinfo(run_path: Path) -> Runinfo:
    """
    Load the run's PA_Logs/runinfo.json
    """
    runinfo_json = run_path / "PA_Logs" / "runinfo.json"
    try:
        with open(runinfo_json, encoding="utf-8") as fin:
            return json.load(fin)
    except FileNotFoundError:
        raise FileNotFoundError(f"Run info file not found: {runinfo_json}") from None
    except json.JSONDecodeError as e:
        raise ValueError(f"Could not parse {runinfo_json}: {e}") from None


def load_column_config() -> dict[str, Any]:
    """
    Load the report table column definitions from runinfo_columns.yaml
    """
    with open(RUNINFO_COLUMNS_FILE, encoding="utf-8") as fin:
        return yaml.safe_load(fin)


def get_lanes(runinfo: Runinfo, project: str) -> list[str]:
    """
    Return the Illumina lanes (as strings, in lane order) holding this project
    """
    lanes = [lane for lane, info in runinfo["lanes"].items() if project in info["projects"]["value"]]
    return sorted(lanes, key=int)


def get_cells(runinfo: Runinfo, project: str) -> list[str]:
    """
    Return the PacBio cell folders (e.g. 1_A01) holding this project
    """
    return [cell for cell, info in runinfo["cells"].items() if project in info["projects"]["value"]]


def get_project_libraries(runinfo: Runinfo, project: str, include_controls: bool = False) -> dict[str, Any]:
    """
    Return {library_id: library} for this project's libraries. With include_controls, also
    return control libraries whose control_projects include this project.
    """
    return {
        lib_id: lib
        for lib_id, lib in runinfo["libraries"].items()
        if lib["project"]["value"] == project or (include_controls and project in lib["control_projects"]["value"])
    }


def get_field(
    entry: dict[str, Any], field: str | list[str], column: str, context: dict[str, str] | None = None
) -> FieldValue:
    """
    Return the value of one runinfo_columns.yaml field for a runinfo.json entry.

    field is a dot-separated path ("library_pool.index_set.index1"), a list of paths whose values
    are joined with "-", or a "_name" looked up in context (e.g. the row's lane or cell).
    """
    if isinstance(field, list):
        values = [get_field(entry, f, column, context) for f in field]
        return "-".join(str(v) for v in values if v != "")
    if field.startswith("_"):
        if context is None or field not in context:
            raise KeyError(f"Column '{column}': '{field}' is not available for this table")
        return context[field]

    def resolve(node: Any, keys: list[str]) -> FieldValue:
        # Lists (e.g. specimen_work_order): collect the value from every item. A single value
        # keeps its type, so a number stays a number.
        if isinstance(node, list):
            values = list(dict.fromkeys(v for v in (resolve(item, keys) for item in node) if v != ""))
            if len(values) == 1:
                return values[0]
            return "; ".join(str(v) for v in values)
        if not keys:
            if not (isinstance(node, dict) and "value" in node):
                raise ValueError(f"Column '{column}': runinfo.json field '{field}' is not a single value")
            value = node["value"]
            # The LIMS export writes "nan" for some empty fields (e.g. strain)
            if value is None or value == "nan":
                return ""
            if isinstance(value, list):
                value = ", ".join(str(v) for v in value)
            if isinstance(value, str):
                value = re.sub(r"\s*[\t\r\n]+\s*", " ", value).strip()
            return value
        if not isinstance(node, dict) or keys[0] not in node:
            raise KeyError(f"Column '{column}': runinfo.json has no field '{field}'")
        return resolve(node[keys[0]], keys[1:])

    return resolve(entry, field.split("."))


def quote_text(value: FieldValue) -> FieldValue:
    """
    Quote a non-empty text value from runinfo.json. MultiQC turns unquoted all-digit cells into
    numbers (dropping leading zeros from IDs like "00123" and adding thousands separators to kit
    lots), but leaves quoted cells as text. JSON numbers stay unquoted so MultiQC formats them as
    numbers, and empty values stay unquoted so empty columns are still dropped.
    """
    return f'"{value}"' if isinstance(value, str) and value != "" else value


def reset_dir(directory: Path) -> None:
    """
    Remove and recreate a generated data folder so files from a previous run can't linger
    """
    shutil.rmtree(directory, ignore_errors=True)
    directory.mkdir(parents=True)


def write_sequencing_info(runinfo: Runinfo, platform: Platform, data_dir: Path, columns: dict[str, Any]) -> None:
    """
    Writes the Sequencing Info table: one row describing the run, with the columns defined in
    runinfo_columns.yaml (paths relative to runinfo.json's "run" entry)
    """
    context = {"_platform": platform.display_name}
    row = {column: quote_text(get_field(runinfo["run"], field, column, context)) for column, field in columns.items()}
    df = pd.DataFrame([row], columns=list(columns))
    df.to_csv(data_dir / "sequencing_info.tsv", index=False, **TSV_ARGS)


def write_sample_info(
    runinfo: Runinfo, project: str, platform: Platform, data_dir: Path, columns: dict[str, Any]
) -> None:
    """
    Writes the Sample Information table: one row per project library per lane (Illumina) or
    cell (PacBio), with the columns defined in runinfo_columns.yaml. Control libraries whose
    control_projects include this project are listed too.
    """
    rows = []
    for unit in platform.get_units(runinfo, project):
        for lib_id in runinfo[platform.units_key][unit]["libraries"]["value"]:
            lib = runinfo["libraries"][lib_id]
            if lib["project"]["value"] != project and project not in lib["control_projects"]["value"]:
                continue
            context = {platform.unit_field: unit, "_library_id": lib_id}
            rows.append(
                {column: quote_text(get_field(lib, field, column, context)) for column, field in columns.items()}
            )

    df = pd.DataFrame(rows, columns=list(columns))
    df.to_csv(data_dir / "runinfo_clean.tsv", **TSV_ARGS)


def write_work_order_info(runinfo: Runinfo, project: str, data_dir: Path, columns: dict[str, Any]) -> None:
    """
    Writes the Work Order table: one row per work order on this project's libraries, with the
    columns defined in runinfo_columns.yaml plus a Sample Count of the project's libraries and
    its control libraries
    """
    rows = []
    for lib in get_project_libraries(runinfo, project).values():
        for specimen_work_order in lib["specimen_work_order"]:
            work_order = specimen_work_order["work_order"]
            rows.append({column: quote_text(get_field(work_order, field, column)) for column, field in columns.items()})

    df = pd.DataFrame(rows, columns=list(columns)).drop_duplicates()
    df["Sample Count"] = len(get_project_libraries(runinfo, project, include_controls=True))
    df.to_csv(data_dir / "work_order_info.tsv", index=False, **TSV_ARGS)


# Columns of all_tophits.csv, in order
MEGABLAST_COLUMNS = [
    "ID",
    "Specimen Name",
    "Expected Species",
    "Contamination Status",
    "Total Reads",
    "Total Reads With Hits",
    "Total Reads With Hits Percent",
    "Genus Match Hits",
    "Genus Match Hits Percent",
    "Other Genus Hits",
    "Other Genus Hits Percent",
    "Species Match Hits",
    "Species Match Hits Percent",
    "Other Species Hits",
    "Other Species Hits Percent",
]

# topHits.txt lines shaped "Label:<tab><count>(<percent>%)", and the all_tophits.csv column they fill.
# The percent goes in the column of the same name plus " Percent".
TOPHITS_COUNT_LINES = {
    "Total reads with hits": "Total Reads With Hits",
    "Genus match hits": "Genus Match Hits",
    "Other genus hits": "Other Genus Hits",
    "Species match hits": "Species Match Hits",
    "Other species hits": "Other Species Hits",
}

# "QC Report (yellow) for:{run}:{library}:lane:{lane}:{project}:..." (PacBio has cell:{cell})
TOPHITS_HEADER = re.compile(r"^QC Report \([^)]*\) for:[^:]+:(?P<library>[^:]+):(?:lane|cell):(?P<unit>[^:]+):")
# Top species table row: read count, % reads, % hit reads, species
TOPHITS_TABLE_ROW = re.compile(r"^(\d+)\t([\d.]+)\t([\d.]+)\t(.+)$")
# PacBio file name: barcoded cells have the barcode after "hifi_reads", e.g.
# m84050_..._s1.hifi_reads.bc2091.subset.fasta.blast.topHits.txt; unbarcoded cells have none, e.g.
# m84050_..._s1.hifi_reads.subset.fasta.blast.topHits.txt
PACBIO_BARCODE = re.compile(r"\.hifi_reads\.(?:(?P<barcode>[^.]+)\.)?subset\.")


def tophits_format_error(tophits_file: Path, problem: str) -> ValueError:
    return ValueError(f"Unexpected topHits format in {tophits_file}: {problem}")


def find_tophits_files(run_path: Path, platform: Platform, project: str, folder: str) -> list[Path]:
    """
    Return the sorted topHits.txt paths for one folder, or an empty list if it has no MegaBlast QC output
    """
    return sorted(
        (
            tophits_file
            for qc_dir in platform.tophits_dirs(run_path, project, folder)
            if qc_dir.is_dir()
            for tophits_file in qc_dir.iterdir()
            if tophits_file.name.endswith("topHits.txt")
        ),
        key=str,
    )


def write_megablast_qc(
    run_path: Path,
    runinfo: Runinfo,
    platform: Platform,
    data_dir: Path,
    project: str,
    folder: str,
    tophits_count: int = 0,
) -> int:
    """
    Writes the MegaBlast QC tables for one folder: all_tophits.csv with one row per topHits file,
    and {specimen}_tophits_{n}.csv with each file's top species. tophits_count is the number of
    topHits files already written by earlier calls; the updated count is returned so numbering
    stays unique across folders. Row n of the combined all_tophits.csv tables matches
    {specimen}_tophits_{n}.csv.
    """
    tophits_files = find_tophits_files(run_path, platform, project, folder)
    if not tophits_files:
        log.warning("No MegaBlast topHits files found for %s; MegaBlast QC is left out for it", folder)
        return tophits_count

    start_count = tophits_count
    rows = []
    for tophits_file in tophits_files:
        # Each file is fully parsed before anything is written, so a bad file can't leave
        # a row without its top species table (which would shift every later table's number)
        row, top_species = parse_tophits_file(tophits_file, runinfo, platform)
        # Specimen names come from the LIMS, so keep only filename-safe characters
        safe_name = re.sub(r"[^\w.-]", "_", str(row["Specimen Name"]))
        top_species.to_csv(data_dir / f"{safe_name}_tophits_{tophits_count}.csv", index=True)
        rows.append(row)
        tophits_count += 1

    # Row numbers continue from earlier folders: MultiQC merges every all_tophits.csv into one
    # table keyed on this column, so repeats would overwrite rows.
    df = pd.DataFrame(rows, columns=MEGABLAST_COLUMNS)
    df.index = pd.RangeIndex(start_count + 1, start_count + 1 + len(df))
    df.to_csv(data_dir / "all_tophits.csv", index=True, index_label="No.")
    return tophits_count


def parse_tophits_file(
    tophits_file: Path, runinfo: Runinfo, platform: Platform
) -> tuple[dict[str, FieldValue], pd.DataFrame]:
    """
    Parse one MegaBlast topHits.txt. Returns its all_tophits.csv row and its top species table.
    Specimen Name and Expected Species come from the library's runinfo.json entry.
    """
    lines = tophits_file.read_text(encoding="utf-8").splitlines()

    header = TOPHITS_HEADER.match(lines[0]) if lines else None
    if header is None:
        raise tophits_format_error(
            tophits_file, "first line is not a 'QC Report (...) for:run:library:lane|cell:...' header"
        )
    library_id = header["library"]
    if library_id not in runinfo["libraries"]:
        raise KeyError(f"Library '{library_id}' from {tophits_file} is not in runinfo.json")
    library = runinfo["libraries"][library_id]

    def value(label: str) -> str:
        """Text after "label:<tab>" on the one line starting with it"""
        found = [line.split("\t", 1)[1].strip() for line in lines if line.startswith(f"{label}:\t")]
        if len(found) != 1:
            raise tophits_format_error(tophits_file, f"expected one '{label}:' line, found {len(found)}")
        return found[0]

    total_reads = value("Total_reads")
    if not total_reads.isdigit():
        raise tophits_format_error(tophits_file, f"'Total_reads:' is not a number: {total_reads!r}")
    row: dict[str, FieldValue] = {
        "ID": platform.tophits_row_id(tophits_file, library_id, header["unit"]),
        "Specimen Name": get_field(library, "sample.specimen_name", "Specimen Name"),
        "Expected Species": get_field(library, "species", "Expected Species"),
        "Contamination Status": value("Contamination status"),
        "Total Reads": total_reads,
    }
    for label, column in TOPHITS_COUNT_LINES.items():
        count = re.fullmatch(r"(\d+)\(([\d.]+%)\)", value(label))
        if count is None:
            raise tophits_format_error(tophits_file, f"'{label}:' is not '<count>(<percent>%)'")
        row[column] = count[1]
        row[f"{column} Percent"] = count[2]

    # Top species table: rows after the "Reads_count ..." header, up to a blank line or the end
    table_start = next((i for i, line in enumerate(lines) if line.startswith("Reads_count\t")), None)
    if table_start is None:
        raise tophits_format_error(tophits_file, "no 'Reads_count' header for the top species table")
    top_species = []
    for line in lines[table_start + 1 :]:
        if not line.strip():
            break
        match = TOPHITS_TABLE_ROW.match(line)
        if match is None:
            raise tophits_format_error(tophits_file, f"unexpected line in the top species table: {line!r}")
        top_species.append(
            {
                "Read Count": match[1],
                "Percent Reads": match[2],
                "Percent Hit Reads": match[3],
                "Species": match[4].replace(",", " ").strip(),
            }
        )
    columns = ["Read Count", "Percent Reads", "Percent Hit Reads", "Species"]
    return row, pd.DataFrame(top_species, columns=columns)


def require_project_column(df: pd.DataFrame, csv_file: Path, block: str | None = None) -> None:
    """
    Raise a clear error if a stats table has no Project column (usually a changed file layout)
    """
    if "Project" not in df.columns:
        where = f"the {block} block of {csv_file}" if block else str(csv_file)
        raise ValueError(f"No 'Project' column in {where}; has the file layout changed?")


def copy_assets(outdir: Path) -> None:
    """
    Copies the report theme files (CSS, logos) next to the MultiQC config
    """
    (outdir / "assets").mkdir(parents=True, exist_ok=True)
    for asset in ASSETS_DIR.iterdir():
        shutil.copy(asset, outdir / "assets")


def tidy_tables(data_dir: Path) -> None:
    """
    Tidies every .csv and .tsv table in the folder: drops duplicate rows and empty columns,
    and turns the pandas index column into a 1-based "No." column. Tables with no rows are
    removed, so MultiQC leaves their section out instead of showing an empty table.
    """
    for table in data_dir.iterdir():
        if table.suffix not in (".csv", ".tsv"):
            continue
        read_args = TSV_ARGS if table.suffix == ".tsv" else {}
        # Every cell is read as text so values are written back exactly as they were
        # (no "00123" -> 123, "50" -> "50.0" or "0.50" -> "0.5")
        df = pd.read_csv(table, dtype=str, keep_default_na=False, **read_args)
        if df.empty:
            log.info("%s has no rows; leaving it out of the report", table)
            table.unlink()
            continue

        df = df.loc[:, ~(df == "").all()]

        # The pandas index column numbers every row, so duplicates are found on the other columns
        # and the rows are renumbered afterwards to leave no gaps
        has_index = df.columns[0] == "Unnamed: 0"
        df = df.drop_duplicates(subset=list(df.columns[1:] if has_index else df.columns))
        if has_index:
            df = df.rename(columns={"Unnamed: 0": "No."})
            df["No."] = [str(n) for n in range(1, len(df) + 1)]

        df.to_csv(table, index=False, **read_args)


# ---------------------------------------------------------------------------
# Illumina
# ---------------------------------------------------------------------------


def find_unaligned_folders(run_path: Path, runinfo: Runinfo, project: str, requested: str | None) -> list[str]:
    """
    Return the Unaligned* folders that contain QA/Project_{project}, sorted by name.
    If requested is given, return only that folder.
    """
    if requested:
        project_qa = run_path / requested / "QA" / f"Project_{project}"
        if not project_qa.is_dir():
            raise NotADirectoryError(f"{project_qa} not found; is {requested} the right Unaligned folder?")
        return [requested]

    unaligned = sorted(p.name for p in run_path.iterdir() if p.is_dir() and p.name.startswith("Unaligned"))
    if not unaligned:
        raise FileNotFoundError(f"No Unaligned folder found in {run_path}")

    with_project = [u for u in unaligned if (run_path / u / "QA" / f"Project_{project}").is_dir()]
    skipped = sorted(set(unaligned) - set(with_project))
    if skipped:
        log.info("Skipping Unaligned folders without Project_%s: %s", project, ", ".join(skipped))
    if not with_project:
        raise FileNotFoundError(f"No Unaligned folder in {run_path} contains QA/Project_{project}")
    log.info("Using Unaligned folders: %s", ", ".join(with_project))
    return with_project


def illumina_module_search_dirs(run_path: Path, project: str, folder: str) -> list[Path]:
    """
    The project's QA folder holds FastQC and other tool output MultiQC reads directly
    """
    return [run_path / folder / "QA" / f"Project_{project}"]


def illumina_tophits_dirs(run_path: Path, project: str, folder: str) -> list[Path]:
    qc_dir = run_path / folder / "QC" / f"Project_{project}"
    if not qc_dir.is_dir():
        return []
    return [sample / "QC" for sample in qc_dir.iterdir() if sample.name.startswith("Sample")]


def illumina_tophits_row_id(tophits_file: Path, library_id: str, lane: str) -> str:
    return f"{library_id}-{lane}"


# PA_Logs/images plots, named {run}_{lane}_{suffix}: suffix -> (section id part, section title).
# Each lane's plots appear in the report in this order.
ILLUMINA_PLOTS = {
    "ClusterCount-by-lane.png": ("cluster_count_by_lane", "Cluster Count By Lane"),
    "Intensity-by-cycle_Intensity.png": ("intensity_by_cycle", "Intensity By Cycle"),
    "flowcell-Intensity.png": ("flowcell_intensity", "Flowcell Intensity"),
    "q-heat-map.png": ("q_heat_map", "Q Heat Map"),
}


def find_lane_images(images_in_dir: list[str], lane: str) -> dict[str, str]:
    """
    Return {suffix: image name} for one lane's ILLUMINA_PLOTS. Names must end in exactly
    "_{lane}_{suffix}", so another lane's plots or a file that happens to contain "_{lane}_"
    are not picked up.
    """
    return {
        suffix: image for suffix in ILLUMINA_PLOTS for image in images_in_dir if image.endswith(f"_{lane}_{suffix}")
    }


def illumina_plot_images(run_path: Path, runinfo: Runinfo, project: str, folders: list[str]) -> list[PlotImage]:
    """
    The run plots found for this project's lanes
    """
    image_dir = run_path / "PA_Logs" / "images"
    images_in_dir = [image.name for image in image_dir.iterdir()]
    plots = []
    for lane in get_lanes(runinfo, project):
        for suffix, image in find_lane_images(images_in_dir, lane).items():
            section_part, title = ILLUMINA_PLOTS[suffix]
            plots.append(PlotImage(f"lane_{lane}_{section_part}", f"{title} (Lane {lane})", image_dir / image, image))
    return plots


def write_illumina_read_stats(run_path: Path, data_dir: Path, project: str) -> None:
    """
    Writes this project's rows of PA_Logs/trimming.logs.csv, with readable column names
    """
    trimming_logs = run_path / "PA_Logs" / "trimming.logs.csv"
    # Read as text so numbers are written back as they appear in the file
    df = pd.read_csv(trimming_logs, dtype=str)
    require_project_column(df, trimming_logs)
    df = df[df["Project"] == project].reset_index(drop=True)
    if df.empty:
        log.warning("No read stats rows for project %s in %s", project, trimming_logs)
    df = df.rename(columns=TRIMMING_LOG_COLUMNS)
    df.to_csv(data_dir / "trimming.logs_clean.csv")


# Readable names for trimming.logs.csv columns
TRIMMING_LOG_COLUMNS = {
    f"{read}_{stat}": f"{read} {title}"
    for read in ("R1", "R2")
    for stat, title in {
        "Bases": "Bases",
        "Reads": "Reads",
        "RL": "RL",
        "TrimBases": "Trim Bases",
        "TrimReads": "Trim Reads",
        "TrimRL": "Trim RL",
    }.items()
}


def write_flowcell_metrics(run_path: Path, data_dir: Path) -> None:
    """
    Writes the run metrics table from PA_Logs/summary.csv: the block after the two title lines,
    up to the first blank line, with "nan" cells left empty. It's written as summary_clean.csv
    because MultiQC's Illumina module would claim a file named summary.csv.
    """
    summary_csv = run_path / "PA_Logs" / "summary.csv"
    lines = summary_csv.read_text(encoding="utf-8").splitlines(True)
    blank_index = next((i for i, line in enumerate(lines) if not line.strip()), None)
    if blank_index is None:
        raise ValueError(f"No blank line ending the run metrics block in {summary_csv}")
    if blank_index <= 2:
        raise ValueError(f"No run metrics before the first blank line in {summary_csv}")
    with open(data_dir / "summary_clean.csv", "w", encoding="utf-8", newline="") as fout:
        writer = csv.writer(fout, lineterminator="\n")
        for row in csv.reader(lines[2:blank_index]):
            writer.writerow("" if cell == "nan" else cell for cell in row)


def copy_illumina_run_files(run_path: Path, runinfo: Runinfo, data_dir: Path, project: str, folders: list[str]) -> None:
    """
    Copies the run plots for this project's lanes and writes the flow cell metrics table
    """
    image_dir = run_path / "PA_Logs" / "images"
    images_in_dir = [image.name for image in image_dir.iterdir()]
    (data_dir / "images").mkdir(exist_ok=True)
    for lane in get_lanes(runinfo, project):
        lane_images = find_lane_images(images_in_dir, lane)
        missing = [suffix for suffix in ILLUMINA_PLOTS if suffix not in lane_images]
        if missing:
            log.warning("Lane %s plots not found in %s: %s", lane, image_dir, ", ".join(missing))
        for image in lane_images.values():
            shutil.copy(image_dir / image, data_dir / "images")
    write_flowcell_metrics(run_path, data_dir)


def copy_illumina_folder_files(run_path: Path, runinfo: Runinfo, data_dir: Path, project: str, folder: str) -> None:
    """
    Copies the tool output MultiQC modules read for one Unaligned folder
    """
    copy_trimming_logs(run_path, data_dir, folder)
    copy_bcl2fastq_stats(run_path, data_dir, folder)


def copy_trimming_logs(run_path: Path, data_dir: Path, folder: str) -> None:
    """
    Copies the Trimmomatic and cutadapt logs for one Unaligned folder. cutadapt logs are
    renamed *.json so MultiQC's cutadapt module reads them.
    """
    log_dir = run_path / folder / "TMP" / "adaptor_trimming_logs"
    if not log_dir.is_dir():
        log.warning("%s not found; Trimmomatic and cutadapt are left out for %s", log_dir, folder)
        return
    copied = 0
    for trimming_log in sorted(log_dir.iterdir()):
        if trimming_log.name.endswith(".adap_trim.err"):
            shutil.copy(trimming_log, data_dir)
            copied += 1
        elif trimming_log.name.endswith(".cutadapt.log"):
            shutil.copy(trimming_log, data_dir / f"{trimming_log.name}.json")
            copied += 1
    if not copied:
        log.warning("No Trimmomatic (*.adap_trim.err) or cutadapt (*.cutadapt.log) logs in %s", log_dir)


def copy_bcl2fastq_stats(run_path: Path, data_dir: Path, folder: str) -> None:
    """
    Copies one Unaligned folder's bcl2fastq Stats.json. MultiQC's bcl2fastq module keys lanes by
    RunId, so two Unaligned folders covering the same lane would overwrite each other; RunId is
    tagged with the folder name.
    """
    stats_json = run_path / folder / "Stats" / "Stats.json"
    if not stats_json.is_file():
        log.warning("%s not found; bcl2fastq is left out for %s", stats_json, folder)
        return
    try:
        stats = json.loads(stats_json.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise ValueError(f"Could not parse {stats_json}: {e}") from None
    if "RunId" not in stats:
        raise ValueError(f"No RunId in {stats_json}")
    stats["RunId"] = f"{stats['RunId']} ({folder})"
    with open(data_dir / "Stats.json", "w", encoding="utf-8") as fout:
        json.dump(stats, fout, indent=2)


# ---------------------------------------------------------------------------
# PacBio
# ---------------------------------------------------------------------------


def find_cells(run_path: Path, runinfo: Runinfo, project: str, requested: str | None) -> list[str]:
    """
    Return the PacBio cell folders listed for the project in runinfo.json.
    If requested is given, return only that cell.
    """
    cells = get_cells(runinfo, project)
    if not cells:
        raise ValueError(f"No cells listed for project {project} in {run_path}/PA_Logs/runinfo.json")
    if requested:
        if requested not in cells:
            raise ValueError(f"Cell {requested} does not hold project {project} (cells: {', '.join(cells)})")
        cells = [requested]
    missing = [c for c in cells if not (run_path / c).is_dir()]
    if missing:
        raise NotADirectoryError(f"Cell folders listed in runinfo.json not found in {run_path}: {', '.join(missing)}")
    log.info("Using cells: %s", ", ".join(cells))
    return cells


def pacbio_tophits_dirs(run_path: Path, project: str, cell: str) -> list[Path]:
    return [run_path / cell / "hifi_reads" / "QC"]


def pacbio_tophits_row_id(tophits_file: Path, library_id: str, cell: str) -> str:
    barcode = PACBIO_BARCODE.search(tophits_file.name)
    if barcode is None:
        raise tophits_format_error(
            tophits_file, "file name has neither '.hifi_reads.<barcode>.subset.' nor '.hifi_reads.subset.'"
        )
    # Unbarcoded cells hold one library, so the cell identifies the row
    return barcode["barcode"] or cell


# Each cell's Filter/dataset-reports/reports plots: file name -> (section id part, section title).
# Each cell's plots appear in the report in this order.
PACBIO_PLOTS = {
    "hexbin_length_plot.png": ("hexbin_length_plot", "Hexbin Length Plot"),
    "base_yield_plot.png": ("base_yield_plot", "Base Yield Plot"),
    "ccs_readlength_qv_hist2d.hexbin.png": ("ccs_readlength_qv_hist2d", "CCS Read Length QV Histogram"),
    "readLenDist0.png": ("readLenDist0", "Read Length Distribution"),
    "ccs_all_readlength_hist_plot.png": ("ccs_all_readlength_hist_plot", "CCS All Read Length Histogram"),
    "ccs_npasses_hist.png": ("ccs_npasses_hist", "CCS Number of Passes Histogram"),
    "ccs_accuracy_hist.png": ("ccs_accuracy_hist", "CCS Accuracy Histogram"),
    "ccs_combined_readlength_hist_plot.png": (
        "ccs_combined_readlength_hist_plot",
        "CCS Combined Read Length Histogram",
    ),
    "readlength_plot.png": ("readlength_plot", "Read Length Plot"),
    "concordance_plot.png": ("concordance_plot", "Concordance Plot"),
    "5mc_detections_hist.png": ("5mc_detections_hist", "5mC Detections Histogram"),
}


def pacbio_plot_images(run_path: Path, runinfo: Runinfo, project: str, cells: list[str]) -> list[PlotImage]:
    """
    Every expected plot for each cell. Copies are prefixed with the cell name so MultiQC can
    tell cells apart.
    """
    return [
        PlotImage(
            f"{section_part}_{cell}",
            f"{title} ({cell})",
            run_path / cell / "Filter" / "dataset-reports" / "reports" / image,
            f"{cell}_{image}",
        )
        for cell in cells
        for image, (section_part, title) in PACBIO_PLOTS.items()
    ]


def write_pacbio_read_stats(run_path: Path, data_dir: Path, project: str) -> None:
    """
    Writes the polymerase and HiFi read stats tables from PA_Logs/run_stats.csv, which has a
    polymerase block (header on line 3) followed by a "HiFi Stats:" block with its own header
    """
    run_stats = run_path / "PA_Logs" / "run_stats.csv"
    lines = run_stats.read_text(encoding="utf-8").splitlines()
    split_index = next((i for i, line in enumerate(lines) if "HiFi Stats:" in line), None)
    if split_index is None:
        raise ValueError(f"'HiFi Stats:' line not found in {run_stats}")

    # Read as text so numbers are written back as they appear in the file
    pf_df = pd.read_csv(run_stats, header=2, skiprows=range(split_index, len(lines)), dtype=str)
    hifi_df = pd.read_csv(run_stats, header=split_index + 1, dtype=str)
    require_project_column(pf_df, run_stats, "polymerase read stats")
    require_project_column(hifi_df, run_stats, "HiFi Stats")

    for df, name, out_file in (
        (pf_df, "polymerase read stats", "pf_polymerase_read_stats.csv"),
        (hifi_df, "HiFi stats", "hifi_stats.csv"),
    ):
        df = df[df["Project"] == project]
        if df.empty:
            log.warning("No %s rows for project %s in %s", name, project, run_stats)
        df.to_csv(data_dir / out_file)


def copy_pacbio_run_files(run_path: Path, runinfo: Runinfo, data_dir: Path, project: str, cells: list[str]) -> None:
    """
    PacBio has no run-wide files beyond the read stats; plots are per cell
    """


def copy_pacbio_folder_files(run_path: Path, runinfo: Runinfo, data_dir: Path, project: str, cell: str) -> None:
    """
    Copies one cell's plots and Iso-Seq reports
    """
    copy_pacbio_plots(run_path, runinfo, data_dir, project, cell)
    copy_isoseq_reports(run_path, data_dir, cell)


def copy_pacbio_plots(run_path: Path, runinfo: Runinfo, data_dir: Path, project: str, cell: str) -> None:
    """
    Copies one cell's plots, prefixed with the cell name so MultiQC can tell cells apart
    """
    image_dir = run_path / cell / "Filter" / "dataset-reports" / "reports"
    if not image_dir.is_dir():
        log.warning("Plot folder %s not found; no sequencing plots for cell %s", image_dir, cell)
        return
    (data_dir / "images").mkdir(exist_ok=True)
    missing = []
    for plot in pacbio_plot_images(run_path, runinfo, project, [cell]):
        if plot.source.is_file():
            shutil.copy(plot.source, data_dir / "images" / plot.file_name)
        else:
            missing.append(plot.source.name)
    if missing:
        log.warning("Cell %s plots not found in %s: %s", cell, image_dir, ", ".join(missing))


def copy_isoseq_reports(run_path: Path, data_dir: Path, cell: str) -> None:
    """
    Copies one cell's Iso-Seq reports, if Iso-Seq was run on it
    """
    report_dir = run_path / cell / "hifi_reads" / "read_seg_iso_seq" / "cromwell_out" / "outputs"
    if not report_dir.is_dir():
        log.info("No Iso-Seq outputs for cell %s", cell)
        return
    (data_dir / "isoseq").mkdir(exist_ok=True)
    for report in report_dir.iterdir():
        if report.suffix in (".csv", ".json"):
            shutil.copy(report, data_dir / "isoseq")


PLATFORMS = {
    "illumina": Platform(
        display_name="Illumina",
        folder_option="unaligned",
        find_folders=find_unaligned_folders,
        get_units=get_lanes,
        units_key="lanes",
        unit_field="_lane",
        module_search_dirs=illumina_module_search_dirs,
        tophits_dirs=illumina_tophits_dirs,
        tophits_row_id=illumina_tophits_row_id,
        plot_images=illumina_plot_images,
        write_read_stats=write_illumina_read_stats,
        copy_run_files=copy_illumina_run_files,
        copy_folder_files=copy_illumina_folder_files,
        software_versions={"Trimmomatic": "0.33", "bcl2fastq": "2.20.0.422"},
    ),
    "pacbio": Platform(
        display_name="PacBio",
        folder_option="cell",
        find_folders=find_cells,
        get_units=get_cells,
        units_key="cells",
        unit_field="_cell",
        module_search_dirs=lambda run_path, project, cell: [],
        tophits_dirs=pacbio_tophits_dirs,
        tophits_row_id=pacbio_tophits_row_id,
        plot_images=pacbio_plot_images,
        write_read_stats=write_pacbio_read_stats,
        copy_run_files=copy_pacbio_run_files,
        copy_folder_files=copy_pacbio_folder_files,
        software_versions={},
    ),
}
