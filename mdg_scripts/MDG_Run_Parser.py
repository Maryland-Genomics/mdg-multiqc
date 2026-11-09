import csv
import json
import logging
import pandas as pd
import os
import shutil
import re
import yaml

log = logging.getLogger(__name__)


def parse(p_path : str, p_project : str, p_sequencing_platform : str, p_folders : list[str], p_outdir : str):
    """
    Generate the parsed data files MultiQC reads.

    p_folders are the Unaligned folders (Illumina) or cell folders (PacBio) to report on.
    Files covering the whole run go to {outdir}/run_data, and files specific to one
    Unaligned folder or cell go to {outdir}/{folder}_data.
    """
    project = p_project
    path = p_path
    sequencing_platform = p_sequencing_platform
    outdir = p_outdir
    runinfo = load_runinfo(path)
    column_config = load_column_config()

    copy_assets(outdir)

    run_data_dir = get_run_data_dir(outdir)
    reset_dir(run_data_dir)
    parse_work_order_info(runinfo, project, run_data_dir, column_config["work_order"])
    generate_sample_info(
        runinfo, run_data_dir, project, sequencing_platform, column_config["sample_info"][sequencing_platform]
    )
    generate_sequencing_info(
        runinfo, run_data_dir, sequencing_platform, column_config["sequencing_info"][sequencing_platform]
    )
    generate_read_stats(path, run_data_dir, project, sequencing_platform)
    if sequencing_platform == "illumina":
        copy_illumina_plots(path, runinfo, run_data_dir, project)
        parse_flowcell_info(path, run_data_dir)
    post_process_csv(run_data_dir)

    # Tophits numbering continues across folders so report section ids stay unique
    tophits_count = 0
    for folder in p_folders:
        folder_data_dir = get_data_dir(outdir, folder)
        reset_dir(folder_data_dir)
        if sequencing_platform == "pacbio":
            copy_pacBio_plots(path, folder_data_dir, folder)
        tophits_count = parse_megablastqc_info(
            path, runinfo, folder_data_dir, project, sequencing_platform, folder, tophits_count
        )
        move_module_files(path, folder_data_dir, folder, sequencing_platform)
        post_process_csv(folder_data_dir)


def get_run_data_dir(outdir):
    return f"{outdir}/run_data"


def get_data_dir(outdir, folder):
    """
    Folder for data specific to one Unaligned folder (Illumina) or cell (PacBio)
    """
    return f"{outdir}/{folder}_data"


RUNINFO_COLUMNS_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)), "runinfo_columns.yaml")

# Tables built from runinfo.json are tab-separated: MultiQC's custom content parser doesn't
# understand CSV quoting, so commas inside LIMS values (contact lists, comments) would split columns.
# get_field() strips tabs and line breaks from values, so no quoting is needed.
TSV_ARGS = {"sep": "\t", "quoting": csv.QUOTE_NONE}


def load_runinfo(path):
    """
    Load the run's PA_Logs/runinfo.json
    """
    runinfo_json = f"{path}/PA_Logs/runinfo.json"
    try:
        with open(runinfo_json, encoding="utf-8") as fin:
            return json.load(fin)
    except FileNotFoundError:
        raise FileNotFoundError(f"Run info file not found: {runinfo_json}") from None
    except json.JSONDecodeError as e:
        raise ValueError(f"Could not parse {runinfo_json}: {e}") from None


def load_column_config():
    """
    Load the report table column definitions from runinfo_columns.yaml
    """
    with open(RUNINFO_COLUMNS_FILE, encoding="utf-8") as fin:
        return yaml.safe_load(fin)


def get_lanes(runinfo, project):
    """
    Return the Illumina lanes (as strings, in lane order) holding this project
    """
    lanes = [lane for lane, info in runinfo["lanes"].items() if project in info["projects"]["value"]]
    return sorted(lanes, key=int)


def get_cells(runinfo, project):
    """
    Return the PacBio cell folders (e.g. 1_A01) holding this project
    """
    return [cell for cell, info in runinfo["cells"].items() if project in info["projects"]["value"]]


def get_project_libraries(runinfo, project, include_controls=False):
    """
    Return {library_id: library} for this project's libraries. With include_controls, also
    return control libraries whose control_projects include this project.
    """
    return {
        lib_id: lib
        for lib_id, lib in runinfo["libraries"].items()
        if lib["project"]["value"] == project
        or (include_controls and project in lib["control_projects"]["value"])
    }


def get_field(entry, field, column, context=None):
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

    def resolve(node, keys):
        # Lists (e.g. specimen_work_order): collect the value from every item
        if isinstance(node, list):
            values = list(dict.fromkeys(str(v) for v in (resolve(item, keys) for item in node) if v != ""))
            return "; ".join(values)
        if not keys:
            if not (isinstance(node, dict) and "value" in node):
                raise ValueError(f"Column '{column}': runinfo.json field '{field}' is not a single value")
            value = node["value"]
            if value is None:
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


def reset_dir(directory):
    """
    Remove and recreate a generated data folder so files from a previous run can't linger
    """
    shutil.rmtree(directory, ignore_errors=True)
    os.makedirs(directory)


PLATFORM_NAMES = {"illumina": "Illumina", "pacbio": "PacBio"}


def generate_sequencing_info(runinfo, mdg_parsed_dir, sequencing_platform, columns):
    """
    Writes the Sequencing Info table: one row describing the run, with the columns defined in
    runinfo_columns.yaml (paths relative to runinfo.json's "run" entry)
    """
    context = {"_platform": PLATFORM_NAMES[sequencing_platform]}
    row = {column: get_field(runinfo["run"], field, column, context) for column, field in columns.items()}
    # Every value here is an identifier, so quote it: MultiQC leaves quoted cells as text
    # instead of turning all-digit values like kit lots into numbers (dropping leading zeros
    # and adding thousands separators). Empty values stay unquoted so empty columns are still dropped.
    row = {column: f'"{value}"' if value != "" else "" for column, value in row.items()}
    df = pd.DataFrame([row], columns=list(columns))
    df.to_csv(f"{mdg_parsed_dir}/sequencing_info.tsv", index=False, **TSV_ARGS)


def parse_flowcell_info(path, mdg_parsed_dir):
    """
    Writes the run metrics table from PA_Logs/summary.csv: the block after the two title lines,
    up to the first blank line. It's written as summary_clean.csv because MultiQC's Illumina
    module would claim a file named summary.csv.
    """
    summary_csv = f"{path}/PA_Logs/summary.csv"
    with open(summary_csv) as fin:
        lines = fin.read().splitlines(True)
    blank_index = next((i for i, line in enumerate(lines) if not line.strip()), None)
    if blank_index is None:
        raise ValueError(f"No blank line ending the run metrics block in {summary_csv}")
    if blank_index <= 2:
        raise ValueError(f"No run metrics before the first blank line in {summary_csv}")
    with open(f"{mdg_parsed_dir}/summary_clean.csv", "w") as fout:
        for line in lines[2:blank_index]:
            fout.write(line.replace("nan", ""))


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
# PacBio barcode in the file name, e.g. m84050_..._s1.hifi_reads.bc2091.subset.fasta.blast.topHits.txt
PACBIO_BARCODE = re.compile(r"\.hifi_reads\.([^.]+)\.subset\.")


def find_tophits_files(path, project, sequencing_platform, folder):
    """
    Return the sorted topHits.txt paths for one Unaligned folder (Illumina) or cell (PacBio),
    or an empty list if the folder has no MegaBlast QC output
    """
    if sequencing_platform == "illumina":
        qc_dir = f"{path}/{folder}/QC/Project_{project}"
        if not os.path.isdir(qc_dir):
            return []
        sample_dirs = [f"{qc_dir}/{sample}/QC" for sample in os.listdir(qc_dir) if sample.startswith("Sample")]
    else:
        sample_dirs = [f"{path}/{folder}/hifi_reads/QC"]
    return sorted(
        f"{sample_dir}/{name}"
        for sample_dir in sample_dirs
        if os.path.isdir(sample_dir)
        for name in os.listdir(sample_dir)
        if name.endswith("topHits.txt")
    )


def parse_megablastqc_info(path, runinfo, mdg_parsed_dir, project, sequencing_platform, folder, tophits_count=0):
    """
    Writes the MegaBlast QC tables for one Unaligned folder (Illumina) or cell (PacBio):
    all_tophits.csv with one row per topHits file, and {specimen}_tophits_{n}.csv with each
    file's top species. tophits_count is the number of topHits files already written by
    earlier calls; the updated count is returned so numbering stays unique across folders.
    Row n of the combined all_tophits.csv tables matches {specimen}_tophits_{n}.csv.
    """
    tophits_files = find_tophits_files(path, project, sequencing_platform, folder)
    if not tophits_files:
        log.warning("No MegaBlast topHits files found for %s; MegaBlast QC is left out for it", folder)
        return tophits_count

    start_count = tophits_count
    rows = []
    for tophits_file in tophits_files:
        # Each file is fully parsed before anything is written, so a bad file can't leave
        # a row without its top species table (which would shift every later table's number)
        row, top_species = parse_tophits_file(tophits_file, runinfo, sequencing_platform)
        # Specimen names come from the LIMS, so keep only filename-safe characters
        safe_name = re.sub(r"[^\w.-]", "_", row["Specimen Name"])
        top_species.to_csv(f"{mdg_parsed_dir}/{safe_name}_tophits_{tophits_count}.csv", index=True)
        rows.append(row)
        tophits_count += 1

    # Row numbers continue from earlier folders: MultiQC merges every all_tophits.csv into one
    # table keyed on this column, so repeats would overwrite rows.
    df = pd.DataFrame(rows, columns=MEGABLAST_COLUMNS)
    df.index = pd.RangeIndex(start_count + 1, start_count + 1 + len(df))
    df.to_csv(f"{mdg_parsed_dir}/all_tophits.csv", index=True, index_label="No.")
    return tophits_count


def parse_tophits_file(tophits_file, runinfo, sequencing_platform):
    """
    Parse one MegaBlast topHits.txt. Returns its all_tophits.csv row and its top species table.
    Specimen Name and Expected Species come from the library's runinfo.json entry.
    """
    with open(tophits_file) as fin:
        lines = fin.read().splitlines()

    def format_error(problem):
        return ValueError(f"Unexpected topHits format in {tophits_file}: {problem}")

    header = TOPHITS_HEADER.match(lines[0]) if lines else None
    if header is None:
        raise format_error("first line is not a 'QC Report (...) for:run:library:lane|cell:...' header")
    library_id = header["library"]
    if library_id not in runinfo["libraries"]:
        raise KeyError(f"Library '{library_id}' from {tophits_file} is not in runinfo.json")
    library = runinfo["libraries"][library_id]

    if sequencing_platform == "illumina":
        row_id = f"{library_id}-{header['unit']}"
    else:
        barcode = PACBIO_BARCODE.search(os.path.basename(tophits_file))
        if barcode is None:
            raise format_error("no barcode ('.hifi_reads.<barcode>.subset.') in the file name")
        row_id = barcode[1]

    def value(label):
        """Text after "label:<tab>" on the one line starting with it"""
        found = [line.split("\t", 1)[1].strip() for line in lines if line.startswith(f"{label}:\t")]
        if len(found) != 1:
            raise format_error(f"expected one '{label}:' line, found {len(found)}")
        return found[0]

    total_reads = value("Total_reads")
    if not total_reads.isdigit():
        raise format_error(f"'Total_reads:' is not a number: {total_reads!r}")
    row = {
        "ID": row_id,
        "Specimen Name": get_field(library, "sample.specimen_name", "Specimen Name"),
        "Expected Species": get_field(library, "species", "Expected Species"),
        "Contamination Status": value("Contamination status"),
        "Total Reads": total_reads,
    }
    for label, column in TOPHITS_COUNT_LINES.items():
        count = re.fullmatch(r"(\d+)\(([\d.]+%)\)", value(label))
        if count is None:
            raise format_error(f"'{label}:' is not '<count>(<percent>%)'")
        row[column] = count[1]
        row[f"{column} Percent"] = count[2]

    # Top species table: rows after the "Reads_count ..." header, up to a blank line or the end
    table_start = next((i for i, line in enumerate(lines) if line.startswith("Reads_count\t")), None)
    if table_start is None:
        raise format_error("no 'Reads_count' header for the top species table")
    top_species = []
    for line in lines[table_start + 1:]:
        if not line.strip():
            break
        match = TOPHITS_TABLE_ROW.match(line)
        if match is None:
            raise format_error(f"unexpected line in the top species table: {line!r}")
        top_species.append({
            "Read Count": match[1],
            "Percent Reads": match[2],
            "Percent Hit Reads": match[3],
            "Species": match[4].replace(",", " ").strip(),
        })
    columns = ["Read Count", "Percent Reads", "Percent Hit Reads", "Species"]
    return row, pd.DataFrame(top_species, columns=columns)


def generate_sample_info(runinfo, mdg_parsed_dir, project, sequencing_platform, columns):
    """
    Writes the Sample Information table: one row per project library per lane (Illumina) or
    cell (PacBio), with the columns defined in runinfo_columns.yaml. Control libraries whose
    control_projects include this project are listed too.
    """
    if sequencing_platform == "illumina":
        placements = [("_lane", lane, runinfo["lanes"][lane]) for lane in get_lanes(runinfo, project)]
    elif sequencing_platform == "pacbio":
        placements = [("_cell", cell, runinfo["cells"][cell]) for cell in get_cells(runinfo, project)]
    else:
        raise ValueError(f"Unsupported sequencing platform '{sequencing_platform}'")

    rows = []
    for context_key, location, info in placements:
        for lib_id in info["libraries"]["value"]:
            lib = runinfo["libraries"][lib_id]
            if lib["project"]["value"] != project and project not in lib["control_projects"]["value"]:
                continue
            context = {context_key: location, "_library_id": lib_id}
            rows.append({column: get_field(lib, field, column, context) for column, field in columns.items()})

    df = pd.DataFrame(rows, columns=list(columns))
    df.to_csv(f"{mdg_parsed_dir}/runinfo_clean.tsv", **TSV_ARGS)


def generate_read_stats(path, mdg_parsed_dir, project, sequencing_platform):
    """
    Writes this project's read stats tables
    """
    if sequencing_platform == "illumina":
        generate_illumina_read_stats(path, mdg_parsed_dir, project)
    elif sequencing_platform == "pacbio":
        generate_pacBio_read_stats(path, mdg_parsed_dir, project)


def require_project_column(df, csv_file, block=None):
    """
    Raise a clear error if a stats table has no Project column (usually a changed file layout)
    """
    if "Project" not in df.columns:
        where = f"the {block} block of {csv_file}" if block else csv_file
        raise ValueError(f"No 'Project' column in {where}; has the file layout changed?")


def generate_pacBio_read_stats(path, mdg_parsed_dir, project):
    """
    Writes the polymerase and HiFi read stats tables from PA_Logs/run_stats.csv, which has a
    polymerase block (header on line 3) followed by a "HiFi Stats:" block with its own header
    """
    run_stats = f"{path}/PA_Logs/run_stats.csv"
    with open(run_stats) as fin:
        line_count = 0
        split_index = None
        for i, line in enumerate(fin):
            line_count = i + 1
            if "HiFi Stats:" in line:
                split_index = i
    if split_index is None:
        raise ValueError(f"'HiFi Stats:' line not found in {run_stats}")

    pf_df = pd.read_csv(run_stats, header=2, skiprows=range(split_index, line_count))
    hifi_df = pd.read_csv(run_stats, header=split_index + 1)
    require_project_column(pf_df, run_stats, "polymerase read stats")
    require_project_column(hifi_df, run_stats, "HiFi Stats")

    for df, name, out_file in (
        (pf_df, "polymerase read stats", "pf_polymerase_read_stats.csv"),
        (hifi_df, "HiFi stats", "hifi_stats.csv"),
    ):
        df = df[df["Project"] == project]
        if df.empty:
            log.warning("No %s rows for project %s in %s", name, project, run_stats)
        df.to_csv(f"{mdg_parsed_dir}/{out_file}")


def generate_illumina_read_stats(path, mdg_parsed_dir, project):
    """
    Writes this project's rows of PA_Logs/trimming.logs.csv, with readable column names
    """
    trimming_logs = f"{path}/PA_Logs/trimming.logs.csv"
    df = pd.read_csv(trimming_logs)
    require_project_column(df, trimming_logs)
    df = df[df["Project"] == project].reset_index(drop=True)
    if df.empty:
        log.warning("No read stats rows for project %s in %s", project, trimming_logs)
    df = df.rename(columns={"R1_Bases": "R1 Bases", "R1_Reads": "R1 Reads", "R1_RL": "R1 RL", "R1_TrimBases": "R1 Trim Bases", "R1_TrimReads": "R1 Trim Reads", "R1_TrimRL": "R1 Trim RL", "R2_Bases": "R2 Bases", "R2_Reads": "R2 Reads", "R2_RL": "R2 RL", "R2_TrimBases": "R2 Trim Bases",
                            "R2_TrimReads": "R2 Trim Reads", "R2_TrimRL": "R2 Trim RL"})
    df.to_csv(f"{mdg_parsed_dir}/trimming.logs_clean.csv")


PACBIO_PLOTS = [
    "hexbin_length_plot.png",
    "base_yield_plot.png",
    "ccs_readlength_qv_hist2d.hexbin.png",
    "readLenDist0.png",
    "ccs_all_readlength_hist_plot.png",
    "ccs_npasses_hist.png",
    "ccs_accuracy_hist.png",
    "ccs_combined_readlength_hist_plot.png",
    "readlength_plot.png",
    "concordance_plot.png",
    "5mc_detections_hist.png",
]


def copy_pacBio_plots(path, mdg_parsed_dir, cell):
    """
    Copies one cell's plots, prefixed with the cell name so MultiQC can tell cells apart
    """
    image_dir = f"{path}/{cell}/Filter/dataset-reports/reports"
    if not os.path.isdir(image_dir):
        log.warning("Plot folder %s not found; no sequencing plots for cell %s", image_dir, cell)
        return
    os.makedirs(f"{mdg_parsed_dir}/images", exist_ok=True)
    missing = []
    for image in PACBIO_PLOTS:
        source = f"{image_dir}/{image}"
        if os.path.isfile(source):
            shutil.copy(source, f"{mdg_parsed_dir}/images/{cell}_{image}")
        else:
            missing.append(image)
    if missing:
        log.warning("Cell %s plots not found in %s: %s", cell, image_dir, ", ".join(missing))


def copy_illumina_plots(path, runinfo, mdg_parsed_dir, project):
    """
    Copies the run plots for this project's lanes
    """
    image_dir = f"{path}/PA_Logs/images"
    images_in_dir = os.listdir(image_dir)
    os.makedirs(f"{mdg_parsed_dir}/images", exist_ok=True)
    for lane in get_lanes(runinfo, project):
        lane_images = [image for image in images_in_dir if f"_{lane}_" in image]
        if not lane_images:
            log.warning("No plots for lane %s in %s", lane, image_dir)
        for image in lane_images:
            shutil.copy(f"{image_dir}/{image}", f"{mdg_parsed_dir}/images")


def parse_work_order_info(runinfo, project, mdg_parsed_dir, columns):
    """
    Writes the Work Order table: one row per work order on this project's libraries, with the
    columns defined in runinfo_columns.yaml plus a Sample Count of the project's libraries and
    its control libraries
    """
    libraries = get_project_libraries(runinfo, project)
    rows = []
    for lib in libraries.values():
        for specimen_work_order in lib["specimen_work_order"]:
            work_order = specimen_work_order["work_order"]
            rows.append({column: get_field(work_order, field, column) for column, field in columns.items()})

    df = pd.DataFrame(rows, columns=list(columns)).drop_duplicates()
    df["Sample Count"] = len(get_project_libraries(runinfo, project, include_controls=True))
    df.to_csv(f"{mdg_parsed_dir}/work_order_info.tsv", index=False, **TSV_ARGS)


def move_module_files(path, mdg_parsed_dir, folder, sequencing_platform):
    """
    Copies tool output MultiQC modules read, for one Unaligned folder (Illumina) or cell (PacBio)
    """
    if sequencing_platform == "illumina":
        copy_trimming_logs(path, mdg_parsed_dir, folder)
        copy_bcl2fastq_stats(path, mdg_parsed_dir, folder)
    elif sequencing_platform == "pacbio":
        copy_isoseq_reports(path, mdg_parsed_dir, folder)


def copy_trimming_logs(path, mdg_parsed_dir, folder):
    """
    Copies the Trimmomatic and cutadapt logs for one Unaligned folder. cutadapt logs are
    renamed *.json so MultiQC's cutadapt module reads them.
    """
    log_dir = f"{path}/{folder}/TMP/adaptor_trimming_logs"
    if not os.path.isdir(log_dir):
        log.warning("%s not found; Trimmomatic and cutadapt are left out for %s", log_dir, folder)
        return
    copied = 0
    for name in sorted(os.listdir(log_dir)):
        if name.endswith(".adap_trim.err"):
            shutil.copy(f"{log_dir}/{name}", mdg_parsed_dir)
            copied += 1
        elif name.endswith(".cutadapt.log"):
            shutil.copy(f"{log_dir}/{name}", f"{mdg_parsed_dir}/{name}.json")
            copied += 1
    if not copied:
        log.warning("No Trimmomatic (*.adap_trim.err) or cutadapt (*.cutadapt.log) logs in %s", log_dir)


def copy_bcl2fastq_stats(path, mdg_parsed_dir, folder):
    """
    Copies one Unaligned folder's bcl2fastq Stats.json. MultiQC's bcl2fastq module keys lanes by
    RunId, so two Unaligned folders covering the same lane would overwrite each other; RunId is
    tagged with the folder name.
    """
    stats_json = f"{path}/{folder}/Stats/Stats.json"
    if not os.path.isfile(stats_json):
        log.warning("%s not found; bcl2fastq is left out for %s", stats_json, folder)
        return
    with open(stats_json) as fin:
        try:
            stats = json.load(fin)
        except json.JSONDecodeError as e:
            raise ValueError(f"Could not parse {stats_json}: {e}") from None
    if "RunId" not in stats:
        raise ValueError(f"No RunId in {stats_json}")
    stats["RunId"] = f"{stats['RunId']} ({folder})"
    with open(f"{mdg_parsed_dir}/Stats.json", "w") as fout:
        json.dump(stats, fout, indent=2)


def copy_isoseq_reports(path, mdg_parsed_dir, cell):
    """
    Copies one cell's Iso-Seq reports, if Iso-Seq was run on it
    """
    report_dir = f"{path}/{cell}/hifi_reads/read_seg_iso_seq/cromwell_out/outputs"
    if not os.path.isdir(report_dir):
        log.info("No Iso-Seq outputs for cell %s", cell)
        return
    os.makedirs(f"{mdg_parsed_dir}/isoseq", exist_ok=True)
    for report in os.listdir(report_dir):
        if report.endswith((".csv", ".json")):
            shutil.copy(f"{report_dir}/{report}", f"{mdg_parsed_dir}/isoseq")


def copy_assets(outdir):
    """
    Copies the report theme files (CSS, logos) next to the MultiQC config
    """
    assets_dir = f"{os.path.dirname(os.path.realpath(__file__))}/assets"
    os.makedirs(f"{outdir}/assets", exist_ok=True)
    for fname in os.listdir(assets_dir):
        shutil.copy(f"{assets_dir}/{fname}", f"{outdir}/assets")


def post_process_csv(mdg_parsed_dir):
    """
    Tidies every .csv and .tsv table in the folder: drops duplicate rows and empty columns,
    and turns the pandas index column into a 1-based "No." column. Tables with no rows are
    removed, so MultiQC leaves their section out instead of showing an empty table.
    """
    for x in os.listdir(mdg_parsed_dir):
        if x.endswith((".csv", ".tsv")):
            table = f"{mdg_parsed_dir}/{x}"
            read_args = TSV_ARGS if x.endswith(".tsv") else {}
            df = pd.read_csv(table, **read_args)
            if df.empty:
                log.info("%s has no rows; leaving it out of the report", table)
                os.remove(table)
                continue

            # drop all duplicates
            df = df.drop_duplicates()

            # drop all empty columns
            df = df.dropna(axis=1, how="all")
            df = df.fillna("")

            # Next, modify index column
            if df.columns[0] == "Unnamed: 0":
                df = df.rename(columns={"Unnamed: 0": "No."})
                df["No."] += 1 - df["No."].iloc[0]

            # Write results
            df.to_csv(table, index=False, **read_args)
