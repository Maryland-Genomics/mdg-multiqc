import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import MDG_Run_Parser
import MultiQC_Config_Generator

log = logging.getLogger("wrapper_multiqc")

SUPPORTED_PLATFORMS = ("illumina", "pacbio")


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MultiQC and generate MDG tables/images")
    parser.add_argument("-i", "--input", 
                        required=True, 
                        type=Path, 
                        help="Path to run folder"
                        )
    parser.add_argument("-p", "--project", 
                        required=True, 
                        help="Project which MultiQC will interpret"
                        )
    parser.add_argument("-u", "--unaligned",
                        help="Only use this Unaligned folder (Illumina only). By default every Unaligned folder containing the project is used",
                        )
    parser.add_argument("-c", "--cell",
                        help="Only use this cell folder, e.g. 1_A01 (PacBio only). By default every cell containing the project is used",
                        )
    parser.add_argument("-a", "--analysis",
                        type=Path,
                        help="File listing analysis files to include, one per line as SOURCE_PATH,DEST_NAME. "
                             "Files are copied to PA_Logs/mdg-reports/<PROJECT>/analysis_data",
                        )
    parser.add_argument("-nr", "--noregen",
                        action="store_true", 
                        help="Don't regenerate parsed files and MultiQC config"
                        )
    parser.add_argument("-nm", "--nomultiqc", 
                        action="store_true", 
                        help="Don't run MultiQC"
                        )
    parser.add_argument("--debug",
                        action="store_true",
                        help="Show the full traceback when the script fails"
                        )
    return parser.parse_args()


def check_project(runinfo: dict, project: str) -> None:
    """Raise if the project is not listed in the run's PA_Logs/runinfo.json."""
    projects = runinfo["run"]["projects"]["value"]
    log.info("Projects in run: %s", projects)
    if project not in projects:
        raise ValueError(f"{project} is not a valid project for this run folder (found: {', '.join(projects)})")


def find_unaligned_folders(run_path: Path, project: str, requested: Optional[str]) -> list[str]:
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


def find_cells(run_path: Path, runinfo: dict, project: str, requested: Optional[str]) -> list[str]:
    """
    Return the PacBio cell folders listed for the project in runinfo.json.
    If requested is given, return only that cell.
    """
    cells = MDG_Run_Parser.get_cells(runinfo, project)
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


def read_analysis_list(list_file: Path) -> list[tuple[Path, str]]:
    """
    Parse the analysis list file: one SOURCE_PATH,DEST_NAME per line.
    Blank lines and lines starting with # are ignored. Every line is validated
    before anything is copied, so a bad list never leaves a half-copied folder.
    """
    if not list_file.is_file():
        raise FileNotFoundError(f"Analysis list file not found: {list_file}")

    entries: list[tuple[Path, str]] = []
    seen: dict[str, int] = {}
    for line_no, raw in enumerate(list_file.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        where = f"{list_file} line {line_no}"
        # Split on the last comma so a comma inside the source path still works
        if "," not in line:
            raise ValueError(f"{where}: expected SOURCE_PATH,DEST_NAME but found no comma: {line}")
        source_str, dest_name = (part.strip() for part in line.rsplit(",", 1))
        if not source_str or not dest_name:
            raise ValueError(f"{where}: source path and destination name must both be set: {line}")
        if "/" in dest_name or dest_name in (".", ".."):
            raise ValueError(f"{where}: destination name must be a plain file name, not a path: {dest_name}")
        if dest_name in seen:
            raise ValueError(f"{where}: destination name {dest_name} already used on line {seen[dest_name]}")
        source = Path(source_str)
        if not source.is_file():
            raise FileNotFoundError(f"{where}: source file not found: {source}")
        seen[dest_name] = line_no
        entries.append((source, dest_name))

    if not entries:
        raise ValueError(f"Analysis list file {list_file} has no entries")
    return entries


def copy_analysis_files(entries: list[tuple[Path, str]], analysis_dir: Path) -> None:
    """Replace analysis_dir with fresh copies of the listed files, so it always matches the list."""
    if analysis_dir.exists():
        shutil.rmtree(analysis_dir)
    analysis_dir.mkdir(parents=True)
    for source, dest_name in entries:
        log.info("Copying %s -> %s", source, analysis_dir / dest_name)
        shutil.copy2(source, analysis_dir / dest_name)


def run_multiqc(search_dirs: list[Path], config_file: Path, outdir: Path, project: str) -> None:
    if shutil.which("multiqc") is None:
        raise FileNotFoundError("multiqc not found on PATH; activate the MultiQC environment first")
    if not config_file.is_file():
        raise FileNotFoundError(f"MultiQC config not found: {config_file} (run without --noregen to generate it)")
    missing = [str(d) for d in search_dirs if not d.is_dir()]
    if missing:
        raise FileNotFoundError(f"MultiQC input folders not found (run without --noregen?): {', '.join(missing)}")
    cmd = [
        "multiqc",
        *(str(d) for d in search_dirs),
        "--config",
        str(config_file),
        "--outdir",
        str(outdir),
        "-n",
        f"multiqc_{project}",
        "--flat",
        "--export",
        "--force",
    ]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main(args: argparse.Namespace) -> None:
    run_path: Path = args.input.resolve()
    project: str = args.project

    if not run_path.is_dir():
        raise NotADirectoryError(run_path)

    runinfo = MDG_Run_Parser.load_runinfo(str(run_path))
    check_project(runinfo, project)

    sequencing_platform = runinfo["run"]["platform"]["value"]
    if sequencing_platform not in SUPPORTED_PLATFORMS:
        raise ValueError(f"Unsupported sequencing platform '{sequencing_platform}' for {run_path}")

    outdir = run_path / "PA_Logs" / "mdg-reports" / project

    # Unaligned folders (Illumina) or cells (PacBio); each gets its own {outdir}/{folder}_data
    if sequencing_platform == "illumina":
        if args.cell:
            raise ValueError("--cell only applies to PacBio runs")
        folders = find_unaligned_folders(run_path, project, args.unaligned)
    else:
        if args.unaligned:
            raise ValueError("--unaligned only applies to Illumina runs")
        folders = find_cells(run_path, runinfo, project, args.cell)

    search_dirs = [Path(MDG_Run_Parser.get_run_data_dir(str(outdir)))]
    for folder in folders:
        if sequencing_platform == "illumina":
            search_dirs.append(run_path / folder / "QA" / f"Project_{project}")
        search_dirs.append(Path(MDG_Run_Parser.get_data_dir(str(outdir), folder)))

    # Analysis files are copied whenever --analysis is given (even with --noregen), and a
    # previously copied analysis_data folder is still included when --analysis is omitted
    analysis_dir = outdir / "analysis_data"
    if args.analysis:
        copy_analysis_files(read_analysis_list(args.analysis), analysis_dir)
    if analysis_dir.is_dir():
        search_dirs.append(analysis_dir)

    if not args.noregen:
        outdir.mkdir(parents=True, exist_ok=True)
        MDG_Run_Parser.parse(str(run_path), project, sequencing_platform, folders, str(outdir))
        MultiQC_Config_Generator.generate_config(str(run_path), project, sequencing_platform, folders, str(outdir))

    if not args.nomultiqc:
        run_multiqc(search_dirs, outdir / f"multiqc_config_{project}.yaml", outdir, project)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cli_args = get_args()
    try:
        main(cli_args)
    except Exception as e:
        # One readable line by default; --debug shows where it failed
        if cli_args.debug:
            log.exception("Report generation failed")
        else:
            log.error("Report generation failed: %s: %s (rerun with --debug for the full traceback)", type(e).__name__, e)
        sys.exit(1)
