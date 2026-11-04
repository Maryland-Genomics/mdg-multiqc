import argparse
import logging
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
    parser.add_argument("-nr", "--noregen", 
                        action="store_true", 
                        help="Don't regenerate parsed files and MultiQC config"
                        )
    parser.add_argument("-nm", "--nomultiqc", 
                        action="store_true", 
                        help="Don't run MultiQC"
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


def run_multiqc(search_dirs: list[Path], config_file: Path, outdir: Path, project: str) -> None:
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


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = get_args()
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

    if not args.noregen:
        outdir.mkdir(parents=True, exist_ok=True)
        MDG_Run_Parser.parse(str(run_path), project, sequencing_platform, folders, str(outdir))
        MultiQC_Config_Generator.generate_config(str(run_path), project, sequencing_platform, folders, str(outdir))

    if not args.nomultiqc:
        run_multiqc(search_dirs, outdir / f"multiqc_config_{project}.yaml", outdir, project)


if __name__ == "__main__":
    main()
