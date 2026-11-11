import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import MDG_Run_Parser
import pandas as pd
import yaml
from MDG_Run_Parser import Platform, Runinfo

log = logging.getLogger(__name__)

# MultiQC reads at most this many bytes of a log file, and lines of an Iso-Seq cluster CSV
LARGE_FILE_LIMIT = 50_000_000_000

MULTIQC_VERSION = (
    "1.36.dev0/citePhilip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis "
    "results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, "
    "October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354"
)


@dataclass(frozen=True)
class Section:
    """
    One custom content report section. Its custom_data entry, search pattern and report
    order are all derived from this.
    """

    id: str
    title: str
    parent_id: str | None = None
    parent_name: str | None = None
    # File name pattern MultiQC finds the section's data with; heading sections have none
    fn: str | None = None
    file_format: str | None = None
    order: int | None = None

    def custom_data(self) -> dict[str, Any]:
        entry: dict[str, Any] = {"section_name": self.title}
        if self.fn is None:
            # A heading with no data of its own. MultiQC drops an html section with empty data,
            # so it gets an empty paragraph.
            entry.update({"description": "", "plot_type": "html", "data": "<p></p>"})
        if self.file_format:
            entry["file_format"] = self.file_format
        if self.parent_id:
            entry["parent_id"] = self.parent_id
        if self.parent_name:
            entry["parent_name"] = self.parent_name
        return entry


FIXED_SECTIONS = [
    Section("mdg_processing_and_qc", "MDG Processing and QC"),
    Section("secondary_analysis", "Secondary Analysis"),
    Section("sequencing_results", "Sequencing Results"),
    Section(
        "sequencing_info",
        "Sequencing Info",
        "sequencing_results",
        "Sequencing Results",
        fn="sequencing_info.tsv",
        file_format="tsv",
    ),
    Section(
        "flow_cell_metrics", "Flow Cell Metrics", "sequencing_results", "Sequencing Results", fn="summary_clean.csv"
    ),
    Section(
        "sample_info",
        "Sample Information",
        "work_order_info",
        "Work Order Information",
        fn="runinfo_clean.tsv",
        file_format="tsv",
    ),
    Section("read_stats", "Read Stats", "sequencing_results", fn="trimming.logs_clean.csv"),
    Section("megablast_qc", "MegaBlast QC", "megablast_qc_parent", "MegaBlast QC", fn="all_tophits.csv"),
    Section(
        "work_order_submission", "Work Order Submission", "work_order_info", fn="work_order_info.tsv", file_format="tsv"
    ),
    Section(
        "polymerase_rs", "Polymerase Read Stats", "read_stats_parent", "Read Stats", fn="pf_polymerase_read_stats.csv"
    ),
    Section("hifi_rs", "HiFi Read Stats", "read_stats_parent", "Read Stats", fn="hifi_stats.csv"),
]

# Placement of fixed sections relative to each other
FIXED_SECTION_ORDER = {
    "work_order_submission": {"before": "sample_info"},
    "sample_info": {"after": "work_order_submission"},
    "read_stats": {"before": "sequencing_info"},
    "sequencing_info": {"after": "read_stats"},
    "flow_cell_metrics": {"after": "sequencing_info"},
    "polymerase_rs": {"before": "hifi_rs"},
}


def get_config_path(outdir: Path, project: str) -> Path:
    return outdir / f"multiqc_config_{project}.yaml"


def generate_config(run_path: Path, project: str, platform_name: str, folders: list[str], outdir: Path) -> Path:
    """
    Write the MultiQC config for the files MDG_Run_Parser.parse() wrote, and return its path
    """
    platform = MDG_Run_Parser.PLATFORMS[platform_name]
    runinfo = MDG_Run_Parser.load_runinfo(run_path)
    sections = [
        *FIXED_SECTIONS,
        *get_tophits_sections(outdir, folders),
        *get_plot_sections(platform, run_path, runinfo, project, folders),
    ]
    config_data = {
        "title": f"{platform.display_name} Run Report",
        "max_table_rows": 2000,
        "log_filesize_limit": LARGE_FILE_LIMIT,
        "show_analysis_paths": False,
        "show_analysis_time": False,
        "skip_generalstats": True,
        "custom_css_files": ["assets/mdg_theme.css"],
        "template_dark_mode": False,
        "custom_favicon": "assets/md-genomics-2.png",
        "ignore_images": False,
        "thousandsSep_format": ",",
        "table_cond_formatting_rules": {
            "all_columns": {
                "green": [{"s_eq": "green"}],
                "yellow": [{"s_eq": "yellow"}],
                "red": [{"s_eq": "red"}],
            }
        },
        "table_cond_formatting_colours": [{"green": "#4CAF50"}, {"yellow": "#FFC107"}, {"red": "#F44336"}],
        "custom_data": {section.id: section.custom_data() for section in sections},
        "sp": get_search_patterns(sections),
        "report_section_order": get_report_section_order(sections),
        "custom_content": {
            "order": [
                "work_order_info",
                "sequencing_results",
                "read_stats_parent",
                "sequencing_plots",
                "mdg_processing_and_qc",
                "megablast_qc_parent",
                "secondary_analysis",
            ]
        },
        "module_order": ["bcl2fastq", "fastqc", "cutadapt", "trimmomatic", "cellranger", "spaceranger", "isoseq"],
        "remove_sections": ["undetermined_by_lane", "fastqc_per_base_sequence_content"],
        "software_versions": get_software_versions(platform, folders, outdir),
    }
    config_path = get_config_path(outdir, project)
    with open(config_path, "w", encoding="utf-8") as fout:
        yaml.dump(config_data, fout, sort_keys=False, allow_unicode=True)
    return config_path


def get_tophits_sections(outdir: Path, folders: list[str]) -> list[Section]:
    """
    One section per {specimen}_tophits_{n}.csv table. Folders are read in the same order
    MDG_Run_Parser wrote them, so n counts up across folders the same way.
    """
    sections: list[Section] = []
    for folder in folders:
        all_tophits = MDG_Run_Parser.get_folder_data_dir(outdir, folder) / "all_tophits.csv"
        if not all_tophits.is_file():
            continue  # MDG_Run_Parser already warned that this folder has no topHits files
        for specimen_name in pd.read_csv(all_tophits, dtype=str, keep_default_na=False)["Specimen Name"]:
            n = len(sections)
            sections.append(
                Section(
                    f"top_hits_{n}",
                    specimen_name,
                    parent_id="megablast_qc_parent",
                    fn=f"*_tophits_{n}.csv",
                    order=(n + 1) * 10,
                )
            )
    return sections


def get_plot_sections(
    platform: Platform, run_path: Path, runinfo: Runinfo, project: str, folders: list[str]
) -> list[Section]:
    """
    One Sequencing Plots section per plot image
    """
    return [
        Section(
            plot.section_id,
            plot.title,
            "sequencing_plots",
            "Sequencing Plots",
            fn=plot.file_name,
            order=(n + 1) * 10,
        )
        for n, plot in enumerate(platform.plot_images(run_path, runinfo, project, folders))
    ]


def get_search_patterns(sections: list[Section]) -> dict[str, dict[str, Any]]:
    sp: dict[str, dict[str, Any]] = {section.id: {"fn": section.fn} for section in sections if section.fn}
    # Iso-Seq cluster reports can be far larger than MultiQC reads by default
    sp["isoseq/cluster-csv"] = {"num_lines": LARGE_FILE_LIMIT, "max_filesize": LARGE_FILE_LIMIT}
    return sp


def get_report_section_order(sections: list[Section]) -> dict[str, dict[str, Any]]:
    order: dict[str, dict[str, Any]] = dict(FIXED_SECTION_ORDER)
    order.update({section.id: {"order": section.order} for section in sections if section.order is not None})
    return order


def get_software_versions(platform: Platform, folders: list[str], outdir: Path) -> dict[str, str]:
    software_versions = {"MultiQC": MULTIQC_VERSION, **platform.software_versions}
    # MDG_Run_Parser only creates an isoseq folder when a cell has Iso-Seq outputs
    if any((MDG_Run_Parser.get_folder_data_dir(outdir, folder) / "isoseq").is_dir() for folder in folders):
        software_versions["Iso-Seq"] = "1.0"
    return software_versions
