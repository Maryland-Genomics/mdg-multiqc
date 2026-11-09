"""MultiQC module to parse output from Space Ranger count"""

import base64
import json
import logging
import os
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

import numpy as np

from multiqc import BaseMultiqcModule
from multiqc.modules.spaceranger.utils import populate_data_and_headers, set_hidden_cols, transform_data
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)

# Report metric name -> column id for the summary table. Matching is case-insensitive,
# because Space Ranger 3.0+ "key metrics" use sentence case ("Mean reads per cell").
SUMMARY_COLUMNS = {
    "Number of Reads": "reads",
    "Number of Spots Under Tissue": "spots under tissue",
    "Mean Reads per Spot": "avg reads/spot",
    "Fraction Reads in Spots Under Tissue": "reads in spots",
    "Fraction Reads in Squares Under Tissue": "reads in spots",
    "Genes Detected": "genes detected",
    "Total Genes Detected": "genes detected",
    "Total Genes Detected in Cells": "genes detected",
    "Median Genes per Spot": "median genes/spot",
    "Median UMI Counts per Spot": "median umi/spot",
    "Median UMIs per Spot": "median umi/spot",
    "Number of Cells": "cells",
    "Number of Cells Under Tissue": "cells",
    "Reads in Cells": "reads in cells",
    "UMIs in Cells": "umis in cells",
    "Mean Reads per Cell": "avg reads/cell",
    "Median Genes per Cell": "median genes/cell",
    "Median UMIs per Cell": "median umi/cell",
    "Valid Barcodes": "valid bc",
    "Valid UMIs": "valid umi",
    "Valid UMI Sequences": "valid umi",
    "Sequencing Saturation": "saturation",
    "Q30 Bases in Barcode": "Q30 bc",
    "Q30 Bases in UMI": "Q30 UMI",
    "Q30 Bases in RNA Read": "Q30 read",
    "Q30 Bases in Probe Read": "Q30 read",
    "Reads Mapped to Probe Set": "reads mapped",
    "Reads Mapped Confidently to Probe Set": "confident reads",
    "Reads Mapped Confidently to Filtered Probe Set": "confident filtered reads",
    "Reads Mapped Confidently to the Filtered Probe Set": "confident filtered reads",
    "Estimated UMIs from Genomic DNA": "genomic umis",
    "Estimated UMIs from Genomic DNA per Unspliced Probe": "genomic umis/unspliced probe",
}

GENERAL_STATS_COLUMNS = {
    "Number of Spots Under Tissue": "spots under tissue",
    "Number of Cells": "cells",
    "Number of Cells Under Tissue": "cells",
    "Mean Reads per Spot": "avg reads/spot",
    "Fraction Reads in Spots Under Tissue": "reads in spots",
    "Fraction Reads in Squares Under Tissue": "reads in spots",
    "Number of Reads": "reads",
    "Valid Barcodes": "valid bc",
}

# Visium HD "Bin Metrics Overview" rows: report metric name -> column id prefix.
# One column is made per bin size, e.g. "bins under tissue 8um".
BIN_COLUMNS = {
    "Number of Bins Under Tissue": "bins under tissue",
    "Mean UMI Counts per Bin": "mean umi/bin",
    "Mean Genes per Bin": "mean genes/bin",
}

_CANONICAL_NAMES = {name.lower(): name for name in SUMMARY_COLUMNS}


def parse_count_html(module: BaseMultiqcModule):
    """
    Space Ranger count report parser. Handles both the pre-3.0 report layout
    (`summary_tab` / `analysis_tab`) and the 3.0+ layout (`tabs.tab_data` cards).
    """

    general_stats_data: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict()
    summary_data_by_sample: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict()
    warnings_data_by_sample: Dict[str, Dict[str, Union[str, float, int, None]]] = defaultdict(lambda: defaultdict())

    plots_data: Dict[str, Dict] = {"saturation": defaultdict(), "genes": defaultdict(), "genomic_dna": defaultdict()}
    plots_conf: Dict[str, Dict] = {"saturation": defaultdict(), "genes": defaultdict(), "genomic_dna": defaultdict()}

    warnings_headers: Dict = dict()
    summary_headers: Dict[str, Dict[str, Union[str, float, int, None]]] = {
        "reads": {
            "rid": "count_data_reads",
            "title": "Reads",
            "description": "Number of reads",
            "shared_key": "read_count",
        }
    }
    general_stats_headers: Dict[str, Dict[str, Union[str, float, int, None]]] = {
        "reads": {
            "rid": "count_genstats_reads",
            "title": "Reads",
            "description": "Number of reads",
            "shared_key": "read_count",
            "namespace": "Space Ranger Count",
        }
    }

    for f in module.find_log_files("spaceranger/count_html", filehandles=True):
        fname = os.path.join(f["root"], f["fn"])
        data = _read_report_json(f["f"])
        if data is None:
            log.error(f"Couldn't find JSON summary data in HTML report, skipping: {fname}")
            continue

        if "summary" in data:
            report = _parse_legacy_report(data["summary"], fname)
        elif "tabs" in data:
            report = _parse_v3_report(data, fname)
        else:
            raise ValueError(
                f"Unrecognised Space Ranger report layout in {fname}: "
                f"expected a 'summary' or 'tabs' key, found {sorted(data)}"
            )

        sample_name = module.clean_s_name(report["sample_id"], f)
        if sample_name in general_stats_data:
            log.debug(f"Duplicate sample name found in {fname}! Overwriting: {sample_name}")

        software_name, software_version = _split_pipeline_version(report["pipeline_version"])
        module.add_software_version(version=software_version, sample=sample_name, software_name=software_name)

        data_rows = [[_CANONICAL_NAMES.get(name.lower(), name), value] for name, value in report["data_rows"]]

        general_stats_data[sample_name] = populate_data_and_headers(
            headers_to_update=general_stats_headers,
            new_data=data_rows,
            new_headers=GENERAL_STATS_COLUMNS,
            colors={
                "spots under tissue": "RdPu",
                "cells": "RdPu",
                "avg reads/spot": "Blues",
                "reads in spots": "PiYG",
                "reads": "YlGn",
                "valid bc": "RdYlGn",
            },
            prefix="Count",
            int_cols=[
                "reads",
                "spots under tissue",
                "cells",
                "avg reads/spot",
            ],
        )
        if len(general_stats_data[sample_name]) == 0:
            continue

        summary_data_by_sample[sample_name] = populate_data_and_headers(
            headers_to_update=summary_headers,
            new_data=data_rows,
            new_headers=SUMMARY_COLUMNS,
            colors={
                "reads": "YlGn",
                "spots under tissue": "RdPu",
                "cells": "RdPu",
                "avg reads/spot": "Blues",
                "avg reads/cell": "Blues",
                "genes detected": "Greens",
                "median genes/spot": "Purples",
                "median genes/cell": "Purples",
                "reads in spots": "PuBuGn",
                "valid bc": "Spectral",
                "valid umi": "RdYlGn",
                "median umi/spot": "YlGn",
                "median umi/cell": "YlGn",
                "saturation": "PRGn",
                "genomic umis": "PuRd",
                "genomic umis/unspliced probe": "YlOrRd",
            },
            prefix="Count",
            int_cols=[
                "reads",
                "spots under tissue",
                "cells",
                "avg reads/spot",
                "median genes/spot",
                "median genes/cell",
                "median umi/cell",
                "genes detected",
                "genomic umis/unspliced probe",
            ],
        )
        if len(summary_data_by_sample[sample_name]) == 0:
            continue

        summary_data_by_sample[sample_name].update(_add_bin_columns(summary_headers, report["bin_rows"]))

        for alarm in report["alarms"]:
            # "Intron mode used" alarm added in Space Ranger 7.0 lacks id
            if "id" not in alarm:
                continue
            alarm_id = alarm["id"]
            warnings_data_by_sample[sample_name][alarm_id] = "FAIL"
            warnings_headers[alarm_id] = {
                "title": alarm_id.replace("_", " ").title(),
                "description": alarm["title"],
                "bgcols": {"FAIL": "#e5001336"},
            }

        for plot_key, plot in report["plots"].items():
            if not plots_conf.get(plot_key):
                plots_conf[plot_key] = {
                    "config": {
                        "id": plot["id"],
                        "title": f"Space Ranger: count: {plot['title']}",
                        "xlab": plot["xlab"],
                        "ylab": plot["ylab"],
                        "ylog": False,
                        "xlog": False,
                        **plot.get("extra_config", {}),
                    },
                    "description": plot["description"],
                    "helptext": plot["helptext"],
                }
            plots_data[plot_key][sample_name] = transform_data(plot["data"])

        module.add_data_source(f, sample_name, module="spaceranger", section="count")

    summary_data_by_sample = module.ignore_samples(summary_data_by_sample)
    warnings_data_by_sample = module.ignore_samples(warnings_data_by_sample)
    general_stats_data = module.ignore_samples(general_stats_data)
    for k in plots_data.keys():
        plots_data[k] = module.ignore_samples(plots_data[k])

    summary_headers = set_hidden_cols(
        summary_headers,
        [
            "Q30 bc",
            "Q30 UMI",
            "Q30 read",
            "reads in spots",
            "avg reads/spot",
            "confident reads",
            "confident filtered reads",
            "genomic unis/unspliced probe",
            "valid umi",
            "saturation",
            "reads in cells",
            "umis in cells",
        ]
        + [col_id for col_id in summary_headers if col_id.startswith(("mean umi/bin", "mean genes/bin"))],
    )

    module.general_stats_addcols(general_stats_data, general_stats_headers)

    # Add sections to the report
    if len(warnings_data_by_sample) > 0:
        module.add_section(
            name="Count - Warnings",
            anchor="spaceranger-count-warnings-section",
            description="Warnings encountered during the analysis",
            plot=table.plot(
                warnings_data_by_sample,
                warnings_headers,
                {
                    "namespace": "Space Ranger Count",
                    "id": "spaceranger-count-warnings",
                    "title": "Space Ranger: Count: Warnings",
                },
            ),
        )

    module.add_section(
        name="Count - Summary stats",
        anchor="spaceranger-count-stats-section",
        description="Summary QC metrics from Space Ranger count",
        plot=table.plot(
            summary_data_by_sample,
            summary_headers,
            {
                "namespace": "Space Ranger Count",
                "id": "spaceranger-count-stats",
                "title": "Space Ranger: Count: Summary stats",
            },
        ),
    )

    # The gDNA plots are only contained for spaceranger workflows with probesets
    if plots_conf.get("genomic_dna"):
        module.add_section(
            name="Count - UMIs from Genomic DNA",
            anchor="spaceranger-count-bcrank-plot-section",
            description=plots_conf["genomic_dna"]["description"],
            helptext=plots_conf["genomic_dna"]["helptext"],
            plot=linegraph.plot(
                plots_data["genomic_dna"],
                plots_conf["genomic_dna"]["config"],
            ),
        )

    if plots_conf.get("genes"):
        module.add_section(
            name="Count - Median genes",
            anchor="spaceranger-count-genes-plot-section",
            description=plots_conf["genes"]["description"],
            helptext=plots_conf["genes"]["helptext"],
            plot=linegraph.plot(plots_data["genes"], plots_conf["genes"]["config"]),
        )

    if plots_conf.get("saturation"):
        module.add_section(
            name="Count - Saturation plot",
            anchor="spaceranger-count-saturation-plot-section",
            description=plots_conf["saturation"]["description"],
            helptext=plots_conf["saturation"]["helptext"],
            plot=linegraph.plot(
                plots_data["saturation"],
                plots_conf["saturation"]["config"],
            ),
        )

    # Write parsed report data to a file
    module.write_data_file(summary_data_by_sample, "multiqc_spaceranger_count")

    return len(general_stats_data)


def _read_report_json(fh) -> Optional[Dict[str, Any]]:
    """Return the `const data = {...}` JSON object embedded in a web_summary.html, or None."""
    for line in fh:
        line = line.strip()
        if line.startswith("const data"):
            return json.loads(line.replace("const data = ", "", 1).rstrip(";"))
    return None


def _split_pipeline_version(pipeline_version: str):
    """
    Split a "Pipeline Version" value into (software name, version).
    Older reports use "spaceranger-2.1.1"; newer ones use a bare "4.1.0".
    """
    pipeline_version = str(pipeline_version).strip()
    if "-" in pipeline_version:
        software_name, software_version = pipeline_version.split("-", 1)
        return software_name, software_version
    return "spaceranger", pipeline_version.rsplit(maxsplit=1)[-1].removeprefix("v")


def _parse_legacy_report(summary: Dict[str, Any], fname: str) -> Dict[str, Any]:
    """Extract metrics and plots from a pre-3.0 Space Ranger report (`summary_tab` / `analysis_tab`)."""
    pipeline_version = next(
        iter(x[1] for x in summary["summary_tab"]["pipeline_info_table"]["rows"] if x[0] == "Pipeline Version")
    )

    # List of data collated from different tables in the report: [metric name, value]
    data_rows = (
        [["Number of Spots Under Tissue", summary["summary_tab"]["filtered_bcs_transcriptome_union"]["metric"]]]
        + summary["summary_tab"]["cells"]["table"]["rows"]
        + summary["summary_tab"]["sequencing"]["table"]["rows"]
        + summary["summary_tab"]["mapping"]["table"]["rows"]
    )
    # This is only contained in spaceranger reports for analyses with probe sets
    try:
        data_rows.extend(summary["analysis_tab"]["gdna"]["gems"]["table"]["rows"])
    except KeyError as e:
        log.debug(f"Could not parse the expected DNA table in the spaceranger report: {fname}: '{e}' field is missing")

    plots: Dict[str, Dict[str, Any]] = {}
    # `analysis_tab` may not be present in the report if there are few reads
    try:
        sat = summary["analysis_tab"]["seq_saturation_plot"]
        plots["saturation"] = {
            "id": "mqc_spaceranger_count_saturation",
            "title": sat["help"]["title"],
            "xlab": sat["plot"]["layout"]["xaxis"]["title"],
            "ylab": sat["plot"]["layout"]["yaxis"]["title"],
            "extra_config": {"ymin": 0, "ymax": 1},
            "description": "Sequencing saturation",
            "helptext": sat["help"]["helpText"],
            "data": sat["plot"]["data"][0],
        }
    except KeyError as e:
        log.debug(f"No saturation plot found in the spaceranger report {fname}: {e}")

    try:
        genes = summary["analysis_tab"]["median_gene_plot"]
        plots["genes"] = {
            "id": "mqc_spaceranger_count_genesXspot",
            "title": genes["help"]["title"],
            "xlab": genes["plot"]["layout"]["xaxis"]["title"],
            "ylab": genes["plot"]["layout"]["yaxis"]["title"],
            "description": "Median gene counts per spot",
            "helptext": genes["help"]["helpText"],
            "data": genes["plot"]["data"][0],
        }
    except KeyError:
        log.debug(f"No median gene plot found in the spaceranger report {fname}")

    # The gDNA plots are only contained for spaceranger workflows with probesets
    try:
        gdna = summary["analysis_tab"]["gdna"]
        plots["genomic_dna"] = {
            "id": "mqc_spaceranger_count_genomic_dna",
            "title": gdna["gems"]["help"]["title"],
            "xlab": gdna["plot"]["layout"]["xaxis"]["title"],
            "ylab": gdna["plot"]["layout"]["yaxis"]["title"],
            "description": "Estimated UMIs from Genomic DNA per Unspliced Probe",
            "helptext": gdna["gems"]["help"]["data"][2][1][0]
            + "\n\nThis summary graphic in the MultiQC report only shows the estimated mean baseline level of unspliced probe counts.",
            "data": gdna["plot"]["data"][2],
        }
    except KeyError:
        log.debug(f"No genomic DNA plot found in the spaceranger report {fname}")

    return {
        "sample_id": summary["sample"]["id"],
        "pipeline_version": pipeline_version,
        "data_rows": data_rows,
        "bin_rows": [],
        "alarms": summary["alarms"].get("alarms", []),
        "plots": plots,
    }


def _parse_v3_report(data: Dict[str, Any], fname: str) -> Dict[str, Any]:
    """
    Extract metrics and plots from a Space Ranger 3.0+ report, where everything lives in
    cards under `tabs.tab_data`. Cards are looked up by key rather than tab position,
    because which tabs exist depends on the assay (Visium vs Visium HD, segmentation, etc).
    """
    tabs: List[Dict[str, Any]] = data["tabs"]["tab_data"]

    def find_card(key: str) -> Optional[Dict[str, Any]]:
        return next((tab[key] for tab in tabs if tab.get(key) is not None), None)

    summary_tab = find_card("top")
    run_summary = find_card("run_summary")
    if summary_tab is None or run_summary is None:
        raise ValueError(f"Space Ranger 3.0+ report {fname} is missing its 'top' or 'run_summary' card")
    left = summary_tab["left"]

    pipeline_version = next(x[1] for x in run_summary["card"]["inner"]["rows"] if x[0] == "Pipeline Version")

    # Key metrics are a dict of {"name": ..., "metric": ...}; the other cards are [name, value] rows
    data_rows: List[List[str]] = [[m["name"], m["metric"]] for m in left["key_metrics"]["card"]["inner"].values()]
    data_rows += left["sequencing_metrics"]["inner"]["rows"]
    data_rows += left["mapping_metrics"]["inner"]["rows"]
    # Cell segmentation is only present for Visium HD runs with segmentation
    if left.get("segmentation_summary") is not None:
        data_rows += left["segmentation_summary"]["card"]["inner"]["rows"]

    # The gDNA card is only present for analyses with probe sets
    gdna_card = find_card("genomic_dna")
    if gdna_card is not None:
        data_rows += gdna_card["card"]["inner"]["table"]["rows"]

    # Visium HD only: a table with one column per bin size
    bin_card = find_card("bin_metrics")
    bin_rows = bin_card["card"]["inner"]["rows"] if bin_card is not None else []

    plots: Dict[str, Dict[str, Any]] = {}
    saturation_card = find_card("sequencing_saturation")
    if saturation_card is not None:
        depth = saturation_card["sequencing_depth"]["inner"]
        if depth.get("seq_saturation_plot") is not None:
            sat = depth["seq_saturation_plot"]
            plots["saturation"] = {
                "id": "mqc_spaceranger_count_saturation",
                "title": sat["title"]["title"],
                "xlab": _axis_title(sat["inner"]["layout"]["xaxis"]),
                "ylab": _axis_title(sat["inner"]["layout"]["yaxis"]),
                "extra_config": {"ymin": 0, "ymax": 1},
                "description": "Sequencing saturation",
                "helptext": sat["title"]["helpText"],
                "data": _plotly_trace_xy(sat["inner"]["data"][0]),
            }
        if depth.get("genes_plot") is not None:
            genes = depth["genes_plot"]
            plots["genes"] = {
                "id": "mqc_spaceranger_count_genesXspot",
                "title": genes["title"]["title"],
                "xlab": _axis_title(genes["inner"]["layout"]["xaxis"]),
                "ylab": _axis_title(genes["inner"]["layout"]["yaxis"]),
                "description": genes["title"]["title"],
                "helptext": genes["title"]["helpText"],
                "data": _plotly_trace_xy(genes["inner"]["data"][0]),
            }

    if gdna_card is not None:
        gdna_plot = gdna_card["card"]["inner"]["plot"]
        gdna_help = {row[0]: row[1][0] for row in gdna_card["card"]["title"]["data"]}
        plots["genomic_dna"] = {
            "id": "mqc_spaceranger_count_genomic_dna",
            "title": gdna_card["card"]["title"]["title"],
            "xlab": _axis_title(gdna_plot["layout"]["xaxis"]),
            "ylab": _axis_title(gdna_plot["layout"]["yaxis"]),
            "description": "Estimated UMIs from Genomic DNA per Unspliced Probe",
            "helptext": gdna_help["Estimated UMIs from Genomic DNA per Unspliced Probe"]
            + "\n\nThis summary graphic in the MultiQC report only shows the estimated mean baseline level of unspliced probe counts.",
            "data": _plotly_trace_xy(next(trace for trace in gdna_plot["data"] if trace.get("name") == "Mean")),
        }

    return {
        "sample_id": data["sample"]["id"],
        "pipeline_version": pipeline_version,
        "data_rows": data_rows,
        "bin_rows": bin_rows,
        # 3.0+ alarms have no id, so derive a stable one from the title
        "alarms": [
            {"id": re.sub(r"[^a-z0-9]+", "_", alarm["title"].lower()).strip("_"), "title": alarm["title"]}
            for alarm in data["alarms"]["alarms"]
        ],
        "plots": plots,
    }


def _add_bin_columns(headers: Dict[str, Dict[str, Union[str, float, int, None]]], bin_rows: List[List[str]]):
    """
    Turn the Visium HD "Bin Metrics Overview" table (first row is the bin sizes, then one
    row per metric with one value per bin size) into one column per metric and bin size.
    """
    if not bin_rows:
        return {}
    bin_sizes = bin_rows[0][1:]
    values: Dict[str, Union[str, float, int, None]] = {}
    for metric_name, *metric_values in bin_rows[1:]:
        if metric_name not in BIN_COLUMNS:
            continue
        for bin_size, value in zip(bin_sizes, metric_values):
            size_id = bin_size.replace(" ", "").replace("µ", "u")
            col_id = f"{BIN_COLUMNS[metric_name]} {size_id}"
            values[col_id] = float(value.replace(",", ""))
            if col_id not in headers:
                headers[col_id] = {
                    "rid": f"Count_{col_id.replace(' ', '_').replace('/', '_')}",
                    # Built by hand: str.title() would turn "µm" into "Μm"
                    "title": f"{BIN_COLUMNS[metric_name].title().replace('Umi', 'UMI')} ({bin_size})",
                    "description": f"{metric_name} ({bin_size} bins)",
                    "namespace": "Space Ranger Count",
                    "scale": "GnBu",
                    "format": "{:,.0f}" if metric_name == "Number of Bins Under Tissue" else "{:,.1f}",
                }
    return values


def _axis_title(axis: Dict[str, Any]) -> str:
    """Plotly axis titles are either a plain string or {"text": ...}."""
    title = axis["title"]
    return title["text"] if isinstance(title, dict) else title


def _plotly_trace_xy(trace: Dict[str, Any]) -> Dict[str, list]:
    """Return {"x": [...], "y": [...]} for a Plotly trace, decoding binary-packed arrays."""
    return {"x": _plotly_array(trace["x"]), "y": _plotly_array(trace["y"])}


def _plotly_array(values: Union[list, Dict[str, str]]) -> list:
    """
    Plotly 6 serialises numeric arrays as {"bdata": <base64>, "dtype": "f8"} (little-endian
    typed array) instead of a JSON list. Decode those; pass plain lists through.
    """
    if isinstance(values, list):
        return values
    if "shape" in values and "," in str(values["shape"]):
        raise ValueError(f"Expected a 1-D Plotly array, got shape {values['shape']}")
    return np.frombuffer(base64.b64decode(values["bdata"]), dtype=np.dtype("<" + values["dtype"])).tolist()
