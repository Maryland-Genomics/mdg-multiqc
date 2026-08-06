import yaml
import os
import pandas as pd
import numpy as np

top_hits_fields : list[str] = []
lane_image_fields : list[str] = []
lane_image_fn : list[str] = []
plot_field_and_image_dict : dict = {}

def generate_config(p_path : str, p_project : str, p_runtype : str, p_unaligned_folder : str, p_outdir : str):
    project = p_project
    path  = p_path
    runtype = p_runtype
    unaligned_folder = p_unaligned_folder
    config_data = {
        "title": "Illumina Run Report" if runtype == "Illumina" else "PacBio Run Report",
        # report header info?
        "max_table_rows": 2000,
        "log_filesize_limit": 50000000000,
        "show_analysis_paths": False,
        "show_analysis_time": False,
        "skip_generalstats": True,
        "custom_css_files": ["assets/mdg_theme.css"],
        "template_dark_mode": False,
        "custom_favicon": "assets/md-genomics-2.png",
        "ignore_images": False,
        "thousandsSep_format": ',',
        "table_cond_formatting_rules": {
            "all_columns": {
                "green": [{"s_eq": "green"}],
                "yellow": [{"s_eq": "yellow"}],
                "red": [{"s_eq": "red"}],
            }
        },
        "table_cond_formatting_colours": [{"green": "#4CAF50"}, {"yellow": "#FFC107"}, {"red": "#F44336"}],
        "custom_data": get_custom_data(path, project, runtype, unaligned_folder),
        "sp": get_search_path(runtype),
        "report_section_order": get_report_section_order(runtype),
        "custom_content": {
            "order": ["work_order_info", "sequencing_results", "read_stats_parent", "sequencing_plots",
                      "mdg_processing_and_qc", "megablast_qc_parent", "secondary_analysis"]
        },
        "module_order": ["bcl2fastq", "fastqc", "cutadapt",
                         "trimmomatic","cellranger", "spaceranger", "isoseq"],
        "remove_sections": ["undetermined_by_lane", "fastqc_per_base_sequence_content"],
        "software_versions": get_software_versions(runtype),
    }
    with open(f'{p_outdir}/multiqc_config_{project}.yaml', 'w') as file:
        yaml.dump(config_data, file)

def get_custom_data(path, project, runtype, unaligned_folder) -> dict:
    custom_data = {
        "mdg_processing_and_qc": {
            "section_name": "MDG Processing and QC",
            "description": "",
            "plot_type": "html",
            "data": "<p><p>"
        },
        "secondary_analysis": {
            "section_name": "Secondary Analysis",
            "description": "",
            "plot_type": "html",
            "data": "<p><p>"
        },
        "sequencing_results": {
            "section_name": "Sequencing Results",
            "description": "",
            "plot_type": "html",
            "data": "<p><p>",
        },
        "sequencing_info": {
            "section_name": "Sequencing Info",
            "parent_id": "sequencing_results",
            "parent_name": "Sequencing Results",
        },
        "flow_cell_metrics": {
            "section_name": "Flow Cell Metrics",
            "parent_id": "sequencing_results",
            "parent_name": "Sequencing Results"
        },
        "sample_info": {
            "section_name": "Sample Information",
            "parent_id": "work_order_info",
            "parent_name": "Work Order Information",
        },
        "read_stats": {
            "section_name": "Read Stats",
            "parent_id": "sequencing_results",
        },
        "megablast_qc": {
            "section_name": "MegaBlast QC",
            "parent_id": "megablast_qc_parent",
            "parent_name": "MegaBlast QC"
        },
        "work_order_submission": {
            "section_name": "Work Order Submission",
            "parent_id": "work_order_info",
        },
        "polymerase_rs": {
            "section_name": "Polymerase Read Stats",
            "parent_id": "read_stats_parent",
            "parent_name": "Read Stats",
        },
        "hifi_rs": {
            "section_name": "HiFi Read Stats",
            "parent_id": "read_stats_parent",
            "parent_name": "Read Stats",
        },
    }

    # Add top hits tables (Contingent on Illumina Parser running first)
    if runtype == "Illumina":
        all_tophits = f"{path}/{unaligned_folder}/QA/Project_{project}/MDG_Run_Parsed_Files/all_tophits.csv"
    elif runtype == "PacBio":
        all_tophits = f"{path}/Project_{project}_MDG_Run_Parsed_Files/all_tophits.csv"
    try:
        df = pd.read_csv(all_tophits)
        libraries = df["ID"]

        for (i, topHits) in enumerate(libraries):
            temp_dict = {
                f"top_hits_{i}": {
                    "section_name": topHits,
                    "parent_id": "megablast_qc_parent",
                }
            }
            custom_data.update(temp_dict) 
            top_hits_fields.append(f"top_hits_{i}")
    except Exception as e:
        print(e)

    try:
        # Add lane images
        if runtype == "Illumina":
            for lane in get_lanes(path):
                for image in os.listdir(f"{path}/PA_Logs/images"):
                    index = image.find(f"_{lane}_")
                    if index != -1:
                        field = None
                        section_name = None
                        index += 3
                        image_type = image[index:]
                        match image_type:
                            case "flowcell-Intensity.png":
                                field = f"lane_{lane}_flowcell_intensity"
                                section_name = f"Lane {lane} Flowcell Intensity"
                            case "Intensity-by-cycle_Intensity.png":
                                field = f"lane_{lane}_intensity_by_cycle"
                                section_name = f"Lane {lane} Intensity By Cycle"
                            case "ClusterCount-by-lane.png":
                                field = f"lane_{lane}_cluster_count_by_lane"
                                section_name = f"Lane {lane} Cluster Count By Lane"
                            case "q-heat-map.png":
                                field = f"lane_{lane}_q_heat_map"
                                section_name = f"Lane {lane} Q Heat Map"
                        if field != None and section_name != None:
                            temp_dict = {
                                field : {
                                    "section_name": section_name,
                                    "parent_id": "sequencing_plots",
                                    "parent_name": "Sequencing Plots",
                                }
                            }
                            custom_data.update(temp_dict)
                            lane_image_fields.append(field)
                            lane_image_fn.append(f"_{lane}_{image_type}")
        elif runtype == "PacBio":
            plots = {
                "hexbin_length_plot": "Hexbin Length Plot",
                "base_yield_plot": "Base Yield Plot",
                "ccs_readlength_qv_hist2d": "CCS Read Length QV Histogram",
                "readLenDist0": "Read Length Distribution",
                "ccs_all_readlength_hist_plot": "CCS All Read Length Histogram",
                "ccs_npasses_hist": "CCS Number of Passes Histogram",
                "ccs_accuracy_hist": "CCS Accuracy Histogram",
                "ccs_combined_readlength_hist_plot": "CCS Combined Read Length Histogram",
                "readlength_plot": "Read Length Plot",
                "concordance_plot": "Concordance Plot",
                "5mc_detections_hist": "5mC Detections Histogram",
            }
            for key, value in plots.items():
                temp_dict = {
                    key: {
                        "section_name": value,
                        "parent_id": "sequencing_plots",
                        "parent_name": "Sequencing Plots",
                    }
                }
                custom_data.update(temp_dict)
            plot_fields = list(plots.keys())
            plot_files = list(plots.keys())
            plot_files[plot_files.index("ccs_readlength_qv_hist2d")] += ".hexbin"
            for i in range(len(plot_files)):
                plot_files[i] += ".png"
            plot_dict = dict(zip(plot_fields, plot_files))
            plot_field_and_image_dict.update(plot_dict)
    except Exception as e:
        print(e)
    return custom_data

def get_search_path(runtype):
    sp = {
        "sequencing_info": {
            "fn": "sequencing_info.csv"
        },
        "flow_cell_metrics": {
            "fn": "summary_clean.csv"
        },
        "sample_info": {
            "fn": "runinfo_clean.csv"
        },
        "read_stats": {
            "fn": "trimming.logs_clean.csv"
        },
        "megablast_qc": {
            "fn": "all_tophits.csv"
        },
        "work_order_submission": {
            "fn": "work_order_info.csv"
        },
        "polymerase_rs": {
            "fn": "pf_polymerase_read_stats.csv"
        },
        "hifi_rs": {
            "fn": "hifi_stats.csv"
        },
        "isoseq/cluster-csv": {
            "num_lines": 50000000000,
            "max_filesize": 50000000000,
        }
    }

    for (i, field) in enumerate(top_hits_fields):
        temp_dict = {
            f"top_hits_{i}": {
                "fn": f"*_tophits_{i}.csv"
            }
        }
        sp.update(temp_dict)

    if (runtype == "Illumina"):
        for (i, field) in enumerate(lane_image_fields):
            temp_dict = {
                field: {
                    "fn": f"*{lane_image_fn[i]}"
                }
            }
            sp.update(temp_dict)
    elif (runtype == "PacBio"):
        for key, value in plot_field_and_image_dict.items():
            temp_dict = {
                key: {
                    "fn": value
                }
            }
            sp.update(temp_dict)
    return sp

def get_report_section_order(runtype):
    report_section_order = {
        "work_order_submission": {
            "before": "sample_info"
        },
        "sample_info": {
            "after": "work_order_submission"
        },
        "read_stats": {
            "before": "sequencing_info"
        },
        "sequencing_info": {
            "after": "read_stats",
        },
        "flow_cell_metrics": {
            "after": "sequencing_info",
        },
        "polymerase_rs": {
            "before": "hifi_rs"
        },
    }
    for (i, top_hits) in enumerate(top_hits_fields):
        temp_dict = {
            top_hits: {
                "order": (i + 1) * 10,
                #"before": f"{top_hits_fields[i + 1]}" if i < (len(top_hits_fields) - 1) else ""
            }
        }
        report_section_order.update(temp_dict)

    if runtype == "Illumina":
        for (i, lane) in enumerate(lane_image_fields):
            temp_dict = {
                lane: {
                    "order": (i + 1) * 10
                }
            }
            report_section_order.update(temp_dict)
    elif runtype == "PacBio":
        for i, key in enumerate(plot_field_and_image_dict):
            temp_dict = {
                key: {
                    "order": (i + 1) * 10
                }
            }
    return report_section_order

def get_software_versions(runtype):
    software_versions = {
        "MultiQC": "1.36.dev0/citePhilip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354",
    }
    temp_dict = {}
    if runtype == "Illumina":
        temp_dict = {
            "bcl2fastq": "2.20.0.422",
            "Trimmomatic": "0.33",
        }
    elif runtype == "PacBio":
        temp_dict = {
            "Iso-Seq": "1.0"
        }
    software_versions.update(temp_dict)
    return software_versions

def get_lanes(path):
    # Get every unique lane in run
    runinfo_df = pd.read_csv(f"{path}/PA_Logs/runinfo.csv", header=1)
    lanes = runinfo_df["lane"]
    unique_lanes = np.unique(lanes)
    return unique_lanes