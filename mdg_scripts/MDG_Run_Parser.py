import csv
import json
import pandas as pd
import os
import shutil
import xml.etree.ElementTree as ET
import re
import yaml

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
    generate_read_stats(path, run_data_dir, project, sequencing_platform)
    if sequencing_platform == "illumina":
        copy_illumina_plots(path, runinfo, run_data_dir, project)
        parse_sequencing_info(path, run_data_dir)
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
            path, folder_data_dir, project, sequencing_platform, folder, tophits_count
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


def parse_sequencing_info(path, mdg_parsed_dir):
    """
    Creates a file to store parsed sequencing info in the illumina folder
    """
    # Copy files
    run_param = shutil.copy(f"{path}/RunParameters.xml", mdg_parsed_dir)
    run_info = shutil.copy(f"{path}/RunInfo.xml", mdg_parsed_dir)

    # Get roots
    param_tree = ET.parse(run_param)
    info_tree = ET.parse(run_info)
    param_root = param_tree.getroot()
    info_root = info_tree.getroot()

    # Set up dataFrame items
    data = {
        "Run ID": [],
        "Platform": [],
        "Instrument": [],
        "Flow Cell": [],
        "Run Recipe": [],
    }

    # Assign dataFrame items
    data["Run ID"].append(info_root.find("Run").attrib.get("Id"))
    data["Instrument"].append(param_root.find("InstrumentType").text)
    if param_root.find("InstrumentSerialNumber").text == "VL00448":
        data["Flow Cell"].append(param_root.find("FlowCellMode").text)
    elif param_root.find("InstrumentSerialNumber").text == "LH01184":
        data["Flow Cell"].append(param_root.find("RecipeVersion").text)
    data["Platform"].append("Illumina")

    # Parse run recipe to string
    run_recipe = ""
    planned_cycles = param_root.find("PlannedCycles")
    if planned_cycles == None:
        planned_cycles = param_root.find("PlannedReads")
        for i, cycle in enumerate(planned_cycles):
            if not i == 0:
                run_recipe += ":"
            run_recipe += cycle.get("Cycles")
    else:
        for i, cycle in enumerate(planned_cycles):
            if not i == 0:
                run_recipe += ":"
            run_recipe += cycle.text
    data["Run Recipe"].append(run_recipe)

    # Remove copies
    os.remove(run_param)
    os.remove(run_info)

    # Convert to csv
    df = pd.DataFrame(data)
    df.to_csv(f"{mdg_parsed_dir}/sequencing_info.csv", index=False)


def parse_flowcell_info(path, mdg_parsed_dir):
    """
    Creates a file to store parsed flowcell info in the illumina folder
    """
    # Copy file (temporarily, since MultiQC recognizes summary.csv and replaces custom content module with Illumina Module)
    flow_cell_info = shutil.copy(f"{path}/PA_Logs/summary.csv", mdg_parsed_dir)

    with open(flow_cell_info, "r") as input:
        with open(f"{mdg_parsed_dir}/temp.csv", "w") as output:
            data = input.read().splitlines(True)
            blank_index = data.index("\n")
            # iterate all lines from file
            for line in range(2, blank_index):
                data[line] = data[line].replace("nan", "")
                output.write(data[line])
            
    # Generate cleaned file, remove copy
    clean_data = f"{mdg_parsed_dir}/summary_clean.csv"
    os.replace(f"{mdg_parsed_dir}/temp.csv", clean_data)
    os.remove(flow_cell_info)


def parse_megablastqc_info(path, mdg_parsed_dir, project, sequencing_platform, folder, tophits_count=0):
    """
    Creates a file to store parsed MegaBlast QC info for one Unaligned folder (Illumina) or cell (PacBio).
    tophits_count is the number of tophits files already written by earlier calls; the
    updated count is returned so numbering stays unique across folders.
    """
    # Set up dataframe
    data = {
        "ID": [],
        "Specimen Name": [],
        "Expected Species": [],
        "Contamination Status": [],
        "Total Reads": [],
        "Total Reads With Hits": [],
        "Total Reads With Hits Percent": [],
        "Genus Match Hits": [],
        "Genus Match Hits Percent": [],
        "Other Genus Hits": [],
        "Other Genus Hits Percent": [],
        "Species Match Hits": [],
        "Species Match Hits Percent": [],
        "Other Species Hits": [],
        "Other Species Hits Percent": [],
        }

    tophits_directory = None
    if (sequencing_platform == "illumina"):
        tophits_directory = f"{path}/{folder}/QC/Project_{project}"
    elif (sequencing_platform == "pacbio"):
        tophits_directory = f"{path}/{folder}/hifi_reads/QC"

    start_count = tophits_count

    try:
        if sequencing_platform == "illumina":
            for sample in os.listdir(tophits_directory):
                if sample.startswith("Sample"):
                    # Assume file called *.topHits.txt, get all files
                    for tophits_file in os.listdir(f"{tophits_directory}/{sample}/QC"):
                        # Get all files in directory that end with topHits.txt
                        if tophits_file.endswith("topHits.txt"):
                            data, tophits_count = get_tophits(data, f"{tophits_directory}/{sample}/QC", tophits_file, mdg_parsed_dir, project, tophits_count, sequencing_platform)
        elif sequencing_platform == "pacbio":
            for tophits_file in os.listdir(tophits_directory):
                if tophits_file.endswith("topHits.txt"):
                    data, tophits_count = get_tophits(data, tophits_directory, tophits_file, mdg_parsed_dir, project, tophits_count, sequencing_platform)

        # Write data to new csv file. Row numbers continue from earlier Unaligned folders: MultiQC
        # merges every all_tophits.csv into one table keyed on this column, so repeats would overwrite rows.
        df = pd.DataFrame(data)
        df.index = range(start_count + 1, start_count + 1 + len(df))
        df.to_csv(f"{mdg_parsed_dir}/all_tophits.csv", index=True, index_label="No.")
    except Exception as e:
        print(e)
    return tophits_count


def get_tophits(data, tophits_directory, tophits_file, mdg_parsed_dir, project, tophits_count, sequencing_platform):
    try:
        with open(f"{tophits_directory}/{tophits_file}") as fin:
            tophits = fin.read().splitlines(True)
    
        # Parse txt contents
        if sequencing_platform == "illumina":
            id_index_start = tophits[0].find(":", tophits[0].find(":") + 1) + 1
            id_index_end = tophits[0].find("lane:") + 6
            id = tophits[0][id_index_start:id_index_end]
            id = id.replace(":lane:", "-")
            data["ID"].append(id)
        elif sequencing_platform == "pacbio":
            id_index_start = tophits_file.find("reads") + 6
            id_index_end = tophits_file.find(".subset")
            id = tophits_file[id_index_start:id_index_end]
            data["ID"].append(id)

        spec_name_index_start = tophits[0].find(project) + len(project) + 1
        spec_name_index_end = tophits[0].find(":", spec_name_index_start)
        spec_name = tophits[0][spec_name_index_start:spec_name_index_end]
        data["Specimen Name"].append(spec_name)

        th_expected_spec_index = index_containing_substring(tophits, "Sp.Description:") # th_ indexes refer to those of tophits
        expected_spec_index = tophits[th_expected_spec_index].find(":") + 2
        expected_spec = tophits[th_expected_spec_index][expected_spec_index:]
        expected_spec = expected_spec.replace("\n", "")
        data["Expected Species"].append(expected_spec)

        th_contamination_status_index = index_containing_substring(tophits, "Contamination status:")
        contaminated_status_index = tophits[th_contamination_status_index].find(":") + 2
        contaminated_status = tophits[th_contamination_status_index][contaminated_status_index:]
        contaminated_status = contaminated_status.replace("\n", "")
        data["Contamination Status"].append(contaminated_status)

        th_total_reads_index = index_containing_substring(tophits, "Total_reads:")
        total_reads_index = tophits[th_total_reads_index].find(":") + 2
        total_reads = tophits[th_total_reads_index][total_reads_index:]
        total_reads = total_reads.replace("\n", "")
        data["Total Reads"].append(total_reads)

        th_total_reads_wh_index = index_containing_substring(tophits, "Total reads with hits:")
        total_reads_wh_index_start = tophits[th_total_reads_wh_index].find(":") + 2
        total_reads_wh_index_end = tophits[th_total_reads_wh_index].find("(")
        total_reads_wh = tophits[th_total_reads_wh_index][total_reads_wh_index_start:total_reads_wh_index_end]
        data["Total Reads With Hits"].append(total_reads_wh)

        total_reads_wh_p_index_start = total_reads_wh_index_end + 1
        total_reads_wh_p_index_end = tophits[th_total_reads_wh_index].find(")")
        total_reads_wh_p = tophits[th_total_reads_wh_index][total_reads_wh_p_index_start:total_reads_wh_p_index_end]
        data["Total Reads With Hits Percent"].append(total_reads_wh_p)

        th_genus_match_hits_index = index_containing_substring(tophits, "Genus match hits:")
        genus_match_hits_index_start = tophits[th_genus_match_hits_index].find(":") + 2
        genus_match_hits_index_end = tophits[th_genus_match_hits_index].find("(")
        genus_match_hits = tophits[th_genus_match_hits_index][genus_match_hits_index_start:genus_match_hits_index_end]
        data["Genus Match Hits"].append(genus_match_hits)

        genus_match_hits_p_index_start = genus_match_hits_index_end + 1
        genus_match_hits_p_index_end = tophits[th_genus_match_hits_index].find(")")
        genus_match_hits_p = tophits[th_genus_match_hits_index][genus_match_hits_p_index_start:genus_match_hits_p_index_end]
        data["Genus Match Hits Percent"].append(genus_match_hits_p)

        th_other_genus_hits_index = index_containing_substring(tophits, "Other genus hits:")
        other_genus_hits_index_start = tophits[th_other_genus_hits_index].find(":") + 2
        other_genus_hits_index_end = tophits[th_other_genus_hits_index].find("(")
        other_genus_hits = tophits[th_other_genus_hits_index][other_genus_hits_index_start:other_genus_hits_index_end]
        data["Other Genus Hits"].append(other_genus_hits)

        other_genus_hits_p_index_start = other_genus_hits_index_end + 1
        other_genus_hits_p_index_end = tophits[th_other_genus_hits_index].find(")")
        other_genus_hits_p = tophits[th_other_genus_hits_index][other_genus_hits_p_index_start:other_genus_hits_p_index_end]
        data["Other Genus Hits Percent"].append(other_genus_hits_p)

        th_species_match_hits_index = index_containing_substring(tophits, "Species match hits:")
        species_match_hits_index_start = tophits[th_species_match_hits_index].find(":") + 2
        species_match_hits_index_end = tophits[th_species_match_hits_index].find("(")
        species_match_hits = tophits[th_species_match_hits_index][species_match_hits_index_start:species_match_hits_index_end]
        data["Species Match Hits"].append(species_match_hits)

        species_match_hits_p_index_start = species_match_hits_index_end + 1
        species_match_hits_p_index_end = tophits[th_species_match_hits_index].find(")")
        species_match_hits_p = tophits[th_species_match_hits_index][species_match_hits_p_index_start:species_match_hits_p_index_end]
        data["Species Match Hits Percent"].append(species_match_hits_p)

        th_other_species_hits_index = index_containing_substring(tophits, "Other species hits:")
        other_species_hits_index_start = tophits[th_other_species_hits_index].find(":") + 2
        other_species_hits_index_end = tophits[th_other_species_hits_index].find("(")
        other_species_hits = tophits[th_other_species_hits_index][other_species_hits_index_start:other_species_hits_index_end]
        data["Other Species Hits"].append(other_species_hits)

        other_species_hits_p_index_start = other_species_hits_index_end + 1
        other_species_hits_p_index_end = tophits[th_other_species_hits_index].find(")")
        other_species_hits_p = tophits[th_other_species_hits_index][other_species_hits_p_index_start:other_species_hits_p_index_end]
        data["Other Species Hits Percent"].append(other_species_hits_p)

        # Get top 50 results
        top_50_data = {
            "Read Count": [],
            "Percent Reads": [],
            "Percent Hit Reads": [],
            "Species": [],
        }

        th_top_50_index = index_containing_substring(tophits, "Top 50 species hits:") + 3
        cur_top_50_index = th_top_50_index
        for i in range(50):
            try:
                top_50_reads_count_index_start = 0
                top_50_reads_count_index_end =  tophits[cur_top_50_index + i].find("	")
                top_50_reads_count = tophits[cur_top_50_index + i][top_50_reads_count_index_start:top_50_reads_count_index_end]
                top_50_data["Read Count"].append(top_50_reads_count.strip())

                top_50_percent_reads_index_start = top_50_reads_count_index_end + 1
                top_50_percent_reads_index_end = tophits[cur_top_50_index + i].find("	", top_50_percent_reads_index_start)
                top_50_percent_reads = tophits[cur_top_50_index + i][top_50_percent_reads_index_start:top_50_percent_reads_index_end]
                top_50_data["Percent Reads"].append(top_50_percent_reads.strip())

                top_50_percent_hit_reads_index_start = top_50_percent_reads_index_end + 1
                top_50_percent_hit_reads_index_end = tophits[cur_top_50_index + i].find("	", top_50_percent_hit_reads_index_start)
                top_50_percent_hit_reads = tophits[cur_top_50_index + i][top_50_percent_hit_reads_index_start:top_50_percent_hit_reads_index_end]
                top_50_data["Percent Hit Reads"].append(top_50_percent_hit_reads.strip())

                top_50_species_index_start = top_50_percent_hit_reads_index_end + 1
                top_50_species_index_end = tophits[cur_top_50_index + i].find("	", top_50_species_index_start)
                top_50_species = tophits[cur_top_50_index + i][top_50_species_index_start:top_50_species_index_end]
                top_50_species = top_50_species.replace(",", " ")
                top_50_data["Species"].append(top_50_species.strip())
            except IndexError:
                    break
        df = pd.DataFrame(top_50_data)
        df.to_csv(f"{mdg_parsed_dir}/{spec_name}_tophits_{tophits_count}.csv", index=True)
        tophits_count += 1
    except Exception as e:
        print(e)
    return data, tophits_count


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
    Copies trimming.logs.csv into the illumina folder and trims file to be read
    """
    if sequencing_platform == "illumina":
        generate_illumina_read_stats(path, mdg_parsed_dir, project)
    elif sequencing_platform == "pacbio":
        generate_pacBio_read_stats(path, mdg_parsed_dir, project)


def generate_pacBio_read_stats(path, mdg_parsed_dir, project):
    # Find HiFi Stats
    splitIndex = None

    with open (f"{path}/PA_Logs/run_stats.csv", "r") as input:
        for i, line in enumerate(input):
            if line.find("HiFi Stats:") != -1:
                splitIndex = i

    pf_df = pd.read_csv(f"{path}/PA_Logs/run_stats.csv", header=2 ,skiprows= lambda x: x >= splitIndex)
    hifi_df = pd.read_csv(f"{path}/PA_Logs/run_stats.csv", header = splitIndex + 1)

    # Drop rows for different projects
    pf_df.drop(pf_df[pf_df["Project"] != project].index, inplace=True)
    hifi_df.drop(hifi_df[hifi_df["Project"] != project].index, inplace=True)

    pf_df.to_csv(f"{mdg_parsed_dir}/pf_polymerase_read_stats.csv")
    hifi_df.to_csv(f"{mdg_parsed_dir}/hifi_stats.csv")


def generate_illumina_read_stats(path, mdg_parsed_dir, project):
    # Copy file
    read_stats = shutil.copy(f"{path}/PA_Logs/trimming.logs.csv", mdg_parsed_dir)

    # Remove all lines of csv that don't include project name
    with open(read_stats, "r") as input:
        with open(f"{mdg_parsed_dir}/temp.csv", "w") as output:
            # iterate all lines from file
            for i, line in enumerate(input):
                # find all lines w/ project name and header row
                if project in line.strip("\n") or i == 0:
                    output.write(line)
            
    # Generate cleaned file, remove copy
    clean_data = f"{mdg_parsed_dir}/trimming.logs_clean.csv"
    os.replace(f"{mdg_parsed_dir}/temp.csv", clean_data)
    os.remove(read_stats)

    # Rename columns
    df = pd.read_csv(clean_data)
    df = df.rename(columns={"R1_Bases": "R1 Bases", "R1_Reads": "R1 Reads", "R1_RL": "R1 RL", "R1_TrimBases": "R1 Trim Bases", "R1_TrimReads": "R1 Trim Reads", "R1_TrimRL": "R1 Trim RL", "R2_Bases": "R2 Bases", "R2_Reads": "R2 Reads", "R2_RL": "R2 RL", "R2_TrimBases": "R2 Trim Bases",
                            "R2_TrimReads": "R2 Trim Reads", "R2_TrimRL": "R2 Trim RL"})

    df.to_csv(clean_data)


def copy_pacBio_plots(path, mdg_parsed_dir, cell):
    """
    Copies one cell's plots, prefixed with the cell name so MultiQC can tell cells apart
    """
    image_dir = f"{path}/{cell}/Filter/dataset-reports/reports"
    images_to_copy = [
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
    try:
        os.makedirs(f"{mdg_parsed_dir}/images")
        print(f"Nested directories '{mdg_parsed_dir}/images' created successfully.")
    except FileExistsError:
        print(f"Existing directory '{mdg_parsed_dir}/images' overwritten.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{mdg_parsed_dir}/images'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    for image in images_to_copy:
        try:
            shutil.copy(f"{image_dir}/{image}", f"{mdg_parsed_dir}/images/{cell}_{image}")
        except FileNotFoundError:
            print(f"{image} not found")


def copy_illumina_plots(path, runinfo, mdg_parsed_dir, project):
    # Get list of all images to copy for this project's lanes
    image_dir = f"{path}/PA_Logs/images/"
    images_in_dir = os.listdir(image_dir)
    images_to_copy = []

    for lane in get_lanes(runinfo, project):
        for image in images_in_dir:
            if image.find(f"_{lane}_") != -1:
                images_to_copy.append(image)
    
    # Create images directory
    try:
        os.makedirs(f"{mdg_parsed_dir}/images")
        print(f"Nested directories '{mdg_parsed_dir}/images' created successfully.")
    except FileExistsError:
        print(f"Existing directory '{mdg_parsed_dir}/images' overwritten.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{mdg_parsed_dir}/images'.")
    except Exception as e:
        print(f"An error occurred: {e}")
    for image in images_to_copy:
        shutil.copy(f"{path}/PA_Logs/images/{image}", f"{mdg_parsed_dir}/images")


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
        # Get trimmomatic files
        try:
            trimmomatic_dir = f"{path}/{folder}/TMP/adaptor_trimming_logs"
            for file in os.listdir(trimmomatic_dir):
                if file.endswith(".adap_trim.err"):
                    shutil.copy(f"{path}/{folder}/TMP/adaptor_trimming_logs/{file}", mdg_parsed_dir)

            #shutil.copy(f"{path}/{folder}/TMP/adaptor_trimming_logs/IL34481-001_S50_L005.trimmomatic.log", mdg_parsed_dir)
        except Exception as e:
            print(e)
        # Get cutadapt files
        try:
            cutadapt_dir = f"{path}/{folder}/TMP/adaptor_trimming_logs"
            for file in os.listdir(cutadapt_dir):
                if file.endswith(".cutadapt.log"):
                    shutil.copy(f"{path}/{folder}/TMP/adaptor_trimming_logs/{file}", f"{mdg_parsed_dir}/{file}.json")
        except Exception as e:
            print(e)
        # Get bcl2fastq files. MultiQC's bcl2fastq module keys lanes by RunId, so two Unaligned
        # folders covering the same lane would overwrite each other; tag RunId with the folder name.
        try:
            with open(f"{path}/{folder}/Stats/Stats.json") as fin:
                stats = json.load(fin)
            stats["RunId"] = f"{stats['RunId']} ({folder})"
            with open(f"{mdg_parsed_dir}/Stats.json", "w") as fout:
                json.dump(stats, fout, indent=2)
        except Exception as e:
            print(e)
    elif sequencing_platform == "pacbio":
        try:
            generated = False
            try:
                os.makedirs(f"{mdg_parsed_dir}/isoseq")
                print(f"Nested directories '{mdg_parsed_dir}/isoseq' created successfully.")
                generated = True
            except FileExistsError:
                print(f"Existing directory '{mdg_parsed_dir}/isoseq' overwritten.")
            except PermissionError:
                print(f"Permission denied: Unable to create '{mdg_parsed_dir}/isoseq'.")
            except Exception as e:
                print(f"An error occurred: {e}")

            if generated:
                report_path = f"{path}/{folder}/hifi_reads/read_seg_iso_seq/cromwell_out/outputs"
                reports = os.listdir(report_path)
                for report in reports:
                    if report.endswith(".csv") or report.endswith(".json"):
                        shutil.copy(f"{report_path}/{report}", f"{mdg_parsed_dir}/isoseq")
        except Exception as e:
            print(e)


def copy_assets(outdir):
    try:
        os.makedirs(f"{outdir}/assets")
        print(f"Nested directories '{outdir}/assets' created successfully.")
    except FileExistsError:
        print(f"Existing directory '{outdir}/assets' overwritten.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{outdir}/assets'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    assets_dir = f"{os.path.dirname(os.path.realpath(__file__))}/assets"
    files = os.listdir(assets_dir)
    for fname in files:
        shutil.copy(f"{assets_dir}/{fname}", f"{outdir}/assets")


def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i
    raise ValueError(f"{substring} not found")


def post_process_csv(mdg_parsed_dir):
    """
    Tidies every .csv and .tsv table in the folder: drops duplicate rows and empty columns,
    and turns the pandas index column into a 1-based "No." column
    """
    for x in os.listdir(mdg_parsed_dir):
        if x.endswith((".csv", ".tsv")):
            read_args = TSV_ARGS if x.endswith(".tsv") else {}
            df = pd.read_csv(f"{mdg_parsed_dir}/{x}", **read_args)

            # drop all duplicates
            df = df.drop_duplicates()

            # drop all empty columns
            df = df.dropna(axis=1, how="all")
            df = df.fillna("")

            # Next, modify index column
            if df.columns[0] == "Unnamed: 0":
                df = df.rename(columns={"Unnamed: 0": "No."})
                offset = 1 - df["No."][0]
                df["No."] += offset
            
            # Write results
            df.to_csv(f"{mdg_parsed_dir}/{x}", index=False, **read_args)