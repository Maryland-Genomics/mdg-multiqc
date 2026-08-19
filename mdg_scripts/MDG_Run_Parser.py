import pandas as pd
import os
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import re

def parse(p_path : str, p_project : str, p_runtype : str, p_unaliged_folder : str, p_outdir : str):
    project = p_project
    path = p_path
    runtype = p_runtype
    unaligned_folder = p_unaliged_folder
    

    mdg_illumina_folder_name = "MDG_Run_Parsed_Files"
    project_dir = f"{path}/{unaligned_folder}/QA/Project_{project}"
    mdg_parsed_dir = None

    if runtype == "Illumina":
        mdg_parsed_dir =  f"{project_dir}/{mdg_illumina_folder_name}"
    elif runtype == "PacBio":
        mdg_parsed_dir = f"{path}/Project_{project}_{mdg_illumina_folder_name}"


    generate_MDG_parsed_folder(mdg_parsed_dir)
    parse_work_order_info(path, project, mdg_parsed_dir)
    generate_sample_info(path, mdg_parsed_dir, project, runtype)
    generate_read_stats(path, mdg_parsed_dir, project, runtype)
    copy_plots(path, mdg_parsed_dir, project, runtype)
    parse_megablastqc_info(path, mdg_parsed_dir, project, runtype, unaligned_folder)
    move_module_files(path, project, mdg_parsed_dir, unaligned_folder, runtype)
    copy_assets(p_outdir)

    if (runtype == "Illumina"):
        parse_sequencing_info(path, mdg_parsed_dir)
        parse_flowcell_info(path, mdg_parsed_dir)
    elif (runtype == "PacBio"):
        pass
    post_process_csv(mdg_parsed_dir)


def generate_MDG_parsed_folder(mdg_parsed_dir):
    """
    Creates a folder in the project directory to store generated data files
    """
    try:
        os.makedirs(mdg_parsed_dir)
        print(f"Nested directories '{mdg_parsed_dir}' created successfully.")
    except FileExistsError:
        print(f"Existing directory '{mdg_parsed_dir}' overwritten.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{mdg_parsed_dir}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


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


def parse_megablastqc_info(path, mdg_parsed_dir, project, runtype, unaligned_folder):
    """
    Creates a file to store parsed MegaBlast QC info in the illumina folder
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
    if (runtype == "Illumina"):
        tophits_directory = f"{path}/{unaligned_folder}/QC/Project_{project}"
    elif (runtype == "PacBio"):
        tophits_directory = f"{path}/{get_cell(path, project)}/hifi_reads/QC"

    tophits_count = 0

    try:
        if runtype == "Illumina":
            for sample in os.listdir(tophits_directory):
                if sample.startswith("Sample"):
                    # Assume file called *.topHits.txt, get all files
                    for tophits_file in os.listdir(f"{tophits_directory}/{sample}/QC"):
                        # Get all files in directory that end with topHits.txt
                        if tophits_file.endswith("topHits.txt"):
                            data, tophits_count = get_tophits(data, f"{tophits_directory}/{sample}/QC", tophits_file, mdg_parsed_dir, project, tophits_count, runtype)
        elif runtype == "PacBio":
            for tophits_file in os.listdir(tophits_directory):
                if tophits_file.endswith("topHits.txt"):
                    data, tophits_count = get_tophits(data, tophits_directory, tophits_file, mdg_parsed_dir, project, tophits_count, runtype)

        # Write data to new csv file
        df = pd.DataFrame(data)
        df.to_csv(f"{mdg_parsed_dir}/all_tophits.csv", index=True)
    except Exception as e:
        print(e)


def get_tophits(data, tophits_directory, tophits_file, mdg_parsed_dir, project, tophits_count, runtype):
    try:
        with open(f"{tophits_directory}/{tophits_file}") as fin:
            tophits = fin.read().splitlines(True)
    
        # Parse txt contents
        if runtype == "Illumina":
            id_index_start = tophits[0].find(":", tophits[0].find(":") + 1) + 1
            id_index_end = tophits[0].find("lane:") + 6
            id = tophits[0][id_index_start:id_index_end]
            id = id.replace(":lane:", "-")
            data["ID"].append(id)
        elif runtype == "PacBio":
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


def generate_sample_info(path, mdg_parsed_dir, project, runtype):
    """
    Copies runinfo.csv into the illumina folder and trims file to be read
    """
    # Copy file
    sample_info = shutil.copy(f"{path}/PA_Logs/runinfo.csv", mdg_parsed_dir)

    # Create new csv file with every line after first (keep copy intact)
    with open(sample_info, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(f"{mdg_parsed_dir}/temp.csv", 'w') as fout:
        fout.writelines(data[1:])

    # Rename cleaned file, remove copy
    clean_data = f"{mdg_parsed_dir}/runinfo_clean.csv"
    os.replace(f"{mdg_parsed_dir}/temp.csv", clean_data)
    os.remove(sample_info)

    # Processing
    df = pd.read_csv(clean_data)
    # filter for rows in the selected project
    df.drop(df[df["project"] != project].index, inplace=True)

    if runtype == "Illumina":
        # Drop some columns
        df = df.drop(labels=["sequencing_sub_sample_id", "project", "mp_insert_size", ], axis=1)
        # Rename columns
        df = df.rename(columns={"lane": "Lane", "master_sample_id": "Sample", "library_id": "Library", "pool_name": "Pool Name", "index_num": "Library Index", "external_id": "External ID", "strain": "Strain", "type": "Library Type", "median_smear_size": "Size (bp)",
                                "species": "Species", "sample_name": "Sample Name", "method": "Method", "loaded": "Loaded", "needed": "Needed", "index_type": "Library Index Type", "index_sequence": "Library Index Sequence", "premade_pool": "Premade Pool",
                                "pool_index_set": "Pool Index", "pool_index_1": "Pool Index 1", "pool_index_2": "Pool Index 2", "library_prep_kit": "Library Prep Kit", "sample_index_name": "Sample Index", "sample_index": "Sample Index Sequence",
                                "read_pairs": "Read Pairs", "control_projects": "Control Projects", "lib_category": "Library Category", "primer": "Primer", "comments": "Comments"})
    elif runtype == "PacBio":
        df = df.drop(labels=["project"], axis=1)

    df.to_csv(clean_data)


def generate_read_stats(path, mdg_parsed_dir, project, runtype):
    """
    Copies trimming.logs.csv into the illumina folder and trims file to be read
    """
    if runtype == "Illumina":
        generate_illumina_read_stats(path, mdg_parsed_dir, project)
    elif runtype == "PacBio":
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


def copy_plots(path, mdg_parsed_dir, project, runtype):
    """
    Copies plots to illumina folder
    """
    if (runtype == "Illumina"):
        copy_illumina_plots(path, mdg_parsed_dir, project)
    elif(runtype == "PacBio"):
        copy_pacBio_plots(path, mdg_parsed_dir, project)


def copy_pacBio_plots(path, mdg_parsed_dir, project):
    cell = get_cell(path, project)
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
            shutil.copy(f"{image_dir}/{image}", f"{mdg_parsed_dir}/images")
        except FileNotFoundError:
            print(f"{image} not found")


def copy_illumina_plots(path, mdg_parsed_dir, project):
    # Get runinfo.csv
    runinfo_df = pd.read_csv(f"{path}/PA_Logs/runinfo.csv", header=1)
    runinfo_df.drop(runinfo_df[runinfo_df["project"] != project].index, inplace=True)

    lanes = runinfo_df["lane"]
    unique_lanes = np.unique(lanes)

    # Get list of all images to copy
    image_dir = f"{path}/PA_Logs/images/"
    images_in_dir = os.listdir(image_dir)
    images_to_copy = []

    for lane in unique_lanes:
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


def parse_work_order_info(path, project, mdg_parsed_dir):
    try:
        df = pd.read_csv(f"{path}/PA_Logs/runinfo.tmp")
    except(FileNotFoundError):
        df = pd.read_csv(f"{path}/PA_Logs/runinfo.tmp.org")

    df = df[["Work Order", "WO PI", "WO Submitter", "WO Additional Contact(s)", "WO Billing Contact", "Quote Number"]]
    df = df[df["Work Order"].str.find(project) != -1]
    df = df.dropna(axis=0, how="all")

    sample_df = pd.read_csv(f"{path}/PA_Logs/runinfo.csv", header=1)
    try:
        unique_lib = len(pd.unique(sample_df["library"]))
    except:
        unique_lib = len(pd.unique(sample_df["library_id"]))

    unique_lib_arr = [unique_lib] * df.shape[0]
    
    df.insert(6, "Sample Count", unique_lib_arr, True)
    
    work_order_info = f"{mdg_parsed_dir}/work_order_info.csv"
    df.to_csv(work_order_info, index=False)


def move_module_files(path, project, mdg_parsed_dir, unaligned_folder, runtype):
    if runtype == "Illumina":
        # Get trimmomatic files
        try:
            shutil.copy(f"{path}/{unaligned_folder}/TMP/adaptor_trimming_logs/sample.adap_trim.err", mdg_parsed_dir)
        except Exception as e:
            print(e)
        # Get cutadapt files
        try:
            cutadapt_dir = f"{path}/{unaligned_folder}/TMP/adaptor_trimming_logs"
            for file in os.listdir(cutadapt_dir):
                if file.endswith(".cutadapt.log"):
                    shutil.copy(f"{path}/{unaligned_folder}/TMP/adaptor_trimming_logs/{file}", f"{mdg_parsed_dir}/{file}.json")
        except Exception as e:
            print(e)
        # Get bcl2fastq files
        try:
            shutil.copy(f"{path}/{unaligned_folder}/Stats/Stats.json", mdg_parsed_dir)
        except Exception as e:
            print(e)
    elif runtype == "PacBio":
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
                report_path = f"{path}/{get_cell(path, project)}/hifi_reads/read_seg_iso_seq/cromwell_out/outputs"
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


def get_cell(path, project):
    df = pd.read_csv(f"{path}/PA_Logs/runinfo.csv", header=1)
    line = df[df["project"] == project].index[0]
    cell = df.loc[line, "cell"]
    return f"1_{cell}"


def post_process_csv(mdg_parsed_dir):
    for x in os.listdir(mdg_parsed_dir):
        if x.endswith(".csv"):
            df = pd.read_csv(f"{mdg_parsed_dir}/{x}")

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
            df.to_csv(f"{mdg_parsed_dir}/{x}", index=False)