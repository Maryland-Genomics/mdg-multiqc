import argparse
import logging
import subprocess
import os
import MultiQC_Config_Generator
import MDG_Run_Parser
import re

# TODO: move necessary module files to directory, make config generate in run folder
# This script calls both Illumina_files_parser and MultiQC_Config_Generator
def main():
    parser = argparse.ArgumentParser(description="Run MultiQC and generate MDG tables/images")
    parser.add_argument('-i', '--input', metavar='', required=True, help='Path to run folder')
    parser.add_argument('-p', '--project', metavar='', default='', required=True, help="Project which MultiQC will interpret")
    parser.add_argument('-nr', '--noregen', action='store_true', help="Don't Regenerate Illumina Files and Config")
    parser.add_argument('-nm', '--nomultiqc', action='store_true', help="Don't Run MultiQC")
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        raise NotADirectoryError(args.input)

    path = os.path.split(args.input)[1]
    if path == "":
        path = os.path.split(args.input)[0]
        path = os.path.split(path)[1] # Get the first folder up if blank

    runtype = "Illumina"
    if path.find("R84050") != -1:
        runtype = "PacBio"

    unaligned_folder = None

    if runtype == "Illumina":
        unaligned_pattern = re.compile("Unaligned*")
        for filename in os.listdir(path):
            if re.search(unaligned_pattern, filename):
                unaligned_folder = filename
                break
        if unaligned_folder == None:
            raise FileNotFoundError("Unaligned folder was expected, but not found")

    if not args.noregen:
        MDG_Run_Parser.parse(path, args.project, runtype, unaligned_folder)
        MultiQC_Config_Generator.generate_config(path, args.project, runtype, unaligned_folder)
    if not args.nomultiqc:
        if runtype == "Illumina":
            subprocess.run(f"multiqc {path}/Unaligned_8bp/QA/Project_{args.project} --config multiqc_config_{args.project}.yaml --outdir {path}/PA_Logs/multiqc -n multiqc_{args.project} --flat --export --force", shell=True)
        elif runtype == "PacBio":
            subprocess.run(f"multiqc {path}/Project_{args.project}_MDG_Run_Parsed_Files --config multiqc_config_{args.project}.yaml --outdir {path}/PA_Logs/multiqc -n multiqc_{args.project} --flat --export --force", shell=True)

if __name__ == '__main__':
    main()