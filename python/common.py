import os
from python import pre_blast_process
from python import blast_runs
from python import post_blast_process
from python import abar_extraction
from python import comM_extraction
import re
import subprocess
import sys
import shutil


def main(input_args):
    """ Main function to run through steps of ABar and comM identification
    """
    print(input_args)
    with open(input_args.seqs, "r") as input_handle:
        seq_lines = input_handle.read().splitlines()

    if input_args.dna_dir is None:
        dna_dir = pre_blast_process.make_dna_db(seq_lines)
    else:
        dna_dir = input_args.dna_dir

    if input_args.contig_bounds is None:
        contig_bounds = pre_blast_process.contig_bounds(seq_lines)
    else:
        contig_bounds = input_args.contig_bounds

    ## BLAST the fragments now
    if input_args.blast:
        python_dir_name = os.path.dirname(os.path.realpath(__file__))
        data_dir = re.sub("python", "data", python_dir_name)
        #tmp_dna_dir = "./tmp_dna_lib"

        blast_runs.blast_runs(dna_dir, data_dir, input_args.comM)

        ## Post-BLAST merging and hit idents
        R_dir = re.sub("python", "R", python_dir_name)
        post_blast_process.merge_blast_hits(contig_dir=contig_bounds,R_dir=R_dir)

        out_name = input_args.output + "_hits.csv"
        post_blast_process.extract_hits(out_name, input_args.no_contigs)

    ## Check if want to search for comM
    if input_args.comM:
        out_name = input_args.output + "_complete_comM.csv"
        post_blast_process.extract_comM(out_name)
        dna_cmd = "cd " + dna_dir + " && ls -d $PWD/*.dna > tot_dna_files.txt"
        dna_list = dna_dir + "/tot_dna_files.txt"
        try:
            subprocess.check_output(dna_cmd, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed making the dna list of isolates for comM extraction")

        comM_dir = input_args.output + "_comM_gap_seqs"
        if os.path.isdir(comM_dir):
            shutil.rmtree(comM_dir)
            os.mkdir(comM_dir)
        else:
            os.mkdir(comM_dir)
        R_dir = re.sub("python", "R", python_dir_name)
        comM_csv = "comM_hits.csv"
        comM_extraction.get_hit_csv(comM_csv, R_dir, contig_dir=contig_bounds)
        hit_csv = comM_extraction.get_file_paths("comM_gaps.csv", dna_list)
        comM_dir = "comM_gaps_seqs"
        comM_extraction.extract_comM(hit_csv, comM_dir)

    if input_args.abar:
        dna_cmd = "cd " + dna_dir + " && ls -d $PWD/*.dna > tot_dna_files.txt"
        dna_list = dna_dir + "/tot_dna_files.txt"
        try:
            subprocess.check_output(dna_cmd, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed making the dna list of isolates for abar extraction")

        abar_dir = input_args.output + "_abar_seqs"
        if os.path.isdir(abar_dir):
            shutil.rmtree(abar_dir)
            os.mkdir(abar_dir)
        else:
            os.mkdir(abar_dir)

        out_name = input_args.output + "_hits.csv"
        hit_csv = abar_extraction.get_file_paths(out_name, dna_list)
        abar_extraction.extract_abars(hit_csv, abar_dir)




    return True

