import os
from python import pre_blast_process
from python import blast_runs
import re

def main(input_args):
    """ Main function to run through steps of ABar and comM identification
    """

    seq_files = open(input_args.seqs, "r")
    seq_lines = seq_files.read().splitlines()

    if input_args.dna_dir is None:
        pre_blast_process.make_dna_db(seq_lines)

    if input_args.contig_bounds is None:
        pre_blast_process.contig_bounds(seq_lines)

    ## BLAST the fragments now
    python_dir_name = os.path.dirname(os.path.realpath(__file__))
    data_dir = re.sub("python", "data", python_dir_name)
    tmp_dna_dir = "./tmp_dna_lib"

    blast_runs.blast_runs(tmp_dna_dir, data_dir)

    