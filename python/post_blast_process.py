import subprocess
import os

def merge_blast_hits(contig_dir, R_dir):
    """ Function to run the R scripts for merging the conserved region hits
        directly after the BLAST search"""

    merge_cmd = "Rscript --vanilla " + R_dir + "/testing_out_functions.R left_end_blast.csv right_end_blast.csv " + contig_dir

    try:
        subprocess.check_output(merge_cmd, shell=True)
    except:
        subprocess.SubprocessError


