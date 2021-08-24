
from pre_blast_process import make_dna_db
from pre_blast_process import contig_bounds

def main(input_args):
    """ Main function to run through steps of ABar and comM identification
    """

    seqs = input_args.seqs
    if input_args.dna_dir is None:
        make_dna_db(seqs)

    if input_args.contig_bounds is None:
        contig_bounds(seqs)



