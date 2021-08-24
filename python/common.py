
from python import pre_blast_process


def main(input_args):
    """ Main function to run through steps of ABar and comM identification
    """

    seqs = input_args.seqs
    if input_args.dna_dir is None:
        pre_blast_process.make_dna_db(seqs)

    if input_args.contig_bounds is None:
        pre_blast_process.contig_bounds(seqs)



