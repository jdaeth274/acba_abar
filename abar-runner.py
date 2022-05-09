# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import argparse
import sys
import time

from python.common import main

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def parse_args():
    purpose = ''' This is a scipt to find the occurence of AbaR elements in a collection of fasta sequences Usage:
        python acba_mlst_runner.py --seqs <list_of_fastas>  --output <output_prefixes> --threads <num_cores_to_use>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='acba_mlst_runner.py')

    parser.add_argument('--seqs', required=True, help='List of seqeuence files (FASTA) (required)', type=str)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--dna-dir', default=None, help= 'Location of directory of DNA files with output.mfa file present',
                        type=str)
    parser.add_argument('--contig-bounds', default=None, help="Location of the contig bounds of files", type=str)
    parser.add_argument('--threads', default=1, help='Number of threads to use for ORF finder', type=int)
    parser.add_argument('--no-contigs', default=True, action='store_false', help='Do not use contig bounds when defining a hit [default = True so use contigs]')
    parser.add_argument('--comM', default=False, action='store_true', help='Find out if comM is disrupted')
    parser.add_argument('--abar-extract', default=False, action='store_true', help='extract fasta seqs of identified abars',
                        dest='abar')
    parser.add_argument('--presence-only', default=False, action='store_true',
                        help="Basing AbaR presence only on if there is a left end and a right end in a sample, no contig or position cutoff further",
                        dest='both_present')
    parser.add_argument('--no-blast',default=True, action='store_false',
                        help="Don't run the BLAST search (useful if just want to extract abars)",
                        dest='blast')


    args = parser.parse_args()

    return args
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    start = time.perf_counter()
    if main(parse_args()):
        end = time.perf_counter()
        print("Took this long to complete: %s (s)" % (end - start))
    sys.exit("Done")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
