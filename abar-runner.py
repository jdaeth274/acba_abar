# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import argparse
import sys

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

    args = parser.parse_args()

    return args
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main(parse_args())
    sys.exit("Done")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
