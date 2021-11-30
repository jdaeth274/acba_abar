'''
Script to extract identified ABAR elements from seqeuences
'''

import pandas
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

def get_file_paths(hit_csv, file_paths):
    # Function to load up the ABAR hit csv and
    # then append the file paths to these isolates to the csv

    with open(file_paths, "r") as paths:
        dna_locs = paths.read().splitlines()

    hits_csv = pandas.read_csv(hit_csv, header=0)
    id_locs = []
    for index,row in hits_csv.iterrows():
        print(index)
        current_row = hits_csv.iloc[index]
        current_id = re.sub("fna","dna",current_row['id'])
        print(current_id)
        res = [i for i in dna_locs if current_id in i]
        res = str(res[0])
        id_locs.append(res)

    hits_csv['dna_loc'] = id_locs

    return hits_csv

def extract_abars(hit_csv, out_df):
    # Function to take a hit_csv with the dna_loc column
    # and extract the Abar seqeuences and write to a fasta
    # file

    for index, row in hit_csv.iterrows():
        current_row = hit_csv.iloc[index]
        current_fasta = current_row['dna_loc']

        with open(current_fasta, "r") as fasta_seq:
            fasta = SeqIO.read(fasta_seq, "fasta")
            if current_row['ori'] == "forward":
                abar = fasta.seq[(current_row['hit_start'] - 1):current_row['hit_end']]
            else:
                abar = fasta.seq[(current_row['hit_end'] - 1):current_row['hit_start']].reverse_complement()

            abar_record = SeqRecord(fasta.id + "_abar_seq_" + str(current_row['start']) + str(current_row['end']))
            abar.append(abar_record)
            abar_loc = out_df + str(current_row['id']) + "_ABAR.fasta"
            SeqIO.write(abar, abar_loc, "fasta")



