import subprocess
import os
import pandas

def merge_blast_hits(contig_dir, R_dir):
    """ Function to run the R scripts for merging the conserved region hits
        directly after the BLAST search"""

    merge_cmd = "Rscript --vanilla " + R_dir + "/testing_out_functions.R left_end_blast.csv right_end_blast.csv " + contig_dir

    try:
        subprocess.check_output(merge_cmd, shell=True)
    except:
        subprocess.SubprocessError

def extract_hits(out_name, contig_use):
    """ Function to take the merged blast dbs and extracate the hits
        using the criteria that a left and right end must be on the
        same contig and witin 150k of each other on that contig"""

    left_end_csv = pandas.read_csv(filepath_or_buffer="./left_end_merged.csv",header=0)
    right_end_csv = pandas.read_csv(filepath_or_buffer="./right_end_merged.csv",header=0)

    ## run through unique values in left end (if they're not in left end they're not
    ## a hit)

    left_end_vals = left_end_csv.subject.unique()
    id = []
    hit_start = []
    hit_end = []
    ori = []
    contig = []

    for subject in left_end_vals:
        current_left_row = left_end_csv[left_end_csv['subject'] == subject]
        current_right_row = right_end_csv[right_end_csv['subject'] == subject]
        if current_right_row.empty:
            continue
        for index, row in current_left_row.iterrows():
            if contig_use:
                right_match = current_right_row[(current_right_row['ori'] == current_left_row.loc[index,'ori']) & \
                                                (current_right_row['contig'] == row.loc['contig'])]
                if right_match.empty:
                    continue
                if row.loc['ori'] == "forward":
                    right_hits = right_match[right_match['sstart'] <= (row.loc['send'] + 150000)]
                    if right_hits.empty:
                        continue
                    right_hits = right_hits.sort_values(by='sstart', ascending=True)
                    right_hits = right_hits.reset_index(drop = True)
                elif row.loc['ori'] == "reverse":
                    right_hits = right_match[right_match['sstart'] >= (row.loc['send'] - 150000)]
                    if right_hits.empty:
                        continue
                    right_hits =  right_hits.sort_values(by='sstart', ascending=False)
                    right_hits = right_hits.reset_index(drop=True)

                id.append(row.loc['subject'])
                hit_start.append(row.loc['sstart'])
                hit_end.append(right_hits.loc[0,"send"])
                ori.append(row.loc['ori'])
                contig.append(row.loc['contig'])
            else:
                right_match = current_right_row[(current_right_row['ori'] == current_left_row.loc[index, 'ori'])]
                if right_match.empty:
                    continue
                if row.loc['ori'] == "forward":
                    right_hits = right_match[right_match['sstart'] <= (row.loc['send'] + 150000)]
                    if right_hits.empty:
                        continue
                    right_hits = right_hits.sort_values(by='sstart', ascending=True)
                    right_hits = right_hits.reset_index(drop=True)
                elif row.loc['ori'] == "reverse":
                    right_hits = right_match[right_match['sstart'] >= (row.loc['send'] - 150000)]
                    if right_hits.empty:
                        continue
                    right_hits = right_hits.sort_values(by='sstart', ascending=False)
                    right_hits = right_hits.reset_index(drop=True)

                id.append(row.loc['subject'])
                hit_start.append(row.loc['sstart'])
                hit_end.append(right_hits.loc[0, "send"])
                ori.append(row.loc['ori'])
                contig.append(row.loc['contig'])


    out_data = {
        'id':id,
        'hit_start':hit_start,
        'hit_end':hit_end,
        'ori':ori,
        'contig':contig
    }

    out_df = pandas.DataFrame(out_data)
    out_df.to_csv(out_name, index=False)








