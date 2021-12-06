#!/bin/bash

set -e
START=$SECONDS
## File to take in a list of fasta files and then run the mcl
## pipeline on them 

if [ $# -lt 3 ]
then
  echo "MCL use pipelin. Usage:
  bash ./isa_pipeline.sh <list of fastas> <file_prefix> <MCL inflation value> <threads>"
  exit
fi

FASTA_FILE=$1
PREFIX=$2
INFLATION=$3
if [ -z $4 ]
then
    THREADS=$(grep -c ^processor /proc/cpuinfo)
    THREADS=$(( THREADS - 1 ))
else
    THREADS=$4
fi
echo "This is the fasta file $FASTA_FILE"
echo "This is the prefix $PREFIX"
echo "This is the inflations value $INFLATION"
echo "This is the number of threads $THREADS"

## First step to create the mfa file 

printf "\n"
printf "Creating .mfa file"
MFA_START=$SECONDS
while read line <&3
do
cat $line >> "${PREFIX}.mfa"
done 3< $FASTA_FILE
MFA_END=$(( SECONDS - MFA_START ))
printf "\rCreating .mfa file $MFA_END (s) \n"
# Now to make the blast db file 
printf "Making BLAST DB"
DB_START=$SECONDS
makeblastdb -dbtype nucl -in "${PREFIX}.mfa" -out "${PREFIX}_db" > /dev/null 2>&1
DB_END=$(( SECONDS - DB_START ))
printf "\rMaking BLAST DB in $DB_END (s) \n"
# Now for the blast run
printf "Running BLAST"
BLAST_START=$SECONDS
blastn -db "${PREFIX}_db" -query "${PREFIX}.mfa" -evalue 1 \
-outfmt 6 -out "${PREFIX}_all_v_all.tsv" -num_alignments 1000000 \
-num_threads $THREADS \
> /dev/null 2>&1
BLAST_END=$(( SECONDS - BLAST_START ))
printf "\rRunning BLAST in $BLAST_END (s) \n"
# Now deblasting the runs 
echo "Now on the MCXdeblast"
mcxdeblast --m9 --out "${PREFIX}_mcl_format.abc" --ecut 1 --line-mode abc \
"${PREFIX}_all_v_all.tsv" > /dev/null 2>&1

#Now loading up the mcl runs
echo "Now on the mcxload"
mcxload -abc "${PREFIX}_mcl_format.abc" --stream-mirror --stream-neg-log10 \
-stream-tf 'ceil(200)' -o "${PREFIX}_mcl.mci" \
-write-tab "${PREFIX}_mcl.tab"

#Now for the mcl runs 
echo ""
echo "Now on MCL"
mcl "${PREFIX}_mcl.mci" -I $INFLATION -use-tab "${PREFIX}_mcl.tab" \
-o "${PREFIX}_mcl_res_I_${INFLATION}"

# Now finally to create the edge list file for cytoscape
echo ""
echo "Creating the tab file"
perl -lane 'while ($#F > 0) {my $first = shift @F; foreach my $second (@F) {print $first."\t".$second;}}' \
 "${PREFIX}_mcl_res_I_${INFLATION}" > "${PREFIX}_mcl_res_I_${INFLATION}.tab"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
END=$(( SECONDS - START ))
echo "Finished in $END (s)"

