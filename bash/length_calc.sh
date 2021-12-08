#! /bin/bash

set -e 

if [ $# -lt 2 ]
then
  echo "length_calc.sh Usage:
  bash ./length_calc.sh <list of fastas> <file_prefix> "
  exit
fi


LIST_O_FILES=$1
PREFIX=$2

OUT_CSV="${PREFIX}.csv"

if [ -f $OUT_CSV ]
then
    printf "\n Removing out csv with same name \n"
    rm $OUT_CSV
    printf "id,length\n" > $OUT_CSV
else
    printf "id,length\n" > $OUT_CSV
fi

NUM_LINES=$(wc -l $LIST_O_FILES)
printf "Beginning read counting on %s files \n \n" $NUM_LINES
COUNTER=0
while read line <&3
do
    LENGTH=$(awk '/^>/ {next} {counter+=length($0)} END {print counter}' $line)
    id=$(basename $line | sed 's/\..*$//g')
    printf "%s,%s\n" $id $LENGTH >> $OUT_CSV
    COUNTER=$(( COUNTER + 1 ))
    printf "\r Finished on file: %s" $COUNTER
done 3< $LIST_O_FILES
printf "\n Done \n"
