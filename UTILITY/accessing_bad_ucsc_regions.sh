#!/bin/bash

CUT_OFF=$1
((CUT_OFF_Y=$CUT_OFF/2))

MERGE_BAD_COVERAGE()
{
awk 'NR>1 gsub (/,/,"\t",$0) gsub (/:/,"\t",$0) {print FILENAME,$1,$2-1,$2,$5}' \
$SAMPLE \
| awk 'BEGIN {OFS="\t"} $5<'"$CUT_OFF"'&&$2~/^[0-9]/ {split($1,SM_TAG,"/"); print $2,$3,$4,$5,SM_TAG[7]} \
$5<'"$CUT_OFF"'&&$2=="X" {split($1,SM_TAG,"/"); print $2,$3,$4,$5,SM_TAG[7]} \
$5<'"$CUT_OFF_Y"'&&$2=="Y" {split($1,SM_TAG,"/"); print $2,$3,$4,$5,SM_TAG[7]}' \
| bedtools merge -i - -c 5 -o distinct \
>> bad_places_$CUT_OFF.txt
}

for SAMPLE in $(ls /isilon/sequencing/Seq_Proj/CGC*CGC*CGC/*/*/REPORTS/DEPTH_OF_COVERAGE/UCSC_CODING_PLUS_10bp/*ALL_UCSC_CODING_10bpFlanks.EveryBase.csv) ;
do
MERGE_BAD_COVERAGE
done

sort -k 1,1 -k2,2n bad_places_$CUT_OFF.txt \
| bedtools merge -i - -c 4,4 -o count_distinct,distinct -delim ";" \
| awk 'BEGIN {print "CHR","START","END","LENGTH","SAMPLE_COUNT_BELOW_""'$CUT_OFF'","SAMPLE_LIST"} {print $1,$2+1,$3,$3-$2,$4,$5}' \
| sed 's/ /,/g' \
>| "ucsc_coding_exons_plus10bp_below_"$CUT_OFF"x_and"$CUT_OFF_Y"x_for_Y_summary.csv"

## in merge bad coverage, can get a list of all depths in bedtools merge if by doing "bedtools merge -i - -c 5,4 -o distinct,collapse"
