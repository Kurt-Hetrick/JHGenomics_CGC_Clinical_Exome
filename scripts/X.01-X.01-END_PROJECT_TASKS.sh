# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q rnd.q,prod.q,test.q

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on
# redirecting stderr/stdout to file as a log.

set

CORE_PATH=$1
DATAMASH=$2

PROJECT=$3

# Sorting concatenated sample meta by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/SAMPLE_META.txt \
| uniq \
| awk 'BEGIN {print "SM_TAG","RG_ID","RG_PU","LIBRARY","LIBRARY_WELL","FAMILY","FATHER","MOTHER","GENDER","PHENOTYPE"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/SAMPLE_META_HEADER.txt

# Sorting concatenated gender check report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/GENDER_CHECK.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","X_AVG_DP","X_NORM_DP","Y_AVG_DP","Y_NORM_DP"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/GENDER_CHECK_HEADER.TXT

# Sorting concatenated verify bam ID report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/VERFIY_BAM_ID.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","VERIFYBAM_FREEMIX","VERIFYBAM_#SNPS","VERIFYBAM_FREELK1","VERIFYBAM_FREELK0","VERIFYBAM_AVG_DP"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/VERFIY_BAM_ID_HEADER.TXT

# Sorting concatenated insert size report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/INSERT_SIZE_METRICS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","MEDIAN_INSERT_SIZE","MEDIAN_ABSOLUTE_DEVIATION_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/INSERT_SIZE_METRICS_HEADER.TXT

# Sorting concatenated read 1 alignment summary metric report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_1_METRICS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","PF_NOISE_READS_R1","PF_READS_ALIGNED_R1","PCT_PF_READS_ALIGNED_R1","PF_ALIGNED_BASES_R1","PF_HQ_ALIGNED_READS_R1",\
"PF_HQ_ALIGNED_BASES_R1","PF_HQ_ALIGNED_Q20_BASES_R1","PF_HQ_MEDIAN_MISMATCHES_R1","PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1",\
"PCT_READS_ALIGNED_IN_PAIRS_R1","BAD_CYCLES_R1","STRAND_BALANCE_R1","PCT_ADAPTER_R1"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_1_METRICS_HEADER.TXT

# Sorting concatenated read 2 alignment summary metric report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_2_METRICS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","PF_NOISE_READS_R2","PF_READS_ALIGNED_R2","PCT_PF_READS_ALIGNED_R2","PF_ALIGNED_BASES_R2","PF_HQ_ALIGNED_READS_R2",\
"PF_HQ_ALIGNED_BASES_R2","PF_HQ_ALIGNED_Q20_BASES_R2","PF_HQ_MEDIAN_MISMATCHES_R2","PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2",\
"PCT_READS_ALIGNED_IN_PAIRS_R2","BAD_CYCLES_R2","STRAND_BALANCE_R2","PCT_ADAPTER_R2"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_2_METRICS_HEADER.TXT

# Sorting concatenated pair alignment summary metric report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_PAIR_METRICS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","TOTAL_READS","PF_NOISE_READS_PAIR","PF_READS_ALIGNED_PAIR","PCT_PF_READS_ALIGNED_PAIR",\
"PF_ALIGNED_BASES_PAIR","PF_HQ_ALIGNED_READS_PAIR","PF_HQ_ALIGNED_BASES_PAIR","PF_HQ_ALIGNED_Q20_BASES_PAIR",\
"PF_HQ_MEDIAN_MISMATCHES_PAIR","PF_MISMATCH_RATE_PAIR","PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR","MEAN_READ_LENGTH",\
"PCT_READS_ALIGNED_IN_PAIRS_PAIR","BAD_CYCLES_PAIR","STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR","PCT_ADAPTER_PAIR"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_PAIR_METRICS_HEADER.TXT

# Sorting concatenated pair mark duplicates report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/MARK_DUPLICATES_METRICS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","UNPAIRED_READS_EXAMINED","READ_PAIRS_EXAMINED","UNMAPPED_READS","UNPAIRED_READ_DUPLICATES",\
"READ_PAIR_DUPLICATES","READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/MARK_DUPLICATES_METRICS_HEADER.TXT

# Sorting concatenated pair hyb selection report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/HYB_SELECTION.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","BAIT_SET","GENOME_SIZE","BAIT_TERRITORY","TARGET_TERRITORY","PF_UNIQUE_READS","PCT_PF_UQ_READS","PF_UQ_READS_ALIGNED",\
"PCT_PF_UQ_READS_ALIGNED","PF_UQ_BASES_ALIGNED","PF_UQ_GIGS_ALIGNED","PCT_SELECTED_BASES","MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE",\
"ZERO_CVG_TARGETS_PCT","PCT_EXC_DUPE","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET","PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X",\
"PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X",\
"HS_LIBRARY_SIZE","AT_DROPOUT","GC_DROPOUT","HET_SNP_SENSITIVITY","HET_SNP_Q"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/HYB_SELECTION_HEADER.TXT

# Sorting concatenated bait bias report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/BAIT_BIAS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","Cref_Q","Gref_Q"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/BAIT_BIAS_HEADER.TXT

# Sorting concatenated pre adapter report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/PRE_ADAPTER.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","DEAMINATION_Q","OxoG_Q"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/PRE_ADAPTER_HEADER.TXT

# Sorting concatenated pre adapter report by SM TAG and adding headers

sort -k 1 $CORE_PATH/$PROJECT/TEMP/QUALITY_YIELD_METRICS.TXT \
| uniq \
| awk 'BEGIN {print "SM_TAG","PCT_Q20_BASES","PCT_Q30_BASES"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/QUALITY_YIELD_METRICS_HEADER.TXT

######################################################################################################

# Joining all of the files together to make a QC report

TIMESTAMP=`date '+%F.%H-%M-%S'`s

join -i -j 1 $CORE_PATH/$PROJECT/TEMP/SAMPLE_META_HEADER.txt $CORE_PATH/$PROJECT/TEMP/GENDER_CHECK_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/VERFIY_BAM_ID_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/INSERT_SIZE_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_1_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_2_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_PAIR_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/MARK_DUPLICATES_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/HYB_SELECTION_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/BAIT_BIAS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/PRE_ADAPTER_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/QUALITY_YIELD_METRICS_HEADER.TXT \
| sed 's/ /,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".QC_REPORT."$TIMESTAMP".csv"

# Concatenate all aneuploidy reports together

( cat $CORE_PATH/$PROJECT/*/*/REPORTS/ANEUPLOIDY_CHECK/*.chrom_count_report.txt | grep "^SM_TAG" | uniq ; \
cat $CORE_PATH/$PROJECT/*/*/REPORTS/ANEUPLOIDY_CHECK/*.chrom_count_report.txt | grep -v "SM_TAG" ) \
| sed 's/\t/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".ANEUPLOIDY_CHECK."$TIMESTAMP".csv"

# Summarize Wall Clock times

sed 's/,/\t/g' $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv" \
| sort -k 1 -k 2 -k 3 \
| awk 'BEGIN {OFS="\t"} {print $0,($6-$5),($6-$5)/60,($6-$5)/3600}' \
| $DATAMASH/datamash -s -g 1,2 max 7 max 8 max 9 | tee $CORE_PATH/$PROJECT/TEMP/WALL.CLOCK.TIMES.BY.GROUP.txt \
| $DATAMASH/datamash -g 1 sum 3 sum 4 sum 5 \
| awk 'BEGIN {print "SAMPLE_PROJECT","WALL_CLOCK_SECONDS","WALL_CLOCK_MINUTES","WALL_CLOCK_HOURS"} {print $0}' \
| sed -r 's/[[:space:]]+/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.BY_SAMPLE.csv"

sed 's/\t/,/g' $CORE_PATH/$PROJECT/TEMP/WALL.CLOCK.TIMES.BY.GROUP.txt \
| awk 'BEGIN {print "SAMPLE_PROJECT","TASK_GROUP","WALL_CLOCK_SECONDS","WALL_CLOCK_MINUTES","WALL_CLOCK_HOURS"} {print $0}' \
| sed -r 's/[[:space:]]+/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.BY_SAMPLE_GROUP.csv"
