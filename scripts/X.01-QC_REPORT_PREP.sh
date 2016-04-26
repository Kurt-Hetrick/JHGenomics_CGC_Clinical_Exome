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

SAMTOOLS_DIR=$1
CORE_PATH=$2
DATAMASH_DIR=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6
FATHER=$7
MOTHER=$8
GENDER=$9
AFFECTED=${10}

# sleeping for a random number of seconds before executing so datastream don't collide on the same file...hopefully

echo $RANDOM | awk '{print $0/1000*5}' | awk '{split ($0,RANDOM_INTEGER,"."); print "sleep",RANDOM_INTEGER[1]}' | bash

# Grabbing the BAM header (for RG ID,PU,LB,etc)

##### THIS IS THE HEADER, NEED TO THINK ABOUT HOW TO GET THIS BACK IN HERE ##########
## | awk 'BEGIN {print "SM_TAG","RG_ID","RG_PU","Library_Name","Library_Well","FAMILY","FATHER","MOTHER","GENDER","PHENOTYPE"}' \
#####################################################################################
## This is based on our library name ##
## | awk 'BEGIN {OFS="\t"} {split($4,Library,"_"); print $0,Library[1],Library[2],Library[3],Library[4],\
#####################################################################################

$SAMTOOLS_DIR/samtools view -H \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
| grep ^@RG \
| awk '{split($5,SMtag,":"); split($2,ID,":"); split($6,PU,":"); split($3,Library,":"); print SMtag[2]"\t"ID[2]"\t"PU[2]"\t"Library[2]}' \
| $DATAMASH_DIR/datamash -s -g 1 collapse 2 collapse 3 unique 4 | sed 's/,/;/g' \
| awk 'BEGIN {OFS="\t"} {split($4,Library,"_"); print $0,Library[1],"'$FAMILY'","'$FATHER'","'$MOTHER'","'$GENDER'","'$AFFECTED'"}' \
| awk 'BEGIN {OFS="\t"} $9=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,"MALE",$10} $9=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"FEMALE",$10}' \
| awk 'BEGIN {OFS="\t"} $10=="-9" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} $10=="0" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} \
$10=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"UNAFFECTED"} $10=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"AFFECTED"}' \
>> $CORE_PATH/$PROJECT/TEMP/SAMPLE_META.txt

#### ADD GENDER CHECK ###
### THIS IS TEMPORARY FOR NOW MAYBE
### THIS IS THE HEADER
## SM_TAG,X_AVG_DP,X_NORM_DP,Y_AVG_DP,Y_NORM_DP
##########################3

awk 'BEGIN {OFS="\t"} $2=="X"&&$3=="whole" {print $6,$7} $2=="Y"&&$3=="whole" {print $6,$7}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ANEUPLOIDY_CHECK/$SM_TAG".chrom_count_report.txt" \
| paste - - \
| awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
>> $CORE_PATH/$PROJECT/TEMP/GENDER_CHECK.TXT

# # GRABBING CONCORDANCE. MULTI-SAMPLE # This will be just for the validation runs.
# 
# ls $CORE_PATH/$PROJECT/REPORTS/CONCORDANCE_MS/*_concordance.csv \
# | awk '{print "awk","1",$0}' \
# | bash \
# | sort -r \
# | uniq \
# | awk 'NR>1' \
# | sort \
# | sed 's/,/\t/g' \
# | awk 'BEGIN {print "SM_TAG","COUNT_DISC_HOM","COUNT_CONC_HOM","PERCENT_CONC_HOM",\
# "COUNT_DISC_HET","COUNT_CONC_HET","PERCENT_CONC_HET",\
# "PERCENT_TOTAL_CONC","COUNT_HET_BEADCHIP","SENSITIVITY_2_HET"} \
# {print $1,$5,$6,$7,$2,$3,$4,$8,$9,$10}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/CONCORDANCE_MS.txt
# 

# GRABBING VERIFY BAM ID #
## THIS IS THE HEADER ##
## awk 'BEGIN {print "SM_TAG""\t""VERIFYBAM_FREEMIX""\t""VERIFYBAM_#SNPS""\t""VERIFYBAM_FREELK1""\t""VERIFYBAM_FREELK0""\t""VERIFYBAM_AVG_DP"} \
#############

awk 'BEGIN {OFS="\t"} NR>1 {print $1,$7,$4,$8,$9,$6}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID/$SM_TAG".selfSM" \
>> $CORE_PATH/$PROJECT/TEMP/VERFIY_BAM_ID.TXT

#### GRABBING INSERT SIZE ####
## THIS IS THE HEADER ##
## BEGIN {print "SM_TAG","MEDIAN_INSERT_SIZE","MEDIAN_ABSOLUTE_DEVIATION_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE"} ##
#################################3

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$1,$2,$5,$6}' $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/METRICS/$SM_TAG".insert_size_metrics.txt" \
>> $CORE_PATH/$PROJECT/TEMP/INSERT_SIZE_METRICS.TXT

# GRABBING ALIGNMENT SUMMARY METRICS FOR READ 1 # good
## THIS THE HEADER ##
## {print "SM_TAG","PF_NOISE_READS_R1","PF_READS_ALIGNED_R1","PCT_PF_READS_ALIGNED_R1","PF_ALIGNED_BASES_R1","PF_HQ_ALIGNED_READS_R1",\
## "PF_HQ_ALIGNED_BASES_R1",\
## "PF_HQ_ALIGNED_Q20_BASES_R1","PF_HQ_MEDIAN_MISMATCHES_R1","PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1",\
## "PCT_READS_ALIGNED_IN_PAIRS_R1","BAD_CYCLES_R1",\
## "STRAND_BALANCE_R1","PCT_ADAPTER_R1"} \
##############################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$18,$19,$20,$22}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_1_METRICS.TXT

# # GRABBING ALIGNMENT SUMMARY METRICS FOR READ 2 # good
## THIS THE HEADER ##
## {print "SM_TAG","PF_NOISE_READS_R2","PF_READS_ALIGNED_R2","PCT_PF_READS_ALIGNED_R2","PF_ALIGNED_BASES_R2","PF_HQ_ALIGNED_READS_R2",\
## "PF_HQ_ALIGNED_BASES_R2",\
## "PF_HQ_ALIGNED_Q20_BASES_R2","PF_HQ_MEDIAN_MISMATCHES_R2","PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2",\
## "PCT_READS_ALIGNED_IN_PAIRS_R2","BAD_CYCLES_R2",\
## "STRAND_BALANCE_R1","PCT_ADAPTER_R2"} \
################################################

awk 'BEGIN {OFS="\t"} NR==9 {print "'$SM_TAG'",$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$18,$19,$20,$22}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_2_METRICS.TXT

# # GRABBING ALIGNMENT SUMMARY METRICS FOR PAIR # good
## THIS THE HEADER ##
## {print "SM_TAG","TOTAL_READS","RAW_GIGS","PF_NOISE_READS_PAIR","PF_READS_ALIGNED_PAIR","PCT_PF_READS_ALIGNED_PAIR",\
## "PF_ALIGNED_BASES_PAIR","PF_HQ_ALIGNED_READS_PAIR",\
## "PF_HQ_ALIGNED_BASES_PAIR",\
## "PF_HQ_ALIGNED_Q20_BASES_PAIR","PF_HQ_MEDIAN_MISMATCHES_PAIR","PF_MISMATCH_RATE_PAIR","PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR",\
## "MEAN_READ_LENGTH","PCT_READS_ALIGNED_IN_PAIRS_PAIR","BAD_CYCLES_PAIR",\
## "STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR","PCT_ADAPTER_PAIR"} \
################################################

awk 'BEGIN {OFS="\t"} NR==10 {print "'$SM_TAG'",$2,($2*$16/1000000000),$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$18,$19,$20,$21,$22}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_PAIR_METRICS.TXT

# # GRABBING MARK DUPLICATES REPORT # good
## THIS IS THE HEADER ##
## {print "SM_TAG","UNPAIRED_READS_EXAMINED","READ_PAIRS_EXAMINED","UNMAPPED_READS","UNPAIRED_READ_DUPLICATES",\
## "READ_PAIR_DUPLICATES","READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE"}
##########################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$2,$3,$4,$5,$6,$7,$8,$9}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" \
>> $CORE_PATH/$PROJECT/TEMP/MARK_DUPLICATES_METRICS.TXT

# # GRABBING HYB SELECTION REPORT
## THIS IS THE HEADER ##
##{print "SM_TAG","BAIT_SET","GENOME_SIZE","BAIT_TERRITORY","TARGET_TERRITORY","PF_UNIQUE_READS","PCT_PF_UQ_READS","PF_UQ_READS_ALIGNED",\
# "PCT_PF_UQ_READS_ALIGNED",\
# "PF_UQ_BASES_ALIGNED","PF_UQ_GIGS_ALIGNED","PCT_SELECTED_BASES","MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE",\
# "ZERO_CVG_TARGETS_PCT","PCT_EXC_DUPE","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET",\
# "PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X",\
# "PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","HS_LIBRARY_SIZE",\
# "AT_DROPOUT","GC_DROPOUT","HET_SNP_SENSITIVITY","HET_SNP_Q"}
##################33

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$1,$2,$3,$4,$8,$10,$11,$12,$14,($14/1000000000),$19,$22,$23,$24,$28,$29,$30,$31,$32,$33,\
$35,$36,$37,$38,$39,$40,$41,$42,$43,$50,$51,$52,$53}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
>> $CORE_PATH/$PROJECT/TEMP/HYB_SELECTION.TXT

## PULLING BAIT BIAS REPORT FOR Cref and Gref
### THIS IS THE HEADER ####
## SM_TAG,Cref_Q,Gref_Q

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/SUMMARY/$SM_TAG".bait_bias_summary_metrics.txt" \
| sed '/^$/d' \
| awk 'BEGIN {OFS="\t"} $12=="Cref"||$12=="Gref"  {print $5}' \
| paste - - \
| awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
>> $CORE_PATH/$PROJECT/TEMP/BAIT_BIAS.TXT


## PULLING PRE-ADAPTER BIAS REPORT FOR Deamination and OxoG
### THIS IS THE HEADER ####
## SM_TAG,Deamination_Q,OxoG_Q

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt" \
| sed '/^$/d' \
| awk 'BEGIN {OFS="\t"} $12=="Deamination"||$12=="OxoG"  {print $5}' \
| paste - - \
| awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
>> $CORE_PATH/$PROJECT/TEMP/PRE_ADAPTER.TXT

# # GRABBING THE QUALITY YIELD REPORT. CHANGE COUNTS TO PERCENTAGES.
## THIS IS THE HEADER ##
## {print "SM_TAG","PCT_Q20_BASES","PCT_Q30_BASES"}
##########################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$7/$5*100,$9/$5*100}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/QUALITY_YIELD/$SM_TAG".quality_yield_metrics.txt" \
>> $CORE_PATH/$PROJECT/TEMP/QUALITY_YIELD_METRICS.TXT

# GENERATE COUNT PCT,IN DBSNP FOR ON BAIT SNVS
## THIS IS THE HEADER ##
# {print "SM_TAG""\t""COUNT_SNV_ON_BAIT""\t""PERCENT_SNV_ON_BAIT_SNP138"} 
#############################################33

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_BAIT/$SM_TAG".SNV.ON_BAIT.PASS.vcf" \
| awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (SNV_COUNT>=l) {print "'$SM_TAG'",SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
else {print "'$SM_TAG'","0","NaN"}}' \
| sed 's/ /\t/g' \
>> $CORE_PATH/$PROJECT/TEMP/BAIT_SNV_METRICS.TXT

# GENERATE COUNT PCT,IN DBSNP FOR ON TARGET SNVS
## THIS IS THE HEADER ##
# {print "SM_TAG""\t""COUNT_SNV_ON_TARGET""\t""PERCENT_SNV_ON_TARGET_SNP138"} 
#############################################33

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_TARGET/$SM_TAG".SNV.ON_TARGET.PASS.vcf" \
| awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (SNV_COUNT>=l) {print "'$SM_TAG'",SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
else {print "'$SM_TAG'","0","NaN"}}' \
| sed 's/ /\t/g' \
>> $CORE_PATH/$PROJECT/TEMP/TARGET_SNV_METRICS.TXT

# GRABBING TI/TV ON UCSC CODING EXONS, ALL
## This is the Header ##
# {print "SM_TAG""\t""ALL_TI_TV_COUNT""\t""ALL_TI_TV_RATIO"}
#################################3

awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
else {print "'$SM_TAG'","0","NaN"}}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_All_.titv.txt" \
>> $CORE_PATH/$PROJECT/TEMP/TITV_ALL.TXT

# GRABBING TI/TV ON UCSC CODING EXONS, KNOWN
## This is the Header ##
# {print "SM_TAG""\t""KNOWN_TI_TV_COUNT""\t""KNOWN_TI_TV_RATIO"}
#################################3

awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
else {print "'$SM_TAG'","0","NaN"}}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Known_.titv.txt" \
>> $CORE_PATH/$PROJECT/TEMP/TITV_KNOWN.TXT

# GRABBING TI/TV ON UCSC CODING EXONS, NOVEL
## This is the Header ##
# {print "SM_TAG""\t""NOVEL_TI_TV_COUNT""\t""NOVEL_TI_TV_RATIO"}
#################################3

awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
else {print "'$SM_TAG'","0","NaN"}}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Novel_.titv.txt" \
>> $CORE_PATH/$PROJECT/TEMP/TITV_NOVEL.TXT

# # Generate a boatload of Indel metrics on bait
#

# INDEL METRICS ON BAIT
# THIS IS THE HEADER
# {print "SM_TAG","COUNT_ALL_INDEL_BAIT","ALL_INDEL_BAIT_PCT_SNP138","COUNT_BIALLELIC_INDEL_BAIT","BIALLELIC_INDEL_BAIT_PCT_SNP138",\
# "I_GT15_COUNT_BAIT","I_GT15_BAIT_PCT_SNP138","D_GT15_COUNT_BAIT","D_GT15_BAIT_PCT_SNP138",\
# "I_LT15_COUNT_BAIT","I_LT15_BAIT_PCT_SNP138","D_LT15_COUNT_BAIT","D_LT15_BAIT_PCT_SNP138",\
# "I_D_RATIO_BAIT","I_D_RATIO_GT15_BAIT","I_D_RATIO_LT15_BAIT","I_1AND2_VS_3_BAIT","D_1AND2_VS_3_BAIT","BIGGEST_INSERTION","BIGGEST_DELETION"}'
########################################################

(grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_BAIT/$SM_TAG".INDEL.ON_BAIT.PASS.vcf" \
| awk '{INDEL_COUNT++NR} \
{INDEL_BIALLELIC+=($5!~",")} \
{DBSNP_COUNT+=($3~"rs")} \
{DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
{BIG_INSERTION_BIALLELIC+=(length($5)-length($4))>15&&$5!~","} \
{BIG_INSERTION_BIALLELIC_DBSNP+=(length($5)-length($4))>15&&$5!~","&&$3~"rs"} \
{BIG_DELETION_BIALLELIC+=(length($5)-length($4))<-15&&$5!~","} \
{BIG_DELETION_BIALLELIC_DBSNP+=(length($5)-length($4))<-15&&$5!~","&&$3~"rs"} \
{SMALL_INSERTION_BIALLELIC+=((length($5)-length($4))<=15&&(length($5)-length($4))>0)&&$5!~","} \
{SMALL_INSERTION_BIALLELIC_DBSNP+=((length($5)-length($4))<=15&&(length($5)-length($4))>0)&&$5!~","&&$3~"rs"} \
{SMALL_DELETION_BIALLELIC+=((length($5)-length($4))>=-15&&(length($5)-length($4))<0)&&$5!~","} \
{SMALL_DELETION_BIALLELIC_DBSNP+=((length($5)-length($4))>=-15&&(length($5)-length($4))<0)&&$5!~","&&$3~"rs"} \
{ONE_AND_TWO_BP_INSERTION_BIALLELIC+=((length($5)-length($4))<=2&&(length($5)-length($4))>0)&&$5!~","} \
{THREE_BP_INSERTION_BIALLELIC+=(length($5)-length($4))==3&&$5!~","} \
{ONE_AND_TWO_BP_DELETION_BIALLELIC+=((length($5)-length($4))>=-2&&(length($5)-length($4))<0)&&$5!~","} \
{THREE_BP_DELETION_BIALLELIC+=(length($5)-length($4))==-3&&$5!~","} \
END {print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100,\
BIG_INSERTION_BIALLELIC,(BIG_INSERTION_BIALLELIC_DBSNP/BIG_INSERTION_BIALLELIC*100),BIG_DELETION_BIALLELIC,\
(BIG_DELETION_BIALLELIC_DBSNP/BIG_DELETION_BIALLELIC*100),\
SMALL_INSERTION_BIALLELIC,(SMALL_INSERTION_BIALLELIC_DBSNP/SMALL_INSERTION_BIALLELIC*100),SMALL_DELETION_BIALLELIC,\
(SMALL_DELETION_BIALLELIC_DBSNP/SMALL_DELETION_BIALLELIC*100),\
((BIG_INSERTION_BIALLELIC+SMALL_INSERTION_BIALLELIC)/(BIG_DELETION_BIALLELIC+SMALL_DELETION_BIALLELIC)),\
(BIG_INSERTION_BIALLELIC/BIG_DELETION_BIALLELIC),(SMALL_INSERTION_BIALLELIC/SMALL_DELETION_BIALLELIC),\
(ONE_AND_TWO_BP_INSERTION_BIALLELIC/THREE_BP_INSERTION_BIALLELIC),\
(ONE_AND_TWO_BP_DELETION_BIALLELIC/THREE_BP_DELETION_BIALLELIC)}' \
; grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_BAIT/$SM_TAG".INDEL.ON_BAIT.PASS.vcf" \
| awk '$5!~","' \
| awk '{INDEL_SIZE=(length($5)-length($4))} \
{if (MIN_INDEL_SIZE==0) {MIN_INDEL_SIZE=MAX_INDEL_SIZE=INDEL_SIZE}; \
if(INDEL_SIZE>MAX_INDEL_SIZE) {MAX_INDEL_SIZE=INDEL_SIZE}; \
if(INDEL_SIZE<MIN_INDEL_SIZE) {MIN_INDEL_SIZE=INDEL_SIZE}} \
END {print MAX_INDEL_SIZE,MIN_INDEL_SIZE}') \
| paste - - \
| sed 's/ /\t/g' \
>> $CORE_PATH/$PROJECT/TEMP/BAIT_INDEL_METRICS.TXT

# ##################################### EXAMPLE OF ABOVE
# 
# (grep -v "^#" NA00857-1Bl-1-3.INDEL.ON_TARGET.PASS.vcf \
# | awk '{INDEL_COUNT++NR} \
# {INDEL_BIALLELIC+=($5!~",")} \
# {DBSNP_COUNT+=($3~"rs")} \
# {DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
# {BIG_INSERTION_BIALLELIC+=(length($5)-length($4))>15&&$5!~","} \
# {BIG_INSERTION_BIALLELIC_DBSNP+=(length($5)-length($4))>15&&$5!~","&&$3~"rs"} \
# {BIG_DELETION_BIALLELIC+=(length($5)-length($4))<-15&&$5!~","} \
# {BIG_DELETION_BIALLELIC_DBSNP+=(length($5)-length($4))<-15&&$5!~","&&$3~"rs"} \
# {SMALL_INSERTION_BIALLELIC+=((length($5)-length($4))<=15&&(length($5)-length($4))>0)&&$5!~","} \
# {SMALL_INSERTION_BIALLELIC_DBSNP+=((length($5)-length($4))<=15&&(length($5)-length($4))>0)&&$5!~","&&$3~"rs"} \
# {SMALL_DELETION_BIALLELIC+=((length($5)-length($4))>=-15&&(length($5)-length($4))<0)&&$5!~","} \
# {SMALL_DELETION_BIALLELIC_DBSNP+=((length($5)-length($4))>=-15&&(length($5)-length($4))<0)&&$5!~","&&$3~"rs"} \
# {ONE_AND_TWO_BP_INSERTION_BIALLELIC+=((length($5)-length($4))<=2&&(length($5)-length($4))>0)&&$5!~","} \
# {THREE_BP_INSERTION_BIALLELIC+=(length($5)-length($4))==3&&$5!~","} \
# {ONE_AND_TWO_BP_DELETION_BIALLELIC+=((length($5)-length($4))>=-2&&(length($5)-length($4))<0)&&$5!~","} \
# {THREE_BP_DELETION_BIALLELIC+=(length($5)-length($4))==-3&&$5!~","} \
# END {print "SM_TAG",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100,\
# BIG_INSERTION_BIALLELIC,(BIG_INSERTION_BIALLELIC_DBSNP/BIG_INSERTION_BIALLELIC*100),BIG_DELETION_BIALLELIC,\
# (BIG_DELETION_BIALLELIC_DBSNP/BIG_DELETION_BIALLELIC*100),\
# SMALL_INSERTION_BIALLELIC,(SMALL_INSERTION_BIALLELIC_DBSNP/SMALL_INSERTION_BIALLELIC*100),SMALL_DELETION_BIALLELIC,\
# (SMALL_DELETION_BIALLELIC_DBSNP/SMALL_DELETION_BIALLELIC*100),\
# ((BIG_INSERTION_BIALLELIC+SMALL_INSERTION_BIALLELIC)/(BIG_DELETION_BIALLELIC+SMALL_DELETION_BIALLELIC)),\
# (BIG_INSERTION_BIALLELIC/BIG_DELETION_BIALLELIC),(SMALL_INSERTION_BIALLELIC/SMALL_DELETION_BIALLELIC),\
# (ONE_AND_TWO_BP_INSERTION_BIALLELIC/THREE_BP_INSERTION_BIALLELIC),\
# (ONE_AND_TWO_BP_DELETION_BIALLELIC/THREE_BP_DELETION_BIALLELIC)}' \
# ; grep -v "^#" NA00857-1Bl-1-3.INDEL.ON_TARGET.PASS.vcf \
# | awk '$5!~","' \
# | awk '{INDEL_SIZE=(length($5)-length($4))} \
# {if (MIN_INDEL_SIZE==0) {MIN_INDEL_SIZE=MAX_INDEL_SIZE=INDEL_SIZE}; \
# if(INDEL_SIZE>MAX_INDEL_SIZE) {MAX_INDEL_SIZE=INDEL_SIZE}; \
# if(INDEL_SIZE<MIN_INDEL_SIZE) {MIN_INDEL_SIZE=INDEL_SIZE}} \
# END {print MAX_INDEL_SIZE,MIN_INDEL_SIZE}') \
# | paste - - \
# | sed 's/ /\t/g' \
# >> $CORE_PATH/$PROJECT/TEMP/BAIT_INDEL_METRICS.TXT
# 
# #####################################################

# # Generate a boatload of Indel metrics on TARGET
#

# INDEL METRICS ON TARGET
# THIS IS THE HEADER
# {print "SM_TAG","COUNT_ALL_INDEL_TARGET","ALL_INDEL_TARGET_PCT_SNP138","COUNT_BIALLELIC_INDEL_TARGET","BIALLELIC_INDEL_TARGET_PCT_SNP138",\
# "I_GT15_COUNT_TARGET","I_GT15_TARGET_PCT_SNP138","D_GT15_COUNT_TARGET","D_GT15_TARGET_PCT_SNP138",\
# "I_LT15_COUNT_TARGET","I_LT15_TARGET_PCT_SNP138","D_LT15_COUNT_TARGET","D_LT15_TARGET_PCT_SNP138",\
# "I_D_RATIO_TARGET","I_D_RATIO_GT15_TARGET","I_D_RATIO_LT15_TARGET","I_1AND2_VS_3_TARGET","D_1AND2_VS_3_TARGET","BIGGEST_INSERTION","BIGGEST_DELETION"}'
########################################################

(grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_TARGET/$SM_TAG".INDEL.ON_TARGET.PASS.vcf" \
| awk '{INDEL_COUNT++NR} \
{INDEL_BIALLELIC+=($5!~",")} \
{DBSNP_COUNT+=($3~"rs")} \
{DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
{BIG_INSERTION_BIALLELIC+=(length($5)-length($4))>15&&$5!~","} \
{BIG_INSERTION_BIALLELIC_DBSNP+=(length($5)-length($4))>15&&$5!~","&&$3~"rs"} \
{BIG_DELETION_BIALLELIC+=(length($5)-length($4))<-15&&$5!~","} \
{BIG_DELETION_BIALLELIC_DBSNP+=(length($5)-length($4))<-15&&$5!~","&&$3~"rs"} \
{SMALL_INSERTION_BIALLELIC+=((length($5)-length($4))<=15&&(length($5)-length($4))>0)&&$5!~","} \
{SMALL_INSERTION_BIALLELIC_DBSNP+=((length($5)-length($4))<=15&&(length($5)-length($4))>0)&&$5!~","&&$3~"rs"} \
{SMALL_DELETION_BIALLELIC+=((length($5)-length($4))>=-15&&(length($5)-length($4))<0)&&$5!~","} \
{SMALL_DELETION_BIALLELIC_DBSNP+=((length($5)-length($4))>=-15&&(length($5)-length($4))<0)&&$5!~","&&$3~"rs"} \
{ONE_AND_TWO_BP_INSERTION_BIALLELIC+=((length($5)-length($4))<=2&&(length($5)-length($4))>0)&&$5!~","} \
{THREE_BP_INSERTION_BIALLELIC+=(length($5)-length($4))==3&&$5!~","} \
{ONE_AND_TWO_BP_DELETION_BIALLELIC+=((length($5)-length($4))>=-2&&(length($5)-length($4))<0)&&$5!~","} \
{THREE_BP_DELETION_BIALLELIC+=(length($5)-length($4))==-3&&$5!~","} \
END {print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100,\
BIG_INSERTION_BIALLELIC,(BIG_INSERTION_BIALLELIC_DBSNP/BIG_INSERTION_BIALLELIC*100),BIG_DELETION_BIALLELIC,\
(BIG_DELETION_BIALLELIC_DBSNP/BIG_DELETION_BIALLELIC*100),\
SMALL_INSERTION_BIALLELIC,(SMALL_INSERTION_BIALLELIC_DBSNP/SMALL_INSERTION_BIALLELIC*100),SMALL_DELETION_BIALLELIC,\
(SMALL_DELETION_BIALLELIC_DBSNP/SMALL_DELETION_BIALLELIC*100),\
((BIG_INSERTION_BIALLELIC+SMALL_INSERTION_BIALLELIC)/(BIG_DELETION_BIALLELIC+SMALL_DELETION_BIALLELIC)),\
(BIG_INSERTION_BIALLELIC/BIG_DELETION_BIALLELIC),(SMALL_INSERTION_BIALLELIC/SMALL_DELETION_BIALLELIC),\
(ONE_AND_TWO_BP_INSERTION_BIALLELIC/THREE_BP_INSERTION_BIALLELIC),\
(ONE_AND_TWO_BP_DELETION_BIALLELIC/THREE_BP_DELETION_BIALLELIC)}' \
; grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_TARGET/$SM_TAG".INDEL.ON_TARGET.PASS.vcf" \
| awk '$5!~","' \
| awk '{INDEL_SIZE=(length($5)-length($4))} \
{if (MIN_INDEL_SIZE==0) {MIN_INDEL_SIZE=MAX_INDEL_SIZE=INDEL_SIZE}; \
if(INDEL_SIZE>MAX_INDEL_SIZE) {MAX_INDEL_SIZE=INDEL_SIZE}; \
if(INDEL_SIZE<MIN_INDEL_SIZE) {MIN_INDEL_SIZE=INDEL_SIZE}} \
END {print MAX_INDEL_SIZE,MIN_INDEL_SIZE}') \
| paste - - \
| sed 's/ /\t/g' \
>> $CORE_PATH/$PROJECT/TEMP/TARGET_INDEL_METRICS.TXT

# BASIC METRICS FOR MIXED VARIANT TYPES ON BAIT

# GENERATE COUNT PCT,IN DBSNP FOR ON BAIT MIXED VARIANT
## THIS IS THE HEADER ##
# {print "SM_TAG""\t""COUNT_MIXED_ON_BAIT""\t""PERCENT_MIXED_ON_BAIT_SNP138"} 
#############################################33

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/FILTERED_ON_BAIT/$SM_TAG".MIXED.ON_BAIT.PASS.vcf" \
| awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (MIXED_COUNT>=l) {print "'$SM_TAG'",MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100} \
else {print "'$SM_TAG'","0","NaN"}}' \
| sed 's/ /\t/g' \
>> $CORE_PATH/$PROJECT/TEMP/BAIT_MIXED_METRICS.TXT

# GENERATE COUNT PCT,IN DBSNP FOR ON TARGET MIXED VARIANT
## THIS IS THE HEADER ##
# {print "SM_TAG""\t""COUNT_MIXED_ON_TARGET""\t""PERCENT_MIXED_ON_TARGET_SNP138"} 
#############################################33

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/FILTERED_ON_TARGET/$SM_TAG".MIXED.ON_TARGET.PASS.vcf" \
| awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (MIXED_COUNT>=l) {print "'$SM_TAG'",MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100} \
else {print "'$SM_TAG'","0","NaN"}}' \
| sed 's/ /\t/g' \
>> $CORE_PATH/$PROJECT/TEMP/TARGET_MIXED_METRICS.TXT

# 
# 
# # # GRABBING ANNOVAR METRICS # Hopefully don't have to do this.
# # 
# # ls $CORE_PATH/$PROJECT/REPORTS/ANNOVAR/*txt \
# # | awk '{split($1,SMtag,"/");print "awk \x27 BEGIN {FS=\x22\x5Ct\x22} \
# # NR>6 \
# # {total_snv+=($10~\x22Snv\x22)} \
# # {total_indel+=($10!~\x22Snv\x22)} \
# # {snv_126+=($10~\x22Snv\x22&&$55~\x22rs\x22)} \
# # {indel_126+=($10!~\x22Snv\x22&&$55~\x22rs\x22)} \
# # {snv_131+=($10~\x22Snv\x22&&$57~\x22rs\x22)} \
# # {indel_131+=($10!~\x22Snv\x22&&$57~\x22rs\x22)} \
# # END {print \x22"SMtag[9]"\x22,\
# # (snv_126/total_snv*100),\
# # (snv_131/total_snv*100),\
# # (indel_126/total_indel*100),\
# # (indel_131/total_indel*100)}\x27",\
# # "'$CORE_PATH'""/""'$PROJECT'""/REPORTS/ANNOVAR/"SMtag[9]}' \
# # | bash \
# # | sed 's/_MS_OnBait_ANNOVAR_REPORT.txt//g' \
# # | sed 's/ /\t/g' \
# # | awk 'BEGIN {print "SM_TAG""\t""PERCENT_SNV_ON_BAIT_SNP126""\t""PERCENT_SNV_ON_BAIT_SNP131""\t"\
# # "PERCENT_INDEL_ON_BAIT_SNP126""\t""PERCENT_INDEL_ON_BAIT_SNP131"} \
# # {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' \
# # >| $CORE_PATH/$PROJECT/TEMP/ANNOVAR_METRICS.TXT
# 

