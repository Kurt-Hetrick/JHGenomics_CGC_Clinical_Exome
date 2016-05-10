#!/bin/bash

module load sge

SAMPLE_SHEET=$1
PED_FILE=$2

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

SCRIPT_DIR="/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts"

JAVA_1_7="/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin"
JAVA_1_8="/isilon/sequencing/Kurt/Programs/Java/jdk1.8.0_73/bin"
CORE_PATH="/isilon/sequencing/Seq_Proj/"
BWA_DIR="/isilon/sequencing/Kurt/Programs/BWA/bwa-0.7.8"
PICARD_DIR="/isilon/sequencing/Kurt/Programs/Picard/picard-tools-2.1.1"
# PICARD_DIR="/isilon/sequencing/VITO/Programs/picard/picard-tools-1.141"
GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.5-0"
# GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.3-0"
VERIFY_DIR="/isilon/sequencing/Kurt/Programs/VerifyBamID/verifyBamID_20120620/bin/"
GENE_LIST="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/RefSeqGene.GRCh37.Ready.txt"
VERIFY_VCF="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
CODING_BED="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/UCSC_hg19_CodingOnly_083013_MERGED_noContigs_noCHR.bed"
SAMTOOLS_DIR="/isilon/sequencing/Kurt/Programs/samtools/samtools-0.1.18"
TABIX_DIR="/isilon/sequencing/Kurt/Programs/TABIX/tabix-0.2.6"
CORE_PATH="/isilon/sequencing/Seq_Proj"
DATAMASH_DIR="/isilon/sequencing/Kurt/Programs/PATH"
CYTOBAND_BED="/isilon/sequencing/Kurt/CGC/GRCh37.Cytobands.bed"
# BEDTOOLS IS v2.22.0
BEDTOOLS_DIR="/isilon/sequencing/Kurt/Programs/PATH"

##### MAKE A DIRECTORY TREE ##### SHOULD BE COMPLETE #####

mkdir -p ~/CGC_PIPELINE_TEMP

MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
PED_PREFIX=`basename $PED_FILE .ped`

##########################################################

SETUP_PROJECT ()
{
FORMAT_MANIFEST
MERGE_PED_MANIFEST
CREATE_SAMPLE_INFO_ARRAY
MAKE_PROJ_DIR_TREE
}

FORMAT_MANIFEST ()
{
sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| sed 's/,/\t/g' \
| sort -k 8 \
>| ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt
}

MERGE_PED_MANIFEST ()
{
awk 1 $PED_FILE \
| sed 's/\r//g' \
| sort -k 2 \
| join -1 8 -2 2 ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt /dev/stdin \
| awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$1,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' \
>| ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt
}

# MAKE AN ARRAY FOR EACH SAMPLE
	## SAMPLE_INFO_ARRAY[0] = PROJECT
	## SAMPLE_INFO_ARRAY[1] = FAMILY
	## SAMPLE_INFO_ARRAY[2] = SM_TAG
		## SAMPLE = SM_TAG

CREATE_SAMPLE_INFO_ARRAY ()
{
SAMPLE_INFO_ARRAY=(`awk '$8=="'$SAMPLE'" {print $1,$19,$8}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

# PROJECT DIRECTORY TREE CREATOR

MAKE_PROJ_DIR_TREE ()
{
mkdir -p $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/LOGS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/BAM \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/HC_BAM \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/INDEL/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/SNV/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/MIXED/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/VCF/FILTERED_ON_BAIT/TABIX \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/GVCF \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,QUALITY_YIELD} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/CONCORDANCE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/COUNT_COVARIATES/{GATK_REPORT,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/DEPTH_OF_COVERAGE/{TARGET,UCSC_CODING_PLUS_10bp} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/INSERT_SIZE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/ANEUPLOIDY_CHECK \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/LOGS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/VCF/{RAW,VQSR} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/{TEMP,FASTQ,REPORTS,LOGS}
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
SETUP_PROJECT
done

############################################################

# to create the qsub cmd line to submit bwa alignments to the cluster
# handle blank lines
# handle something else too

awk '{split($18,INDEL,";");split($8,smtag,"[@-]"); \
print "qsub","-N","A.001_BWA_"$8"_"$2"_"$3"_"$4,\
"-o","'$CORE_PATH'/"$1"/"$19"/"$8"/LOGS/"$8"_"$2"_"$3"_"$4".BWA.log",\
"'$SCRIPT_DIR'""/A.001_BWA.sh",\
"'$BWA_DIR'","'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$19,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12"\n""sleep 3s"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt

# create a hold job id qsub command line based on the number of
# submit merging the bam files created by bwa mem above
# only launch when every lane for a sample is done being processed by bwa mem

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1,2,3 collapse 4 collapse 5 \
| awk 'BEGIN {FS="\t"} \
gsub(/,/,",A.001_BWA_"$3"_",$4) \
gsub(/,/,",INPUT=/isilon/sequencing/Seq_Proj/"$1"/TEMP/",$5) \
{print "qsub","-N","B.001_MERGE_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".MERGE.BAM.FILES.log",\
"-hold_jid","A.001_BWA_"$3"_"$4, \
"'$SCRIPT_DIR'""/B.001_MERGE_SORT_AGGRO.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,$3,"INPUT=/isilon/sequencing/Seq_Proj/"$1"/TEMP/"$5"\n""sleep 3s"}'

# Mark duplicates on the bam file above. Create a Mark Duplicates report which goes into the QC report

awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","C.001_MARK_DUPLICATES_"$3"_"$1,\
"-hold_jid","B.001_MERGE_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".MARK_DUPLICATES.log",\
"'$SCRIPT_DIR'""/C.001_MARK_DUPLICATES.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

# Generate a list of places that could be potentially realigned.

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@-]"); \
print "qsub","-N","D.001_REALIGNER_TARGET_CREATOR_"$3"_"$1,\
"-hold_jid","C.001_MARK_DUPLICATES_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".REALIGNER_TARGET_CREATOR.log",\
"'$SCRIPT_DIR'""/D.001_REALIGNER_TARGET_CREATOR.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2]"\n""sleep 3s"}'

# With the list generated above walk through the BAM file and realign where necessary
# Write out a new bam file

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@-]"); \
print "qsub","-N","E.001_INDEL_REALIGNER_"$3"_"$1,\
"-hold_jid","D.001_REALIGNER_TARGET_CREATOR_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".INDEL_REALIGNER.log",\
"'$SCRIPT_DIR'""/E.001_INDEL_REALIGNER.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2]"\n""sleep 3s"}'

# Run Base Quality Score Recalibration

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$18,$17,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@-]"); \
print "qsub","-N","F.001_PERFORM_BQSR_"$3"_"$1,\
"-hold_jid","E.001_INDEL_REALIGNER_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PERFORM_BQSR.log",\
"'$SCRIPT_DIR'""/F.001_PERFORM_BQSR.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2],$6,$7"\n""sleep 3s"}'

# write Final Bam file

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","G.001_FINAL_BAM_"$3"_"$1,\
"-hold_jid","F.001_PERFORM_BQSR_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".FINAL_BAM.log",\
"'$SCRIPT_DIR'""/G.001_FINAL_BAM.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

##### ALL H.00X SERIES OF SCRIPTS CAN BE RUN IN PARALLEL SINCE THEY ARE DEPENDENT ON FINAL BAM FILE GENERATION #####

# Run Haplotype Caller in GVCF mode

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.001_HAPLOTYPE_CALLER_"$1"_"$3,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".HAPLOTYPE_CALLER.log",\
"'$SCRIPT_DIR'""/H.001_HAPLOTYPE_CALLER.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

# Run POST BQSR TABLE

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$18,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@-]"); \
print "qsub","-N","H.002_POST_BQSR_TABLE_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".POST_BQSR_TABLE.log",\
"'$SCRIPT_DIR'""/H.002_POST_BQSR_TABLE.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2],$6"\n""sleep 3s"}'

# Run ANALYZE COVARIATES

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@-]"); \
print "qsub","-N","H.002-A.001_ANALYZE_COVARIATES_"$3"_"$1,\
"-hold_jid","H.002_POST_BQSR_TABLE_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".ANALYZE_COVARIATES.log",\
"'$SCRIPT_DIR'""/H.002-A.001_ANALYZE_COVARIATES.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

# RUN DOC CODING PLUS 10 BP FLANKS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.003_DOC_CODING_10bpFLANKS_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_CODING_10bpFLANKS.log",\
"'$SCRIPT_DIR'""/H.003_DOC_CODING_10bpFLANKS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'","'$CODING_BED'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 3s"}'

# RUN DOC TARGET BED

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.005_DOC_TARGET_BED_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_TARGET_BED.log",\
"'$SCRIPT_DIR'""/H.005_DOC_TARGET_BED.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

# RUN ANEUPLOIDY_CHECK AFTER DOC TARGET BED FINISHES

awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.005-A.001_DOC_CHROM_DEPTH_"$3"_"$1,\
"-hold_jid","H.005_DOC_TARGET_BED_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".ANEUPLOIDY_CHECK.log",\
"'$SCRIPT_DIR'""/H.005-A.001_CHROM_DEPTH.sh",\
"'$CORE_PATH'","'$CYTOBAND_BED'","'$DATAMASH_DIR'","'$BEDTOOLS_DIR'",$1,$2,$3"\n""sleep 3s"}'

# RUN COLLECT MULTIPLE METRICS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$17,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.006_COLLECT_MULTIPLE_METRICS_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".COLLECT_MULTIPLE_METRICS.log",\
"'$SCRIPT_DIR'""/H.006_COLLECT_MULTIPLE_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5,$6"\n""sleep 3s"}'

# RUN COLLECT HS METRICS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$15,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.007_COLLECT_HS_METRICS_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".COLLECT_HS_METRICS.log",\
"'$SCRIPT_DIR'""/H.007_COLLECT_HS_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5,$6"\n""sleep 3s"}'

# RUN SELECT VERIFYBAM ID VCF

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.008_SELECT_VERIFYBAMID_VCF_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".SELECT_VERIFYBAMID_VCF.log",\
"'$SCRIPT_DIR'""/H.008_SELECT_VERIFYBAMID_VCF.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'","'$VERIFY_VCF'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

# RUN VERIFYBAMID

awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.008-A.001_VERIFYBAMID_"$3"_"$1,\
"-hold_jid","H.008_SELECT_VERIFYBAMID_VCF_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".VERIFYBAMID.log",\
"'$SCRIPT_DIR'""/H.008-A.001_VERIFYBAMID.sh",\
"'$CORE_PATH'","'$VERIFY_DIR'",$1,$2,$3"\n""sleep 3s"}'

#### JOINT CALLING AND VQSR ####

### CREATE A GVCF ".list" file for each sample

CREATE_GVCF_LIST ()
{
awk 'BEGIN {OFS="/"} $19=="'$FAMILY'" {print "'$CORE_PATH'",$1,$19,$8,"GVCF",$8".g.vcf"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort \
| uniq \
>| $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/$FAMILY/$FAMILY".gvcf.list"
}

CREATE_FAMILY_INFO_ARRAY ()
{
FAMILY_INFO_ARRAY=(`awk '$19=="'$FAMILY'" {print $1,$19}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CREATE_FAMILY_SAMPLE_LIST ()
{
awk '$19=="'$FAMILY'" {print $8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort \
| uniq \
>| $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/$FAMILY/$FAMILY".sample.list"
}

for FAMILY in $( awk '{print $19}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq ) ;
do
CREATE_FAMILY_INFO_ARRAY
CREATE_GVCF_LIST
CREATE_FAMILY_SAMPLE_LIST
done

### Run GenotypeGVCF per Family

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1,2 collapse 3 unique 4 unique 5 \
| awk 'BEGIN {FS="\t"}
gsub (/,/,",H.001_HAPLOTYPE_CALLER_"$1"_",$3) \
$3~"," \
{print "qsub","-N","H.001-A.001_GENOTYPE_GVCF_"$2"_"$1,\
"-hold_jid","H.001_HAPLOTYPE_CALLER_"$1"_"$3,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".GENOTYPE_GVCF.log",\
"'$SCRIPT_DIR'""/H.001-A.001_GENOTYPE_GVCF.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$4,$5"\n""sleep 3s"} \
$3!~"," \
{print "qsub","-N","H.001-A.001_GENOTYPE_GVCF_"$2"_"$1,\
"-hold_jid","H.001_HAPLOTYPE_CALLER_"$1"_"$3,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".GENOTYPE_GVCF.log",\
"'$SCRIPT_DIR'""/H.001-A.001_GENOTYPE_GVCF.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$4,$5"\n""sleep 3s"}'

### Run Variant Recalibrator for the SNP model, this is done in parallel with the INDEL model

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","H.001-A.001-A.001_VARIANT_RECALIBRATOR_SNP_"$2"_"$1,\
"-hold_jid","H.001-A.001_GENOTYPE_GVCF_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_RECALIBRATOR_SNP.log",\
"'$SCRIPT_DIR'""/H.001-A.001-A.001_VARIANT_RECALIBRATOR_SNP.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### Run Variant Recalibrator for the INDEL model, this is done in parallel with the SNP model

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","H.001-A.001-A.002_VARIANT_RECALIBRATOR_INDEL_"$2"_"$1,\
"-hold_jid","H.001-A.001_GENOTYPE_GVCF_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_RECALIBRATOR_INDEL.log",\
"'$SCRIPT_DIR'""/H.001-A.001-A.002_VARIANT_RECALIBRATOR_INDEL.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### Run Apply Recalbration with the SNP model to the VCF file

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","H.001-A.001-A.001-A.001_APPLY_RECALIBRATION_SNP_"$2"_"$1,\
"-hold_jid","H.001-A.001-A.001_VARIANT_RECALIBRATOR_SNP_"$2"_"$1",""H.001-A.001-A.002_VARIANT_RECALIBRATOR_INDEL_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".APPLY_RECALIBRATION_SNP.log",\
"'$SCRIPT_DIR'""/H.001-A.001-A.001-A.001_APPLY_RECALIBRATION_SNP.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### Run Apply Recalibration with the INDEL model to the VCF file.

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","H.001-A.001-A.001-A.001-A.001_APPLY_RECALIBRATION_INDEL_"$2"_"$1,\
"-hold_jid","H.001-A.001-A.001-A.001_APPLY_RECALIBRATION_SNP_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".APPLY_RECALIBRATION_INDEL.log",\
"'$SCRIPT_DIR'""/H.001-A.001-A.001-A.001-A.001_APPLY_RECALIBRATION_INDEL.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### Add all possible GATK annotations to the VCF file.

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-hold_jid","H.001-A.001-A.001-A.001-A.001_APPLY_RECALIBRATION_INDEL_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".ADD_MORE_ANNOTATION.log",\
"'$SCRIPT_DIR'""/M.01_ADD_MORE_ANNOTATION.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'","'$PED_FILE'",$1,$2,$3"\n""sleep 3s"}'

##### DOING VCF BREAKOUTS #####

### SUBSETTING FROM COHORT (FAMILY PLUS CONTROL SET) VCF ###

# FILTER TO JUST VARIANT SITES

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","M.01-A.01_FILTER_COHORT_VARIANT_ONLY_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_VARIANT_ONLY.log",\
"'$SCRIPT_DIR'""/M.01-A.01_FILTER_COHORT_VARIANT_ONLY.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

# FILTER TO JUST PASSING VARIANT SITES

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","M.01-A.02_FILTER_COHORT_VARIANT_ONLY_PASS_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_VARIANT_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.02_FILTER_COHORT_VARIANT_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### SUBSETTING TO FAMILY ###

## SUBSETTING TO FAMILY ALL SITES ###

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","M.01-A.03_FILTER_TO_FAMILY_ALL_SITES_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_TO_FAMILY_ALL_SITES.log",\
"'$SCRIPT_DIR'""/M.01-A.03_FILTER_TO_FAMILY_ALL_SITES.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSETTING TO FAMILY ALL VARIANT SITES ##

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","M.01-A.04_FILTER_TO_FAMILY_VARIANTS_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_TO_FAMILY_VARIANTS.log",\
"'$SCRIPT_DIR'""/M.01-A.04_FILTER_TO_FAMILY_VARIANTS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSETTING TO FAMILY PASSING VARIANTS ##

awk 'BEGIN {OFS="\t"} {print $1,$19,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","M.01-A.05_FILTER_TO_FAMILY_VARIANTS_PASS_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_TO_FAMILY_VARIANTS_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.05_FILTER_TO_FAMILY_VARIANTS_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### SUBSETTING TO SAMPLE VCFS ###

## SUBSET TO SAMPLE VCF ALL SITES ##

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.06_FILTER_TO_SAMPLE_ALL_SITES_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_ALL_SITES.log",\
"'$SCRIPT_DIR'""/M.01-A.06_FILTER_TO_SAMPLE_ALL_SITES.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE VARIANTS ONLY 

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.07_FILTER_TO_SAMPLE_VARIANTS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_VARIANTS.log",\
"'$SCRIPT_DIR'""/M.01-A.07_FILTER_TO_SAMPLE_VARIANTS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING VARIANTS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.08_FILTER_TO_SAMPLE_VARIANTS_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_VARIANTS_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.08_FILTER_TO_SAMPLE_VARIANTS_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING SNVS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09_FILTER_TO_SNV_ONLY_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_SNV_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.09_FILTER_TO_SAMPLE_SNV_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING INDELS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.10_FILTER_TO_INDEL_ONLY_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_INDEL_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.10_FILTER_TO_SAMPLE_INDEL_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING MIXED

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.11_FILTER_TO_MIXED_ONLY_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_MIXED_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.11_FILTER_TO_SAMPLE_MIXED_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO TARGET SNV ONLY PASS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_SNV_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

## SUBSET TO TARGET INDEL ONLY PASS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_INDEL_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

## SUBSET TO TARGET MIXED ONLY PASS

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_MIXED_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/M.01-A.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

## SUBSET TO SAMPLE VCF ALL SITES ON TARGET##

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_ALL_SITES_TARGET.log",\
"'$SCRIPT_DIR'""/M.01-A.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

## SUBSET TO SAMPLE VARIANTS ONLY ON TARGET

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.16_FILTER_TO_SAMPLE_VARIANTS_TARGET_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_VARIANTS_TARGET.log",\
"'$SCRIPT_DIR'""/M.01-A.16_FILTER_TO_SAMPLE_VARIANTS_TARGET.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING VARIANTS ON TARGET

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.17_FILTER_TO_SAMPLE_VARIANTS_PASS_TARGET_"$3"_"$2"_"$1,\
"-hold_jid","M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_VARIANTS_PASS_TARGET.log",\
"'$SCRIPT_DIR'""/M.01-A.17_FILTER_TO_SAMPLE_VARIANTS_PASS_TARGET.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

### TITV SECTION ###

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$14}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"$3"_"$2"_"$1,\
"-hold_jid","M.01-A.09_FILTER_TO_SNV_ONLY_PASS_"$3"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF.log",\
"'$SCRIPT_DIR'""/M.01-A.09-A.01_FILTER_TO_SAMPLE_TITV_VCF.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND OVERLAP WITH DBSNP 129

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$14}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"$3"_"$2"_"$1,\
"-hold_jid","M.01-A.09_FILTER_TO_SNV_ONLY_PASS_"$3"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF_KNOWN.log",\
"'$SCRIPT_DIR'""/M.01-A.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND DO NOT OVERLAP WITH DBSNP 129

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$14}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"$3"_"$2"_"$1,\
"-hold_jid","M.01-A.09_FILTER_TO_SNV_ONLY_PASS_"$3"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF_NOVEL.log",\
"'$SCRIPT_DIR'""/M.01-A.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

### RUN TITV FOR THE PASSING SNVS THAT FALL IN UCSC CODING REGIONS THAT TOUCH EITHER THE BED OR TARGET FILE

## ALL SNVS TITV

awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09-A.01-A.01_TITV_ALL_"$3"_"$2"_"$1,\
"-hold_jid","M.01-A.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"$3"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_ALL.log",\
"'$SCRIPT_DIR'""/M.01-A.09-A.01-A.01_TITV_ALL.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## ALL KNOWN SNVS TITV

awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09-A.02-A.01_TITV_KNOWN_"$3"_"$2"_"$1,\
"-hold_jid","M.01-A.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"$3"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_KNOWN.log",\
"'$SCRIPT_DIR'""/M.01-A.09-A.02-A.01_TITV_KNOWN.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## ALL NOVEL SNVS TITV

awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{print "qsub","-N","M.01-A.09-A.03-A.01_TITV_NOVEL_"$3"_"$2"_"$1,\
"-hold_jid","M.01-A.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"$3"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_NOVEL.log",\
"'$SCRIPT_DIR'""/M.01-A.09-A.03-A.01_TITV_NOVEL.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

### QC REPORT PREP ###

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$20,$21,$22,$23}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk 'BEGIN {FS="\t"}
{print "qsub","-N","X.01-QC_REPORT_PREP_"$1"_"$3,\
"-hold_jid","M.01-A.03_FILTER_TO_FAMILY_ALL_SITES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".QC_REPORT_PREP.log",\
"'$SCRIPT_DIR'""/X.01-QC_REPORT_PREP.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'","'$DATAMASH_DIR'",$1,$2,$3,$4,$5,$6,$7"\n""sleep 3s"}'

### END PROJECT TASKS ###

# awk 'BEGIN {OFS="\t"} {print $1,$19,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | $DATAMASH_DIR/datamash -s -g 1,2 collapse 3 \
# | awk 'BEGIN {FS="\t"}
# gsub (/,/,",X.01-QC_REPORT_PREP_"$1"_",$3) \
# {print "qsub","-N","X.01-X.01-END_PROJECT_TASKS_"$1,\
# "-hold_jid","X.01-QC_REPORT_PREP_"$1"_"$3",""M.01_ADD_MORE_ANNOTATION_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$1".END_PROJECT_TASKS.log",\
# "'$SCRIPT_DIR'""/X.01-X.01-END_PROJECT_TASKS.sh",\
# "'$CORE_PATH'","'$DATAMASH_DIR'",$1"\n""sleep 3s"}'

awk 'BEGIN {OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1 collapse 2 \
| awk 'BEGIN {FS="\t"}
gsub (/,/,",X.01-QC_REPORT_PREP_"$1"_",$2) \
{print "qsub","-N","X.01-X.01-END_PROJECT_TASKS_"$1,\
"-hold_jid","X.01-QC_REPORT_PREP_"$1"_"$2,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".END_PROJECT_TASKS.log",\
"'$SCRIPT_DIR'""/X.01-X.01-END_PROJECT_TASKS.sh",\
"'$CORE_PATH'","'$DATAMASH_DIR'",$1"\n""sleep 3s"}'

### kEY FOR BLAH ###
#
#      1  CGC_160212_HJLWVBCXX_CGCDev1_TEST
#      2  HJLWVBCXX
#      3  1
#      4  ATGCCTAA
#      5  ILLUMINA
#      6  A01_NA12878
#      7  2/12/2016
#      8  NA12878
#      9  CGC
#     10  HiSeq2500_RapidRun
#     11  HJLWVBCXX_1_ATGCCTAA_A01_NA12878
#     12  /isilon/sequencing/GATK_resource_bundle/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.fasta
#     13  MBS
#     14  /isilon/sequencing/data/Work/BED/Production_BED_files/TsTv_BED_File_Agilent_ClinicalExome_S06588914_OnExon_merged_021015_noCHR.bed
#     15  /isilon/sequencing/data/Work/BED/Production_BED_files/ALLBED_BED_File_Agilent_ClinicalExome_S06588914_ALLBed_merged_021015_noCHR.bed
#     16  /isilon/sequencing/data/Work/BED/Production_BED_files/Targets_BED_File_Agilent_ClinicalExome_S06588914_OnTarget_merged_noCHR_013015.bed
#     17  /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.vcf
#     18  /isilon/sequencing/GATK_resource_bundle/2.2/b37/1000G_phase1.indels.b37.vcf;/isilon/sequencing/GATK_resource_bundle/2.2/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
#     19  XC01463
#     20  NA12891
#     21  NA12892
#     22  2
#     23  2
#######


###### SAMPLE MANIFEST KEY...NOT SURE WHAT I AM GOING TO END UP DOING HERE ######

# PROJECT=$1 # the Seq Proj folder name. 1st column in sample manifest
# FLOWCELL=$2 # flowcell that sample read group was performed on. 2nd column of sample manifest
# LANE=$3 # lane of flowcell that sample read group was performed on. 3rd column of the sample manifest
# INDEX=$4 # sample barcode. 4th column of the sample manifest
# PLATFORM=$5 # type of sequencing chemistry matching SAM specification. 5th column of the sample manifest.
# LIBRARY_NAME=$6 # library group of the sample read group.
# 								# Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not
# 								# 6th column of the sample manifest
# RUN_DATE=$7 # should be the run set up date to match the seq run folder name, but it has been arbitrarily populated. field X of manifest.
# SM_TAG=$8 # sample ID. sample name for all files, etc. field X of manifest
# CENTER=$9 # the center/funding mechanism. field X of manifest.
# DESCRIPTION=${10} # Generally we use to denote the sequencer setting (e.g. rapid run). field X of manifest.
# REF_GENOME=${11} # the reference genome used in the analysis pipeline. field X of manifest.
# TI_TV_BED=${12} # populated from sample manifest. where ucsc coding exons overlap with bait and target bed files
# BAIT_BED=${13} # populated from sample manifest. a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
# 								# Used for limited where to run base quality score recalibration on where to create gvcf files.
# TARGET_BED=${14} # populated from sample manifest. bed file acquired from manufacturer of their targets. field X of sample manifest.
# DBSNP=${15} # populated from sample manifest. used to annotate ID field in VCF file. masking in base call quality score recalibration.
# KNOWN_INDEL_1=${16} # populated from sample manifest. used for BQSR masking, sensitivity in local realignment.
# KNOWN_INDEL_2=${17} # populated from sample manifest. used for BQSR masking, sensitivity in local realignment.
#
# RIS_ID=${SM_TAG%@*} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
# BARCODE_2D=${SM_TAG#*@} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
#
####################################################################################

#### BOILERPLATE...I HAVE NOT DECIDED WHAT I AM GOING TO DO HERE######

# function GRAB_MANIFEST {
# sed 's/\r//g' $SAMPLE_SHEET \
# | awk 'BEGIN {FS=","} NR>1 \
# {split($19,INDEL,";");split($8,smtag,"@");print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2]}'
# }
#
# function GRAB_PROJECT_NAMES {
# PROJECT_NAMES=`sed 's/\r//g' $SAMPLE_SHEET \
# | awk 'BEGIN {FS=","} NR>1 print $1}'`
# }
######################################################################
