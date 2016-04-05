#!/bin/bash

SAMPLE_SHEET=$1
PED_FILE=$2

JAVA_1_7="/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin"
JAVA_1_8="/isilon/sequencing/Kurt/Programs/Java/jdk1.8.0_73/bin"
CORE_PATH="/isilon/sequencing/Seq_Proj/"
BWA_DIR="/isilon/sequencing/Kurt/Programs/BWA/bwa-0.7.8/"
PICARD_DIR="/isilon/sequencing/Kurt/Programs/Picard/picard-tools-2.1.1"
# PICARD_DIR="/isilon/sequencing/VITO/Programs/picard/picard-tools-1.141"
GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.5-0"
# GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.3-0"
VERIFY_DIR="/isilon/sequencing/Kurt/Programs/VerifyBamID/verifyBamID_20120620/bin/"
GENE_LIST="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/RefSeqGene.GRCh37.Ready.txt"
VERIFY_VCF="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
CODING_BED="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/UCSC_hg19_CodingOnly_083013_MERGED_noContigs_noCHR.bed"
SAMTOOLS_DIR="/isilon/sequencing/Kurt/Programs/samtools/samtools-0.1.18/"
TABIX_DIR="/isilon/sequencing/Kurt/Programs/TABIX/tabix-0.2.6/"
CORE_PATH="/isilon/sequencing/Seq_Proj"
DATAMASH_DIR="/isilon/sequencing/Kurt/Programs/PATH/"
CYTOBAND_BED="/isilon/sequencing/Kurt/CGC/GRCh37.Cytobands.bed"
# BEDTOOLS IS v2.22.0
BEDTOOLS_DIR="/isilon/sequencing/Kurt/Programs/PATH/"

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
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/{TEMP,FASTQ,REPORTS} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/LOGS
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/A.001_BWA.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/B.001_MERGE_SORT_AGGRO.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/C.001_MARK_DUPLICATES.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/D.001_REALIGNER_TARGET_CREATOR.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/E.001_INDEL_REALIGNER.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/F.001_PERFORM_BQSR.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/G.001_FINAL_BAM.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.001_HAPLOTYPE_CALLER.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

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

for FAMILY in $( awk '{print $19}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq ) ;
do
CREATE_FAMILY_INFO_ARRAY
CREATE_GVCF_LIST
done

### Run GenotypeGVCF per Family

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1,2 collapse 3 unique 4 \
| awk 'BEGIN {FS="\t"}
gsub (/,/,",H.001_HAPLOTYPE_CALLER_"$1"_",$3) \
{print "qsub","-N","H.001-A.0001_GENOTYPE_GVCF_"$2"_"$1,\
"-hold_jid","H.001_HAPLOTYPE_CALLER_"$1"_"$3,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".GENOTYPE_GVCF.log",\
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.001-A.001_GENOTYPE_GVCF.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$4"\n""sleep 3s"}'

# Run POST BQSR TABLE

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$18,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@-]"); \
print "qsub","-N","H.002_POST_BQSR_TABLE_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".POST_BQSR_TABLE.log",\
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.002_POST_BQSR_TABLE.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.002-A.001_ANALYZE_COVARIATES.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.003_DOC_CODING_10bpFLANKS.sh",\
"'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'","'$CODING_BED'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 3s"}'

# # RUN DOC SUPERSET BED ### TAKING THIS OUT
# 
# awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@-]"); \
# print "qsub","-N","H.004_DOC_SUPERSET_BED_"$3"_"$1,\
# "-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_SUPERSET_BED.log",\
# "/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.004_DOC_SUPERSET_BED.sh",\
# "'$JAVA_1_7'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4,$5"\n""sleep 3s"}'

# RUN DOC TARGET BED

awk 'BEGIN {OFS="\t"} {print $1,$19,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@-]"); \
print "qsub","-N","H.005_DOC_TARGET_BED_"$3"_"$1,\
"-hold_jid","G.001_FINAL_BAM_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_TARGET_BED.log",\
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.005_DOC_TARGET_BED.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.005-A.001_CHROM_DEPTH.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.006_COLLECT_MULTIPLE_METRICS.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.007_COLLECT_HS_METRICS.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.008_SELECT_VERIFYBAMID_VCF.sh",\
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
"/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome/scripts/H.008-A.001_VERIFYBAMID.sh",\
"'$CORE_PATH'","'$VERIFY_DIR'",$1,$2,$3"\n""sleep 3s"}'


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
