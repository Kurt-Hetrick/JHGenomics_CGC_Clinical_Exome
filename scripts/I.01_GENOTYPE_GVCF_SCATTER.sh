# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q prod.q,rnd.q,test.q

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

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
REF_GENOME=$6
DBSNP=$7
CHROMOSOME=$8
CONTROL_GVCF_LIST=$9

# Genotype the gvcf files

START_GENOTYPE_GVCF=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--standard_min_confidence_threshold_for_calling 30 \
--standard_min_confidence_threshold_for_emitting 0 \
--annotation AS_BaseQualityRankSumTest \
--annotation AS_FisherStrand \
--annotation AS_InbreedingCoeff \
--annotation AS_MappingQualityRankSumTest \
--annotation AS_RMSMappingQuality \
--annotation AS_ReadPosRankSumTest \
--annotation AS_StrandOddsRatio \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--annotation LikelihoodRankSumTest \
-L $CHROMOSOME \
--variant $CONTROL_GVCF_LIST \
--variant $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".gvcf.list" \
-o $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW."$CHROMOSOME".vcf"

END_GENOTYPE_GVCF=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",I.01,GENOTYPE_GVCF_$CHROMOSOME,"$HOSTNAME","$START_GENOTYPE_GVCF","$END_GENOTYPE_GVCF \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--standard_min_confidence_threshold_for_calling 30 \
--standard_min_confidence_threshold_for_emitting 0 \
--annotation AS_BaseQualityRankSumTest \
--annotation AS_FisherStrand \
--annotation AS_InbreedingCoeff \
--annotation AS_MappingQualityRankSumTest \
--annotation AS_RMSMappingQuality \
--annotation AS_ReadPosRankSumTest \
--annotation AS_StrandOddsRatio \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--annotation LikelihoodRankSumTest \
-L $CHROMOSOME \
--variant $CONTROL_GVCF_LIST \
--variant $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".gvcf.list" \
-o $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW."$CHROMOSOME".vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY.COMMAND.LINES.txt
