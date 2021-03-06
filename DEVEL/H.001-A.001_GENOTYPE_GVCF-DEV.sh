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

JAVA_1_7=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
REF_GENOME=$6
DBSNP=$7

RIS_ID=${SM_TAG%@*}
BARCODE_2D=${SM_TAG#*@}

START_GENOTYPE_GVCF=`date '+%s'`

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--standard_min_confidence_threshold_for_calling 30 \
--standard_min_confidence_threshold_for_emitting 0 \
--variant $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".gvcf.list" \
-o $CORE_PATH/$PROJECT/$FAMILY/VCF/RAW/$FAMILY".RAW.vcf"

END_GENOTYPE_GVCF=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",I.001,GENOTYPE_GVCF,"$HOSTNAME","$START_GENOTYPE_GVCF","$END_GENOTYPE_GVCF \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--standard_min_confidence_threshold_for_calling 30 \
--standard_min_confidence_threshold_for_emitting 0 \
--variant $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".gvcf.list" \
-o $CORE_PATH/$PROJECT/$FAMILY/VCF/RAW/$FAMILY".RAW.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY.COMMAND.LINES.txt

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/RAW/$FAMILY".RAW.vcf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
