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
PED_FILE=$4

PROJECT=$5
FAMILY=$6
REF_GENOME=$7


START_ADD_MORE_ANNOTATION=`date '+%s'`

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--annotation AlleleBalance \
--annotation AlleleBalanceBySample \
--annotation AlleleCountBySample \
--annotation BaseCounts \
--annotation BaseCountsBySample \
--annotation GCContent \
--annotation GenotypeSummaries \
--annotation HardyWeinberg \
--annotation HomopolymerRun \
--annotation InbreedingCoeff \
--annotation LikelihoodRankSumTest \
--annotation LowMQ \
--annotation MVLikelihoodRatio \
--annotation NBaseCount \
--annotation SampleList \
--annotation TandemRepeatAnnotator \
--annotation TransmissionDisequilibriumTest \
--annotation VariantType \
--pedigree $PED_FILE \
--pedigreeValidationType SILENT \
--variant $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.vcf" \
-L $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.vcf" \
-o $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.vcf"

END_ADD_MORE_ANNOTATION=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",M.001,ADD_MORE_ANNOTATION,"$HOSTNAME","$START_ADD_MORE_ANNOTATION","$END_ADD_MORE_ANNOTATION \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--useAllAnnotations \
--annotation AlleleBalance \
--annotation AlleleBalanceBySample \
--annotation AlleleCountBySample \
--annotation BaseCounts \
--annotation BaseCountsBySample \
--annotation GCContent \
--annotation GenotypeSummaries \
--annotation HardyWeinberg \
--annotation HomopolymerRun \
--annotation InbreedingCoeff \
--annotation LikelihoodRankSumTest \
--annotation LowMQ \
--annotation MVLikelihoodRatio \
--annotation NBaseCount \
--annotation SampleList \
--annotation TandemRepeatAnnotator \
--annotation TransmissionDisequilibriumTest \
--annotation VariantType \
--pedigree $PED_FILE \
--pedigreeValidationType SILENT \
--variant $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.vcf" \
-L $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.vcf" \
-o $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.vcf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.vcf.idx" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
