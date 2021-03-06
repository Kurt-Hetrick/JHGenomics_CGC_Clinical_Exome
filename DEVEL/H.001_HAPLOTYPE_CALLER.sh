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
SM_TAG=$6
REF_GENOME=$7
BAIT_BED=$8

RIS_ID=${SM_TAG%@*}
BARCODE_2D=${SM_TAG#*@}

## -----Haplotype Caller-----

## Call on Bait (padded or superset)

START_HAPLOTYPE_CALLER=`date '+%s'`

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $REF_GENOME \
--input_file $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
-L $BAIT_BED \
--emitRefConfidence BP_RESOLUTION \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-A DepthPerSampleHC \
-A ClippingRankSumTest \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A QualByDepth \
-pairHMM VECTOR_LOGLESS_CACHING \
-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf"

END_HAPLOTYPE_CALLER=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.001,HAPLOTYPE_CALLER,"$HOSTNAME","$START_HAPLOTYPE_CALLER","$END_HAPLOTYPE_CALLER \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $REF_GENOME \
--input_file $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
-L $BAIT_BED \
--emitRefConfidence BP_RESOLUTION \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-A DepthPerSampleHC \
-A ClippingRankSumTest \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A QualByDepth \
-pairHMM VECTOR_LOGLESS_CACHING \
-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf.idx" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
