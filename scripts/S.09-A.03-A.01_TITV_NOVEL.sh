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

SAMTOOLS_DIR=$1
CORE_PATH=$2

PROJECT=$3
FAMILY=$4
SM_TAG=$5

# Filter to just on all of the variants all

START_TITV_NOVEL=`date '+%s'`

$SAMTOOLS_DIR/bcftools/vcfutils.pl qstats \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".SNV.TITV_NOVEL.vcf" \
>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Novel_.titv.txt"

END_TITV_NOVEL=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",S.01,TITV_NOVEL,"$HOSTNAME","$START_TITV_NOVEL","$END_TITV_NOVEL \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $SAMTOOLS_DIR/bcftools/vcfutils qstats \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".SNV.TITV_NOVEL.vcf" \
\>\| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Novel_.titv.txt" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Novel_.titv.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
