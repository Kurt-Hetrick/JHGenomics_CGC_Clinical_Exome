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

JAVA_1_7=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6
REF_GENOME=$7

START_VARIANT_TO_TABLE_SAMPLE=`date '+%s'`

# not doing --splitMultiallelic here...maybe do one as an example and discuss with Molly
# do an example of molten output to look at/show molly

( cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.1.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.2.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.3.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.4.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.5.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.6.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.7.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.8.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.9.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.10.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.11.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.12.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.13.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.14.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.15.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.16.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.17.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.18.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.19.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.20.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.21.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.22.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.X.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.Y.txt \
| grep "^CHROM" ; \
cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.1.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.2.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.3.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.4.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.5.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.6.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.7.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.8.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.9.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.10.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.11.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.12.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.13.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.14.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.15.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.16.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.17.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.18.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.19.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.20.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.21.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.22.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.X.txt \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.Y.txt \
| grep -v "^CHROM" ) \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.txt"

END_VARIANT_TO_TABLE_SAMPLE=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",T.01,GATHER_VARIANT_TO_TABLE_"$SAMPLE"_ALL_SITES,"$HOSTNAME","$START_VARIANT_TO_TABLE_SAMPLE","$END_VARIANT_TO_TABLE_SAMPLE \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# echo ( cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.1.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.2.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.3.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.4.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.5.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.6.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.7.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.8.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.9.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.10.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.11.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.12.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.13.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.14.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.15.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.16.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.17.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.18.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.19.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.20.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.21.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.22.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.X.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.Y.txt \
# \| grep "^CHROM" ; \
# cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.1.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.2.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.3.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.4.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.5.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.6.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.7.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.8.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.9.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.10.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.11.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.12.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.13.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.14.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.15.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.16.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.17.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.18.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.19.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.20.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.21.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.22.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.X.txt \
# $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.Y.txt \
# \| grep -v "^CHROM" ) \
# \>\| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.txt" \
# >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
