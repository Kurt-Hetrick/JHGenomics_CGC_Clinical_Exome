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

TABIX_DIR=$1
CORE_PATH=$2

PROJECT=$3
FAMILY=$4

# Filter to just on all of the variants all

START_VARIANT_TO_TABLE_BGZIP_FAMILY=`date '+%s'`

# bgzip the per family vcf converted into a table file

$TABIX_DIR/bgzip \
-c /home/sandbox/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.txt" \
>| $CORE_PATH/$PROJECT/$FAMILY/VCF/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.txt.gz"

# delete the uncompressed version of the table since it is being written to the sandbox.
# not doing this on the cgc side since the file i/o is not as bad as on cidr's production side and thus everything is being read/written on isilon

rm -rvf /home/sandbox/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.txt"

END_VARIANT_TO_TABLE_BGZIP_FAMILY=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",U.01,VARIANT_TO_TABLE_BGZIP_"$FAMILY"_ALL_SITES,"$HOSTNAME","$START_VARIANT_TO_TABLE_BGZIP_FAMILY","$END_VARIANT_TO_TABLE_BGZIP_FAMILY \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $TABIX_DIR/bgzip \
-c /home/sandbox/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.txt" \
\>\| $CORE_PATH/$PROJECT/$FAMILY/VCF/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.txt.gz" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.txt.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
