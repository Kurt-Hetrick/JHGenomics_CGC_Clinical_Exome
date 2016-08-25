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

START_VARIANT_TO_TABLE_COHORT=`date '+%s'`

# hopefully this is does everything in one iteration.
# it's still heavily bound by I/O.
# This might be worth it to write to sys temp and then in the next step (bgzip) write back out to isilon.
## Except that sys temp is very small on the c6220s...very small...sigh...I'm actually surprised that things are not crashing...

( cat $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.1.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.2.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.3.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.4.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.5.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.6.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.7.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.8.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.9.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.10.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.11.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.12.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.13.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.14.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.15.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.16.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.17.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.18.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.19.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.20.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.21.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.22.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.X.txt ; \
awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.Y.txt ) \
>| /home/sandbox/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.txt

END_VARIANT_TO_TABLE_COHORT=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",T.01,GATHER_VARIANT_TO_TABLE_COHORT_ALL_SITES,"$HOSTNAME","$START_VARIANT_TO_TABLE_COHORT","$END_VARIANT_TO_TABLE_COHORT \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# ( cat $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.1.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.2.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.3.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.4.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.5.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.6.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.7.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.8.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.9.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.10.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.11.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.12.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.13.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.14.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.15.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.16.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.17.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.18.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.19.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.20.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.21.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.22.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.X.txt ; \
# awk 'NR>1' $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.Y.txt ) \
# >| $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY.VQSR.ANNOTATED.txt

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
