#!/bin/sh

PROBLEMS_DIR=../problems/pure_int_bench
OUTPUT_SCRIPT=CPlexConsole_run.sh

echo "Problems path is: $PROBLEMS_DIR"
echo "Output script is: $OUTPUT_SCRIPT"

rm -f $OUTPUT_SCRIPT

echo "#!/bin/sh" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "cd ../build" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

for PROB_PATH in $PROBLEMS_DIR/*.gz
do
    FNAME=`basename $PROB_PATH`
    CUTS_TL=${1:-0.0} #cuts gen time limit
    CPLEX_TL=${2:-0.0} #cplex time limit
    MIN_DIM=${3:-0} #minimal problem dim
    MAX_DIM=${4:-0} #maximal problem dim
    CPLEX_WS=${5:-8} #workers for CPlex
    LOGF=${FNAME}_log.txt #out log file name
    LOGBF=${FNAME}_b_log.txt #out log file name for baseline
    CUTSF=${FNAME}_cuts.txt #out cuts file name
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo "#$FNAME" >> $OUTPUT_SCRIPT
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo ./CPlexConsole -i $PROB_PATH -o $LOGF --cuts_output $CUTSF --wcgn 1 --wcpl $CPLEX_WS \
    --cplex_tl $CPLEX_TL --cuts_tl $CUTS_TL --min_dim $MIN_DIM --max_dim $MAX_DIM >> $OUTPUT_SCRIPT
    echo ./CPlexConsole -i $PROB_PATH -o $LOGBF -b --wcpl $CPLEX_WS >> $OUTPUT_SCRIPT
    #echo >> $OUTPUT_SCRIPT
done
