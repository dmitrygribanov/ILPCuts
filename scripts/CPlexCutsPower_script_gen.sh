#!/bin/sh

PROBLEMS_DIR=../problems/pure_int_bench
OUTPUT_SCRIPT=CPlexCutsPower_run.sh

echo "Problems path is: $PROBLEMS_DIR"
echo "Output script is: $OUTPUT_SCRIPT"
echo

rm -f $OUTPUT_SCRIPT

echo "#!/bin/sh" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "cp -r -f ../primes_table ../build/primes_table" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "OUT_DIR=CPlexCutsPower_\`date +\"%d.%m.%Y_%a_%H-%M\"\`" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "rm -r -f ../build/\$OUT_DIR" >> $OUTPUT_SCRIPT
echo "mkdir ../build/\$OUT_DIR" >> $OUTPUT_SCRIPT
echo "cd ../build/\$OUT_DIR" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

for PROB_PATH in $PROBLEMS_DIR/*.gz
do 
    FNAME=`basename $PROB_PATH`
    echo "Processing: $FNAME ..."

    CPLEX_WS=${1:-8} #workers for CPlex
    MIN_DIM=${2:-0} #minimal problem dim
    MAX_DIM=${3:-0} #maximal problem dim
    CUTS_TL=${4:-0.0} #cuts gen time limit in minutes
    CPLEX_TL=${5:-0.0} #cplex time limit in minutes
    MODEL_CUTS_SWITCH=${6:-on} #switch for model cuts option
    NODE_LIM=${7:-1} #node limit for CPlex
    CUTS_TL=`echo "$CUTS_TL * 60" | bc` #transforming into TL in seconds
    CPLEX_TL=`echo "$CPLEX_TL * 60" | bc` #transforming into TL in seconds
    LOGF=${FNAME}_log.txt #out log file name
    RESF=${FNAME}_results.txt #results file name
    CUTSF=${FNAME}_cuts.txt #out cuts file name
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo "#$FNAME" >> $OUTPUT_SCRIPT
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo ../CPlexCutsPower -i ../$PROB_PATH -o $LOGF --cuts_output $CUTSF --comp_result $RESF \
    --wcgn 1 --wcpl $CPLEX_WS --node_lim $NODE_LIM --model_cuts $MODEL_CUTS_SWITCH \
    --cplex_tl $CPLEX_TL --cuts_tl $CUTS_TL --min_dim $MIN_DIM --max_dim $MAX_DIM >> $OUTPUT_SCRIPT
    #echo >> $OUTPUT_SCRIPT
done

chmod 755 $OUTPUT_SCRIPT