#!/bin/sh

PROBLEMS_DIR=../problems/pure_int_bench
OUTPUT_SCRIPT=run_rootNodeSolve.sh

echo "Problems path is: $PROBLEMS_DIR"
echo "Output script is: $OUTPUT_SCRIPT"
echo

rm -f $OUTPUT_SCRIPT

echo "#!/bin/sh" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "cp -r -f ../primes_table ../build/primes_table" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "OUT_DIR=rootNodeSolve_\`date +\"%d.%m.%Y_%a_%H-%M\"\`" >> $OUTPUT_SCRIPT
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
    SUBMAT_SIZE_LIM=${6:-0.0} #maximal size of square submat to generate conic cover cuts
    ROW_TRIALS_NUM=${7:-1} #number of row trials to generate cuts
    CUTS_GEN_ROUNDS=${8:-1} #number of cuts generation iterations
    CUTS_TL=`echo "$CUTS_TL * 60" | bc` #transforming into TL in seconds
    CPLEX_TL=`echo "$CPLEX_TL * 60" | bc` #transforming into TL in seconds
    LOGF=${FNAME}_log.txt #out log file name
    LOGBF=${FNAME}_b_log.txt #out log file name for baseline
    CUTSF=${FNAME}_cuts.txt #out cuts file name
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo "#$FNAME" >> $OUTPUT_SCRIPT
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo ../CPlexConsole -i ../$PROB_PATH -o $LOGF --cuts_output $CUTSF --wcgn 1 --wcpl $CPLEX_WS --iters $CUTS_GEN_ROUNDS --trials_num $ROW_TRIALS_NUM \
    --cplex_tl $CPLEX_TL --cuts_tl $CUTS_TL --min_dim $MIN_DIM --max_dim $MAX_DIM --submat_limit $SUBMAT_SIZE_LIM >> $OUTPUT_SCRIPT
    echo >> $OUTPUT_SCRIPT
    echo ../CPlexConsole -i ../$PROB_PATH -o $LOGBF -b --wcpl $CPLEX_WS \
    --cplex_tl $CPLEX_TL --min_dim $MIN_DIM --max_dim $MAX_DIM >> $OUTPUT_SCRIPT
    #echo >> $OUTPUT_SCRIPT
done

chmod 755 $OUTPUT_SCRIPT
