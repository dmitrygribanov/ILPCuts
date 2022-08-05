#!/bin/sh

PROBLEMS_PATH=../problems/pure_int_bench
OUTPUT_SCRIPT=CPlexConsole_run.sh

echo "Problems path is: $PROBLEMS_PATH"
echo "Output script is: $OUTPUT_SCRIPT"

rm -f $OUTPUT_SCRIPT

echo "#!/bin/sh" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

echo "cd ../build" >> $OUTPUT_SCRIPT
echo >> $OUTPUT_SCRIPT

for prob_path in $PROBLEMS_PATH/*.gz
do
    fname=`basename $prob_path`
    tl=$1 #time limit
    ws=$2 #workers for CPlex
    logf=${fname}_log.txt #out log file name
    logbf=${fname}_b_log.txt #out log file name for baseline
    cutsf=${fname}_cuts.txt #out cuts file name
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo "#$fname" >> $OUTPUT_SCRIPT
    echo "#-------------------------------------------------------" >> $OUTPUT_SCRIPT
    echo timeout $tl ./CPlexConsole -i $prob_path -o $logf --cuts_output $cutsf --wcgn 1 --wcpl $ws >> $OUTPUT_SCRIPT
    echo timeout $tl ./CPlexConsole -i $prob_path -o $logbf -b --wcpl $ws >> $OUTPUT_SCRIPT
    #echo >> $OUTPUT_SCRIPT
done
