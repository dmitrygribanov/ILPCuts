#!/bin/sh

cd ../build

#-------------------------------------------------------
#air05.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/air05.mps.gz -o air05.mps.gz_log.txt --cuts_output air05.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/air05.mps.gz -o air05.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#bley_xs1.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/bley_xs1.mps.gz -o bley_xs1.mps.gz_log.txt --cuts_output bley_xs1.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/bley_xs1.mps.gz -o bley_xs1.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#bley_xs1noM.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/bley_xs1noM.mps.gz -o bley_xs1noM.mps.gz_log.txt --cuts_output bley_xs1noM.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/bley_xs1noM.mps.gz -o bley_xs1noM.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#bnatt400.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/bnatt400.mps.gz -o bnatt400.mps.gz_log.txt --cuts_output bnatt400.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/bnatt400.mps.gz -o bnatt400.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#bnatt500.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/bnatt500.mps.gz -o bnatt500.mps.gz_log.txt --cuts_output bnatt500.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/bnatt500.mps.gz -o bnatt500.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#cdc7-4-3-2.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/cdc7-4-3-2.mps.gz -o cdc7-4-3-2.mps.gz_log.txt --cuts_output cdc7-4-3-2.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/cdc7-4-3-2.mps.gz -o cdc7-4-3-2.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#circ10-3.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/circ10-3.mps.gz -o circ10-3.mps.gz_log.txt --cuts_output circ10-3.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/circ10-3.mps.gz -o circ10-3.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#cod105_vc1024.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/cod105_vc1024.mps.gz -o cod105_vc1024.mps.gz_log.txt --cuts_output cod105_vc1024.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/cod105_vc1024.mps.gz -o cod105_vc1024.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#cvs16r128-89.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/cvs16r128-89.mps.gz -o cvs16r128-89.mps.gz_log.txt --cuts_output cvs16r128-89.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/cvs16r128-89.mps.gz -o cvs16r128-89.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#d20200.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/d20200.mps.gz -o d20200.mps.gz_log.txt --cuts_output d20200.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/d20200.mps.gz -o d20200.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#eil33-2.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/eil33-2.mps.gz -o eil33-2.mps.gz_log.txt --cuts_output eil33-2.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/eil33-2.mps.gz -o eil33-2.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#enlight_hard.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/enlight_hard.mps.gz -o enlight_hard.mps.gz_log.txt --cuts_output enlight_hard.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/enlight_hard.mps.gz -o enlight_hard.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#f2000.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/f2000.mps.gz -o f2000.mps.gz_log.txt --cuts_output f2000.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/f2000.mps.gz -o f2000.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#glass-sc.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/glass-sc.mps.gz -o glass-sc.mps.gz_log.txt --cuts_output glass-sc.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/glass-sc.mps.gz -o glass-sc.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#graph20-20-1rand.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/graph20-20-1rand.mps.gz -o graph20-20-1rand.mps.gz_log.txt --cuts_output graph20-20-1rand.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/graph20-20-1rand.mps.gz -o graph20-20-1rand.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#hypothyroid-k1.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/hypothyroid-k1.mps.gz -o hypothyroid-k1.mps.gz_log.txt --cuts_output hypothyroid-k1.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/hypothyroid-k1.mps.gz -o hypothyroid-k1.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#k1mushroom.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/k1mushroom.mps.gz -o k1mushroom.mps.gz_log.txt --cuts_output k1mushroom.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/k1mushroom.mps.gz -o k1mushroom.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#mzzv11.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/mzzv11.mps.gz -o mzzv11.mps.gz_log.txt --cuts_output mzzv11.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/mzzv11.mps.gz -o mzzv11.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3024952-loue.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3024952-loue.mps.gz -o neos-3024952-loue.mps.gz_log.txt --cuts_output neos-3024952-loue.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3024952-loue.mps.gz -o neos-3024952-loue.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3083819-nubu.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3083819-nubu.mps.gz -o neos-3083819-nubu.mps.gz_log.txt --cuts_output neos-3083819-nubu.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3083819-nubu.mps.gz -o neos-3083819-nubu.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3355120-tarago.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3355120-tarago.mps.gz -o neos-3355120-tarago.mps.gz_log.txt --cuts_output neos-3355120-tarago.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3355120-tarago.mps.gz -o neos-3355120-tarago.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3355323-arnon.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3355323-arnon.mps.gz -o neos-3355323-arnon.mps.gz_log.txt --cuts_output neos-3355323-arnon.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3355323-arnon.mps.gz -o neos-3355323-arnon.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3381206-awhea.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3381206-awhea.mps.gz -o neos-3381206-awhea.mps.gz_log.txt --cuts_output neos-3381206-awhea.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3381206-awhea.mps.gz -o neos-3381206-awhea.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3426085-ticino.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3426085-ticino.mps.gz -o neos-3426085-ticino.mps.gz_log.txt --cuts_output neos-3426085-ticino.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3426085-ticino.mps.gz -o neos-3426085-ticino.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3426132-dieze.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3426132-dieze.mps.gz -o neos-3426132-dieze.mps.gz_log.txt --cuts_output neos-3426132-dieze.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3426132-dieze.mps.gz -o neos-3426132-dieze.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3603137-hoteo.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3603137-hoteo.mps.gz -o neos-3603137-hoteo.mps.gz_log.txt --cuts_output neos-3603137-hoteo.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3603137-hoteo.mps.gz -o neos-3603137-hoteo.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-3682128-sandon.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-3682128-sandon.mps.gz -o neos-3682128-sandon.mps.gz_log.txt --cuts_output neos-3682128-sandon.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-3682128-sandon.mps.gz -o neos-3682128-sandon.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-4360552-sangro.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-4360552-sangro.mps.gz -o neos-4360552-sangro.mps.gz_log.txt --cuts_output neos-4360552-sangro.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-4360552-sangro.mps.gz -o neos-4360552-sangro.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-4738912-atrato.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-4738912-atrato.mps.gz -o neos-4738912-atrato.mps.gz_log.txt --cuts_output neos-4738912-atrato.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-4738912-atrato.mps.gz -o neos-4738912-atrato.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos-950242.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos-950242.mps.gz -o neos-950242.mps.gz_log.txt --cuts_output neos-950242.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos-950242.mps.gz -o neos-950242.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#neos859080.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/neos859080.mps.gz -o neos859080.mps.gz_log.txt --cuts_output neos859080.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/neos859080.mps.gz -o neos859080.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#nu25-pr12.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/nu25-pr12.mps.gz -o nu25-pr12.mps.gz_log.txt --cuts_output nu25-pr12.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/nu25-pr12.mps.gz -o nu25-pr12.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#opm2-z10-s4.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/opm2-z10-s4.mps.gz -o opm2-z10-s4.mps.gz_log.txt --cuts_output opm2-z10-s4.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/opm2-z10-s4.mps.gz -o opm2-z10-s4.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#pb-market-split8-70-4.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/pb-market-split8-70-4.mps.gz -o pb-market-split8-70-4.mps.gz_log.txt --cuts_output pb-market-split8-70-4.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/pb-market-split8-70-4.mps.gz -o pb-market-split8-70-4.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#peg-solitaire-a3.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/peg-solitaire-a3.mps.gz -o peg-solitaire-a3.mps.gz_log.txt --cuts_output peg-solitaire-a3.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/peg-solitaire-a3.mps.gz -o peg-solitaire-a3.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#pythago7825.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/pythago7825.mps.gz -o pythago7825.mps.gz_log.txt --cuts_output pythago7825.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/pythago7825.mps.gz -o pythago7825.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#qap10.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/qap10.mps.gz -o qap10.mps.gz_log.txt --cuts_output qap10.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/qap10.mps.gz -o qap10.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#ramos3.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/ramos3.mps.gz -o ramos3.mps.gz_log.txt --cuts_output ramos3.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/ramos3.mps.gz -o ramos3.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#reblock115.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/reblock115.mps.gz -o reblock115.mps.gz_log.txt --cuts_output reblock115.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/reblock115.mps.gz -o reblock115.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#rococoB10-011000.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/rococoB10-011000.mps.gz -o rococoB10-011000.mps.gz_log.txt --cuts_output rococoB10-011000.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/rococoB10-011000.mps.gz -o rococoB10-011000.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#rococoC10-001000.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/rococoC10-001000.mps.gz -o rococoC10-001000.mps.gz_log.txt --cuts_output rococoC10-001000.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/rococoC10-001000.mps.gz -o rococoC10-001000.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#seymour.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/seymour.mps.gz -o seymour.mps.gz_log.txt --cuts_output seymour.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/seymour.mps.gz -o seymour.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#sorrell3.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/sorrell3.mps.gz -o sorrell3.mps.gz_log.txt --cuts_output sorrell3.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/sorrell3.mps.gz -o sorrell3.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#sorrell4.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/sorrell4.mps.gz -o sorrell4.mps.gz_log.txt --cuts_output sorrell4.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/sorrell4.mps.gz -o sorrell4.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#sorrell7.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/sorrell7.mps.gz -o sorrell7.mps.gz_log.txt --cuts_output sorrell7.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/sorrell7.mps.gz -o sorrell7.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#splice1k1.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/splice1k1.mps.gz -o splice1k1.mps.gz_log.txt --cuts_output splice1k1.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/splice1k1.mps.gz -o splice1k1.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#supportcase2.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/supportcase2.mps.gz -o supportcase2.mps.gz_log.txt --cuts_output supportcase2.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/supportcase2.mps.gz -o supportcase2.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#supportcase22.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/supportcase22.mps.gz -o supportcase22.mps.gz_log.txt --cuts_output supportcase22.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/supportcase22.mps.gz -o supportcase22.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#supportcase30.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/supportcase30.mps.gz -o supportcase30.mps.gz_log.txt --cuts_output supportcase30.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/supportcase30.mps.gz -o supportcase30.mps.gz_b_log.txt -b --wcpl 15
#-------------------------------------------------------
#wachplan.mps.gz
#-------------------------------------------------------
./CPlexConsole -i ../problems/pure_int_bench/wachplan.mps.gz -o wachplan.mps.gz_log.txt --cuts_output wachplan.mps.gz_cuts.txt --wcgn 1 --wcpl 15 --cplex_tl 1800 --cuts_tl 1800 --min_dim 0 --max_dim 4000
./CPlexConsole -i ../problems/pure_int_bench/wachplan.mps.gz -o wachplan.mps.gz_b_log.txt -b --wcpl 15
