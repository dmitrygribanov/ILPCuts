import subprocess
import os

os.system("./build/CPlexConsole -i ./problems/pure_int_bench/sorrell3.mps.gz -o build/sorrell3_log.txt --cuts_output build/sorrell3_cuts.txt --wcgn 1 --wcpl 15")

#subprocess.run(["./build/CPlexConsole", "-i", "./problems/pure_int_bench/sorrell3.mps.gz", "-o", "sorrell3_log.txt", 
#"--cuts_output", "sorrell3_cuts.txt", "--wcgn", "1", "--wcpl", "15"])