#!/usr/bin/env python

# script to run test-run of SGA and PGA by varying population size and save 
# running time in a file

import sys
import subprocess
import os
from datetime import datetime

ROOT_DIR = os.path.abspath(".")
BIN_DIR_NAME = "bin"
DATA_DIR_NAME = "stats_data"
OUTFILE_PREF = "data_xPop"

BIN_DIR = ROOT_DIR + "/" + BIN_DIR_NAME
DATA_DIR = ROOT_DIR + "/" + DATA_DIR_NAME
SGA_PREF = "sga_P"
PGA_PREF = "pga_partial_P"

POPLIST = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]


def main():
    data_file = DATA_DIR + "/" + OUTFILE_PREF + "_6Nodes" + ".csv"
    
    with open(data_file, "a+") as f:
        
        print "Execution time comparison by varying population size."
        print "SGA vs. PGA (using 6 processors)"

        datestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        f.write("# " + datestamp + "\n")
        f.write("# Using: " + PGA_PREF + " and " + SGA_PREF + " binary\n")
        f.write("# POPSIZE,sga time,pga_time,speedup\n")

        os.chdir(BIN_DIR)
        
        for popsize in POPLIST:
            sga_file = BIN_DIR + "/" + SGA_PREF + str(popsize)
            pga_file = BIN_DIR + "/" + PGA_PREF + str(popsize)

            if os.path.isfile(sga_file) and os.path.isfile(pga_file):
                ps = subprocess.Popen(sga_file, stdout=subprocess.PIPE)
                ps.wait()
                val_str = ps.stdout.readline().rstrip()
                sga_time = float(val_str)

                ps = subprocess.Popen("mpiexec -f ~/hosts -n 6 " + pga_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                val_str = ps.stdout.readline().rstrip()
                pga_time = float(val_str)

                speedup =  sga_time / pga_time
                
                f.write(str(popsize) + "," + "{0:.4f}".format(sga_time) + "," + "{0:.4f}".format(pga_time) + "," + "{0:.4f}".format(speedup) + "\n")
                print "Popsize: %i SGA: %3.4f PGA: %3.4f  Speedup: %3.4f" %(popsize, sga_time, pga_time, speedup)
            else:
                print " Both version of executable file not found"
            
            print "Done!"
    # # Get filelist from bin dir
    # os.chdir(BIN_DIR)
    # ps = subprocess.Popen("ls", stdout=subprocess.PIPE)
    # binfiles = ps.stdout.readlines()
    #
    # os.chdir(DATA_DIR)
    # subprocess.call(["ls", "-l"])


if __name__ == "__main__":
    main();