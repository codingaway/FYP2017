#!/usr/bin/env python

# Python script to generate Run test programs, collect execution time, create data file for gnuplot,
# generate plot and save them.

import sys
import subprocess
import os
from datetime import datetime

ROOT_DIR = os.path.abspath(".")
BIN_DIR_NAME = "bin"
DATA_DIR_NAME = "stats_data"
OUTFILE_PREF = "Amdahls_xNodes"

BIN_DIR = ROOT_DIR + "/" + BIN_DIR_NAME
DATA_DIR = ROOT_DIR + "/" + DATA_DIR_NAME
SGA_PREF = "sga_P"
PGA_PREF = "pga_global_P"

NODES = 6
POPSIZE = 5000


def main():
    data_file = DATA_DIR + "/" + OUTFILE_PREF + "_P"+ str(POPSIZE) + ".csv"
    
    with open(data_file, "a+") as f:
        sga_file = BIN_DIR + "/" + SGA_PREF + str(POPSIZE)
        pga_file = BIN_DIR + "/" + PGA_PREF + str(POPSIZE)
        
        if os.path.isfile(sga_file) and os.path.isfile(pga_file):
            datestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
            f.write("# " + datestamp + " Population size: " + str(POPSIZE) + "\n")
            f.write("# Using: " + PGA_PREF + str(POPSIZE)
                + " and " + SGA_PREF + str(POPSIZE) + " binary\n")
            f.write("# Nodes,sga time,pga_time,speedup\n")
            
            os.chdir(BIN_DIR) 
            ps = subprocess.Popen(sga_file, stdout=subprocess.PIPE)
            ps.wait()
            val_str = ps.stdout.readline().rstrip()
            sga_time = float(val_str)

            for i in range(1, NODES + 1):
                ps = subprocess.Popen("mpiexec -f ~/hosts -n "+ str(i) + " "+ pga_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                val_str = ps.stdout.readline().rstrip()
                pga_time = float(val_str)

                speedup =  sga_time / pga_time

                print "PGA time:", pga_time

                f.write(str(i) + "," + "{0:.4f}".format(sga_time) + "," + "{0:.4f}".format(pga_time) + "," + "{0:.4f}".format(speedup) + "\n")
        else:
            print "Both version of executable file not found"


    # # Get filelist from bin dir
    # os.chdir(BIN_DIR)
    # ps = subprocess.Popen("ls", stdout=subprocess.PIPE)
    # binfiles = ps.stdout.readlines()
    #
    # os.chdir(DATA_DIR)
    # subprocess.call(["ls", "-l"])


if __name__ == "__main__":
    main();