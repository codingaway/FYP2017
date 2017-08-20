#!/usr/bin/env python

# Python script to generate Run test programs, collect execution time, create data file for gnuplot,
# generate plot and save them.

import sys
import subprocess
import os
from datetime import datetime

ROOT_DIR = os.path.abspath(".")
BIN_DIR_NAME = "bin"
#DATA_DIR_NAME = "stats_data"
#OUTFILE_PREF = "data_xNodes"

BIN_DIR = ROOT_DIR + "/" + BIN_DIR_NAME
#DATA_DIR = ROOT_DIR + "/" + DATA_DIR_NAME
SGA_PREF = "sga_Report_P"
PGA_PREF = "pga_partial_Report_P"

NODES = 6
POPSIZE = 5000


def main():
    
    sga_file = BIN_DIR + "/" + SGA_PREF + str(POPSIZE)
    pga_file = BIN_DIR + "/" + PGA_PREF + str(POPSIZE)

    print "GA report analysis"
    print "SGA vs. PGA (using population size: %i)." %POPSIZE

    if os.path.isfile(sga_file) and os.path.isfile(pga_file):

        os.chdir(BIN_DIR) 
        ps = subprocess.Popen(sga_file)
        ps.wait()

        ps = subprocess.Popen("mpiexec -f ~/hosts -n "+ str(NODES) + " "+ pga_file, shell=True)
        ps.wait()
        print "Done!"
    else:
        print "Both version of executables are not found"


if __name__ == "__main__":
    main();