#!/usr/bin/env python

# Python script to generate Run test programs, collect execution time, create data file for gnuplot,
# generate plot and save them.

import sys
import subprocess
import os

ROOT_DIR = os.path.abspath(".")
BIN_DIR_NAME = "bin"
DATA_DIR_NAME = "stats_data"
OUTFILE_PREF = "Exec_times_xNodes"

BIN_DIR = ROOT_DIR + "/" + BIN_DIR_NAME
DATA_DIR = ROOT_DIR + "/" + DATA_DIR_NAME

NODES = 6
POPSIZE = 5000


def main():
    data_file = DATA_DIR + "/" + OUTFILE_PREF + "_P"+ str(POPSIZE) + ".csv"
    
    with open(data_file, "a+") as f:
        f.write("Nodes | exec_time\n")
        
        for i in range(1, NODES):
            
            pgab_file = BIN_DIR + "/pga_bcast_P" + str(POPSIZE)

            if os.path.isfile(pgab_file):
                os.chdir(BIN_DIR)
                
                ps = subprocess.Popen("mpiexec -f ~/hosts -n "+ str(i) + " "+ pgab_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                pgab_time = ps.stdout.readline()
                
                print "PGA_B time:", pgab_time
                
                f.write(str(i) + " " + str(pgab_time).rstrip() + "\n")

            else:
                print " Both version of executable file not found"


    # # Get filelist from bin dir
    # os.chdir(BIN_DIR)
    # ps = subprocess.Popen("ls", stdout=subprocess.PIPE)
    # binfiles = ps.stdout.readlines()
    #
    # os.chdir(DATA_DIR)
    # subprocess.call(["ls", "-l"])


if __name__ == "__main__":
    main();