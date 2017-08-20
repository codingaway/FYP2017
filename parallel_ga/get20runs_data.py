#!/usr/bin/env python

# script to run test-run of all three PGA's for 20 runs and save
# running time in a file.

import sys
import subprocess
import os
from datetime import datetime

ROOT_DIR = os.path.abspath(".")
BIN_DIR_NAME = "bin"
DATA_DIR_NAME = "stats_data"
OUTFILE_PREF = "Exec_times_x"

BIN_DIR = ROOT_DIR + "/" + BIN_DIR_NAME
DATA_DIR = ROOT_DIR + "/" + DATA_DIR_NAME
PGA1_F = "/pga_v2_P"
PGA2_F = "/pga_global_P"
PGA3_F = "/pga_partial_P"

RUNS = 20
POPSIZE = 5000


def main():
    data_file = DATA_DIR + "/" + OUTFILE_PREF + str(RUNS) + "_P"+ str(POPSIZE) + ".csv"
    
    with open(data_file, "a+") as f:
        datestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        f.write("# " + datestamp + " Population size: " + str(POPSIZE) + "\n")
        f.write("# Using: " + PGA1_F + str(POPSIZE)
                + " " + PGA2_F + str(POPSIZE) 
                + " and " + PGA3_F + str(POPSIZE) + " binary\n")
        f.write("# Nodes,sga time,pga_time,speedup\n")
        
        for i in range(0, RUNS):

            pga1_file =  BIN_DIR + "/" + PGA1_F + str(POPSIZE)
            pga2_file = BIN_DIR + "/"  + PGA2_F + str(POPSIZE)
            pga3_file = BIN_DIR + "/"  + PGA3_F + str(POPSIZE)

            if os.path.isfile(pga1_file) and os.path.isfile(pga2_file) and os.path.isfile(pga3_file):
                os.chdir(BIN_DIR)
                ps = subprocess.Popen("mpiexec -f ~/hosts -n 6 " + pga1_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                pga1_time = ps.stdout.readline()

                ps = subprocess.Popen("mpiexec -f ~/hosts -n 6 " + pga2_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                pga2_time = ps.stdout.readline()
                
                ps = subprocess.Popen("mpiexec -f ~/hosts -n 6 " + pga3_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                pga3_time = ps.stdout.readline()

                print "PGA time:", pga1_time
                print "PGA global time:", pga2_time
                print "PGA partial time:", pga3_time
                
                f.write(str(i) + "," + str(pga1_time).rstrip() + "," 
                    + str(pga2_time).rstrip() + "," + str(pga3_time).rstrip() 
                    + "\n")

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