#!/usr/bin/env python

# Python script to generate Run test programs, collect execution time, create data file for gnuplot,
# generate plot and save them.

import sys
import subprocess
import os

ROOT_DIR = os.path.abspath(".")
BIN_DIR_NAME = "bin"
DATA_DIR_NAME = "stats_data"
GRAPH_DIR_NAME = "graphs"

BIN_DIR = ROOT_DIR + "/" + BIN_DIR_NAME
DATA_DIR = ROOT_DIR + "/" + DATA_DIR_NAME
GRAPH_DIR = ROOT_DIR + "/" + GRAPH_DIR_NAME

POPLIST = [100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]


def main():
    data_file = DATA_DIR + "/popsize_vs_execution.data"
    
    with open(data_file, "a+") as f:
        f.write("POPSIZE | SGA | PGA | PGA_BCAST\n")
        
        for popsize in POPLIST:

            sga_file =  BIN_DIR + "/sga_P" + str(popsize)
            pga_file = BIN_DIR + "/pga_P" + str(popsize)
            pgab_file = BIN_DIR + "/pga_bcast_P" + str(popsize)

            if os.path.isfile(sga_file) and os.path.isfile(pga_file):
                os.chdir(BIN_DIR)
                ps = subprocess.Popen(sga_file, stdout=subprocess.PIPE)
                ps.wait()
                sga_time = ps.stdout.readline()

                ps = subprocess.Popen("mpiexec -f ~/hosts -n 5 " + pga_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                pga_time = ps.stdout.readline()
                
                ps = subprocess.Popen("mpiexec -f ~/hosts -n 5 " + pga_file, stdout=subprocess.PIPE, shell=True)
                ps.wait()
                pgab_time = ps.stdout.readline()

                print "SGA time:", sga_time
                print "PGA time:", pga_time
                print "PGA_B time:", pgab_time
                
                f.write(str(popsize) + " " + str(sga_time).rstrip() + " " 
                    + str(pga_time).rstrip() + " " + str(pgab_time).rstrip() 
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