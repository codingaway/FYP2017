#!/bin/bash
#

# GLOBAL PARAMS for the GA programs
POPSIZE=50
MAXGENS=1000
NVARS=3
PXOVER=0.8
PMUTATION=0.15
IN_FILE=simple_ga_input.txt
OUT_FILE=simple_ga_output.txt

# Create Input file for Genes value range

# Run SGA with given parameters

# Run PGA with given parameters

# Save statistics into a file


~/bincpp/$ARCH/simple_ga > simple_ga_output.txt
#
if [ $? -ne 0 ]; then
  echo "Errors running simple_ga"
  exit
fi
#
echo "simple_ga read the input file simple_ga_input.txt, created simple_ga_output.txt"
