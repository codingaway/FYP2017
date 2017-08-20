#!/bin/sh
# 
# File:   running_time.sh
# Author: Abdul Halim <13029096@studentmail.ul.ie>
#
# Created on 08-Mar-2017, 14:00:44
#
# get time in nano seconds
start=`date +%s%N`
end=`date +%s%N`
runtime=$((end-start))
./sga > /dev/null 2>&1
echo "Start time : $start"
echo "End Time   : $end" 
echo $runtime