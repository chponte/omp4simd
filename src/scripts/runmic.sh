#!/bin/bash

#*************************** Arguments ***************************
# $1: Xeon Phi address
# $2: Path to folder containing all programs
# $3: Path to results folder

if [ $# -lt 3 ]; then
    echo "Missing args"
    echo "Usage: runphi.sh <address> <path to program folder> \
            <path to results folder> [label]"
    exit 0
fi

#************************** Definitions **************************
# MIC libraries path
LIB_PATH="~/compilers/intel/lib/mic"
# Task script
TASK="~/sources/omp4simd/src/scripts/task.sh"
# Log file name
LOGFILE="info.log"
# Unique name for multiple purposes
uni=$(date "+%d_%b_%Y_%T")

#*****************************************************************

# Check if results directory exists
if ! [ -d $3 ]; then
    mkdir -p $3
fi

# Upload files into Xeon Phi
ssh $1 "mkdir -p ~/bins/$uni && mkdir -p ~/results/$uni" && \
    scp $(find $2) $1:~/bins/$uni/ && scp $TASK $1:~/bins/$uni/

# Check if any library is missing
locallibs=ls $LIB_PATH | sort | tr -d '\n'
remotelibs=ssh $1 'ls ~/libs/mic' | sort | tr -d '\n'
if [ "x$(echo $locallibs$remotelibs | grep -e "^\(.*\)\1$")" = "x" ]; then
    ssh $1 "rm -rf ~/libs/mic; mkdir ~/libs" > /dev/null 2>&1
    scp -r $LIB_PATH $1:~/libs/mic
fi

# # Check if any compiled math function is missing
# comp=$(echo "$(ls $ST_PATH/util | grep -e "^.*\.c$" | sort | tr -d '.c' | tr -d '\n')$(ssh $1 'ls ~/libs/util' 2>/dev/null | sort | tr -d '\n')" | grep -e "^\(.*\)\1$")
# if [ "x$comp" = "x" ]; then
#     ssh $1 "rm -rf ~/libs/util; mkdir -p ~/libs/util/" > /dev/null 2>&1
#     cd $ST_PATH && make util_func
#     cd util && scp $(ls | grep -v "Makefile" | grep -v -e "^.*\.c$") $1:~/libs/util/ && make clean
# fi

# Print info associated with current execution
if [ "x$3" != "x" ]; then
    echo "[$1] $uni: $4" >> $3/$LOGFILE
else
    echo "[$1] $uni" >> $3/$LOGFILE
fi
echo "$(ls $2 | tr '\n' ' ')" >> $3/$LOGFILE

# Run codes and extract results
screen -S $uni -d -m ssh $1 "export LD_LIBRARY_PATH=~/libs/mic && \
    ~/bins/run.sh $uni ~/bins/$uni 10 ~/results/$uni \"$3\"" \
    && scp -r $1:~/results/$uni $RES_PATH
