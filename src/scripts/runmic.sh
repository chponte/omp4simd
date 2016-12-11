#!/bin/bash

#*************************** Licensing ***************************
# This code is distributed under the GNU GPLv3 license. 
# Author: Christian Ponte Fernández
# Modified: 10 November 2016

#*************************** Arguments ***************************
# $1: Xeon Phi address
# $2: Path to folder containing all programs
# $3: Path to results folder
# $4: Repetitions
# $5: Label [Optional]

if [ $# -lt 4 ]; then
    echo "Missing args"
    echo "Usage: $0 <address> <path to program folder> \
<path to results folder> <repetitions per program> [label]"
    exit 0
fi

#************************** Definitions **************************
# MIC libraries path
LIB_PATH=~/compilers/intel/lib/mic
# Task script
TASK=~/sources/omp4simd/src/scripts/task.sh
# Log file name
LOGFILE=info.log
# Unique name for multiple purposes
uni=$(date "+%d_%b_%Y_%T")

# Args:
#   $1: Values separated by spaces
function mean {
    _sum=0
    for value in $1; do
        _sum=$(echo "$_sum+$value" | bc -l)
    done
    echo "$_sum/$(wc -w <<< $1)" | bc -l
}

# Args:
#   $1: Values separated by spaces
#   $2: Mean
function standard_deviation {
    local _sum=0
    for value in $1; do
        _sum=$(echo "$_sum + ($value-$2)^2" | bc -l)
    done
    echo "sqrt($_sum/$(wc -w <<< $1))" | bc -l
}

# Args:
#   $1: Results file path
#   $2: Label
function results_to_csv {

    # Write label at the beginning of the csv file
    echo "$2,,,,,"
    echo ",,,,,"

    echo "Executable,Thread number,Avg. runtime,S.d.,Runtimes,"
    while IFS='' read -r line || [[ -n "$line" ]]; do
        _e=$(echo $line | cut -d' ' -f1 | sed -e "s/.*\/\([^\/]*\)$/\1/g")
        _tn=$(echo $line | cut -d' ' -f2)
        _rts=$(echo $line | cut -d' ' -f'1 2' --complement)
        _m=$(mean "$_rts")
        _sd=$(standard_deviation "$_rts" "$_m")
        echo "$_e,$_tn,$_m,$_sd,$_rts,"
    done < $1
    unset IFS
}

#*****************************************************************

# Check if results directory exists
if ! [ -d $3/$uni ]; then
    mkdir -p $3/$uni
fi

# Upload files into Xeon Phi
findcount=$(find $2 | wc -l)
ssh $1 "mkdir -p ~/bins/$uni && mkdir -p ~/results/$uni" && \
    scp $(find $2 | tail -n $((findcount-1))) $1:~/bins/$uni/ && \
    scp $TASK $1:~/bins/$uni/

# Check if any library is missing
locallibs=$(ls $LIB_PATH | sort | tr -d '\n')
remotelibs=$(ssh $1 'ls ~/libs/mic' 2>/dev/null | sort | tr -d '\n')
if [ "x$(echo $locallibs$remotelibs | grep -e "^\(.*\)\1$")" = "x" ]; then
    ssh $1 "rm -rf ~/libs/mic; mkdir ~/libs" > /dev/null 2>&1
    scp -r $LIB_PATH $1:~/libs/mic
fi

# Print info associated with current execution
if [ "x$3" != "x" ]; then
    echo "[$1] $uni: $5" >> $3/$LOGFILE
else
    echo "[$1] $uni" >> $3/$LOGFILE
fi
echo "$(find $2 | tail -n $((findcount-1)) | tr '\n' ' ')" >> $3/$LOGFILE

# Export auxiliary functions to screen subshell
export -f mean standard_deviation results_to_csv

# Run codes and extract results
screen -S "runmic.sh $uni" -d -m bash -c \
    "ssh $1 \"export LD_LIBRARY_PATH=~/libs/mic && \
        ~/bins/$uni/task.sh $uni ~/bins/$uni $4 ~/results/$uni results.out\" && \
    scp -r $1:~/results/$uni $3 && \
    results_to_csv $3/$uni/results.out \"$5\" > $3/$uni/results.csv"
