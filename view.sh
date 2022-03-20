#!/bin/bash
n="1000"
if [ "$1" == "kd" ]
then
logfile="./logs/kd.out"
tail -f -n $n ${logfile}
fi


