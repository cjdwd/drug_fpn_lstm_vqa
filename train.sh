#!/bin/bash
trainfile="train.py"
if [ "$1" == "kd" ]
then
logfile="./logs/kd.out"

if [ "$1" == "kiba" ]
then
logfile="./logs/kiba.out"

if [ "$1" == "dude" ]
then
logfile="./logs/dude.out"

nohup python -u ${trainfile} $1 > ${logfile} 2>&1 &

fi


