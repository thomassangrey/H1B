#!/bin/bash
#
# Use this shell script to compile (if necessary) your code and then execute it. Below is an example of what might be found in this file if your program was written in Python
#
cd  src/modules/ 
python setup.py build_ext --inplace
cd ../../
declare -r N="11"
head -${N} input/h1b_input.csv > input/year_2016_H1B.csv


python3.6 ./src/event_profiler.py 
