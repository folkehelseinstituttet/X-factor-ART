#!/bin/bash

file_in="DATA/yin_science_2017_TF_types.txt"
file_out="DATA/extracted_lines.txt"
file_out2="DATA/extracted_lines_clean.txt"

# get only the lines with the names of transcription factors and their classes
# then, remove all the extraneous whitespace in the beginning of lines
# and replace many whitespace characters with new line
# then, I can read it into R and create a nice tibble
grep " : " $file_in | sed -re 's/^\s+//g' | sed -re 's/  +/\n/g' > $file_out

# check whether all lines have two columns
in2csv $file_out -d ":" -f csv | csvlook | trim 30

# there are some that have only one column - are these all the same?
# (i.e., containing only 'Index'?)
grep -v ":" -n $file_out

grep -v ":" -c $file_out
# YES - these are all the same and I can remove those
grep ":" $file_out > $file_out2

