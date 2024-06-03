#!/bin/bash

file=$(ls ICSout.txt)

IFS="\n"
for line in $(cat "$file")
do
  echo $line
  IFS=" "
  set -o noglob
  arr=($line)
done
echo 
echo ${arr[0]}
string_print=$(echo "${arr[0]} ${arr[2]} ${arr[4]} ${arr[6]} ${arr[1]} ${arr[3]} ${arr[5]} ${arr[7]}" | tr -d "\n")
echo $string_print >> collatet.txt
