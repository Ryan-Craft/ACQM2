#!/bin/bash

mkdir "energiesplot"
rm -r collated_ICS.txt
cd energiesplot
energies=($(seq 0 0.5 50))

rm -r energy_*
for i in ${!energies[@]}
do 
  arr=() 
  mkdir "energy_${energies[$i]}"
  sed s/EEEE/${energies[$i]}/g ../data.in > "energy_${energies[$i]}"/data.in
  cd "energy_${energies[$i]}"  
  ../../main
  file=$(ls ICSout.txt)

  file=$(ls ICSout.txt) 
  allvals=($(cat $file | tr ',\n' ' '))
  string_print=$(echo "${energies[$i]} ${allvals[0]} ${allvals[2]} ${allvals[4]} ${allvals[6]} ${allvals[1]} ${allvals[3]} ${allvals[5]} ${allvals[7]}")
  echo $string_print >> ../../collated_ICS.txt 

  cd ..

done
