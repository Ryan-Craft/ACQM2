#!/bin/bash


energies=($(seq 5.0 2.0 9.0))

echo ${energies[*]}
for i in ${!energies[@]}
do 
   
  mkdir "energy_${energies[$i]}"
  sed s/EEEE/${energies[$i]}/g data.in > "energy_${energies[$i]}"/data.in
  cd "energy_${energies[$i]}"  
  ../main
  cd ..

done
