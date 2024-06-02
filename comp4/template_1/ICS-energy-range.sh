#!/bin/bash


energies=($(seq 5 2.5 50))

echo ${energies[*]}
for i in ${!energies[@]}
do  
  cd "energiesfolder"
  mkdir energy=energy_${energies[$i]}
  cd "energy_${energies[$i]}"  
  echo ${energies[$i]} >> proout_${energies[$i]}.txt
  cd ..

done
