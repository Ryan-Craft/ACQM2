#!/bin/bash

#rm -r dx=*
for dx in 0.5 0.1 0.01
do
    mkdir dx=${dx}        #make a directory for that value of R
    sed s/dddd/${dx}/ NumerovsQHOParams.txt > dx=${dx}/NumerovsQHOParams.txt
    #replace the string "nnnn" in the data file with the value of R from the list

    cd dx=${dx} #go to the directory for this specific R
    #../NUmerovsQHO
    #run the code ,which will use the parameters in this folder to generate an output file in this directory
        n=`awk '{print $NF}' NumerovsQHOParams.txt`
        echo $n 
        ../NumerovsQHO > Numerov_${dx}_${n}.data 
        mv psi.txt psi_${n}.txt


    cd ..

done
