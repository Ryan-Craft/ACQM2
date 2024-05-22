#!/bin/bash
ls PEC.* && rm PEC.* 
rm -r R=*
for R in 0.01 0.05 0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0
do  
    mkdir R=${R}        #make a directory for that value of R
    sed s/RRRR/${R}/ LaguerreParams.txt > R=${R}/LaguerreParams.txt   #replace the string "RRRR" in the data file with the value of R from the list
    cd R=${R} #go to the directory for this specific R
    ../main #run the code ,which will use the parameters in this folder to generate an output file in this directory

    for i in $(seq 1 4) #loop over the sequence 1,2,3,4
    do
        echo ${R} $(awk '$1=='${i} E.out) >> ../PEC.${i} 
        #idk wtf this does
    done
    cd ..
    
done