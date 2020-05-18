#!/bin/sh
echo This shell script copies the *dhdl.xvg file from each folder to a newly made folder called dhdl_files and rename all the *dhdl.xvg files.

set -e     # exit upon error
read -p "Please input the common prefix of the files: " p
read -p "Please input the number of replicas: " N
read -p "Please input the number of *dhdl.xvg files in each folder: " n

mkdir dhdl_files

echo Copying/renaming the *dhdl.xvg files ...

for (( i=0; i<$N; i=i+1 ))
do
    cp state_${i}/*dhdl.xvg dhdl_files/${p}_dhdl_${i}.xvg

    for (( j=2; j<$n+1; j=j+1 ))
    do
        cp state_${i}/*0${j}.xvg dhdl_files/${p}_dhdl_${i}_part${j}.xvg
    done
done

echo Complete!