#!/bin/bash
# Launches init.hoc with increased stack size to enable generation of model neurons with large, myelinated axons using recursive functions
# AUTHOR: Aman Aberra
# CONTACT: aman.aberra@duke.edu
sim_dir="$(dirname "$0")"
cd $sim_dir # ensures .mod files are loaded properly

special_file="${sim_dir}/x86_64/special"
#nrngui -NSTACK 100000 -NFRAME 20000 init.hoc 
if [ -f ${special_file} ]; then
    #./${special_file} -NSTACK 100000 -NFRAME 20000 init.hoc
    nrngui -NSTACK 100000 -NFRAME 20000 init.hoc 
else
    echo "Compiling mod files"
    nrnivmodl mechanisms
    #./${special_file} -NSTACK 100000 -NFRAME 20000 init.hoc
    nrngui -NSTACK 100000 -NFRAME 20000 init.hoc 
fi