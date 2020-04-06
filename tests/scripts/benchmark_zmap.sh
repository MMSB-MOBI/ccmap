#!/bin/bash



for size in 100 500 1000 1500 2000 5000 15000 30000 500000 
    do
    for ncpu in $(seq 1 32)
        do
        python tests/scripts/threadsTest.py zmap $ncpu $size --lib build/lib.macosx-10.9-x86_64-3.8 \
                                                             --pdbI tests/1A2K_r_u.pdb \
                                                             --pdbJ tests/1A2K_l_u.pdb \
                                                             --encode --dist 5.0 \
                                                             --inputs tests/1A2K_poses_specs_50K.json
    done |  awk -v size="$size" '$4 ~ /finished/ {print size,$1,$6}' 
done  