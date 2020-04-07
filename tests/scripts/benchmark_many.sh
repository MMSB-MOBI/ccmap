#!/bin/bash

PY_LIB=build/lib.macosx-10.9-x86_64-3.8
PY_SCRIPT=tests/scripts/threadsTest.py

if [ $# -lt 2 ] || [ $# -gt 3 ]
then
echo "Usage: benchmark.sh [cmap|lcmap|zmap|lzmap] POSES_SPECS.JSON PREFIX_TO_PDB_DIR[optional]"
exit 1
fi

TYPE=$1
SPECS_FILE=$2
PREFIX=""
if [ $# -eq 3 ];
then
    PREFIX=$3
fi

if ! [ -e "$SPECS_FILE" ]
then
    echo "$SPECS_FILE not found"
    exit 1
fi

FLAGS=""
if ! [ -z "$PREFIX" ]
    then
        if ! [ -d "$PREFIX" ]
        then
            echo "$PREFIX is not a valid folder"
            exit 1   
        fi
        FLAGS="$FLAGS --pref $PREFIX "
fi

if [ "$TYPE" == 'cmap' ] || [ "$TYPE" == 'lcmap' ];
    then
    FLAGS="$FLAGS --gen 10 "
fi    


for size in 100 500 1000 1500 2000 5000 15000 30000 50000 
    do
    for ncpu in $(seq 1 32)
        do
#        echo "python $PY_SCRIPT $TYPE $ncpu $size --lib $PY_LIB --inp $SPECS_FILE --encode --dist 5.0 $FLAGS"
              python $PY_SCRIPT $TYPE $ncpu $size --lib $PY_LIB --inp $SPECS_FILE --encode --dist 5.0 $FLAGS
                                                            # --inputs tests/1A2K_poses_specs_50K.json
    done |  awk -v size="$size" '$3 ~ /finished/ {print size,$1,$5}' 
done  