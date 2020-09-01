#!/bin/bash

# Clean and compile files
make clean
wait
make
wait

wait
./3d/propagation-sphere-3d > 3d/out-sphere-3d.txt
wait
./3d/propagation-reference-3d > 3d/out-reference-3d.txt
wait
diff 3d/out-reference-3d.txt 3d/out-sphere-3d.txt

