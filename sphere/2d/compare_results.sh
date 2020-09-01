#!/bin/bash

# Clean and compile files
make clean
wait
make
wait

wait
./propagation-sphere > out-sphere.txt
wait
diff out-reference.txt out-sphere.txt

