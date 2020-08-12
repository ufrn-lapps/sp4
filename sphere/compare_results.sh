#!/bin/bash

# Clean and compile files
make clean
wait
make
wait

./propagation-reference > out-reference.txt
 wait
./propagation-sphere > out-sphere.txt
wait

diff out-reference.txt out-sphere.txt

