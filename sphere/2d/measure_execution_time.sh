#!/bin/bash

# Clean and compile files
make clean
wait
make
wait


for i in {1..5}
do
  echo "Execution $i"
  echo "Reference" 
  $(/usr/bin/time -f'%e' ./propagation-reference)
  wait
  echo "Sphere" 
  $(/usr/bin/time -f'%e' ./propagation-sphere)
  wait
done

