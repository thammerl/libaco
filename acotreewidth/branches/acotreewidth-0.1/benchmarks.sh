#!/bin/bash

BENCHMARKS=( '-f benchmarks/games120.col -a 2 -b 20 -e 0.05 -m 20' 'results/games120.sol'
             '-f benchmarks/anna.col -a 2 -b 20 -e 0.05 -m 20' 'results/anna.sol' 
             '-f benchmarks/DSJC125.9.col -a 2 -b 20 -e 0.05 -m 20' 'results/DSJC125.9.sol' 
             '-f benchmarks/le450_5a.col -a 2 -b 20 -e 0.05 -m 5' 'results/le450_5a.sol' 
             '-f benchmarks/le450_5b.col -a 2 -b 20 -e 0.05 -m 5' 'results/le450_5b.sol' 
             '-f benchmarks/miles500.col -a 2 -b 20 -e 0.05 -m 20' 'results/miles500.sol' 
             '-f benchmarks/miles750.col -a 2 -b 20 -e 0.05 -m 20' 'results/miles750.sol' 
             '-f benchmarks/myciel7.col -a 2 -b 20 -e 0.05 -m 20' 'results/myciel7.sol' 
             '-f benchmarks/queen8_8.col -a 2 -b 20 -e 0.05 -m 20' 'results/queen8_8.sol' 
             '-f benchmarks/school1.col -a 2 -b 20 -e 0.05 -m 5' 'results/school1.sol' )

for (( i=0;i<${#BENCHMARKS[@]};i=i+2)); do
  ./acotreewidth ${BENCHMARKS[$i]} -i 10000 > ${BENCHMARKS[$i+1]}.1
  ./acotreewidth ${BENCHMARKS[$i]} -i 10000 > ${BENCHMARKS[$i+1]}.2
  ./acotreewidth ${BENCHMARKS[$i]} -i 10000 > ${BENCHMARKS[$i+1]}.3
  ./acotreewidth ${BENCHMARKS[$i]} -i 10000 > ${BENCHMARKS[$i+1]}.4
  ./acotreewidth ${BENCHMARKS[$i]} -i 10000 > ${BENCHMARKS[$i+1]}.5
done
