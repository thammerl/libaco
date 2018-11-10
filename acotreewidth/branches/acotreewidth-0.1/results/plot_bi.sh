#!/bin/bash

files=`ls *.plot | sed 's/\.sol\..\.plot//' | uniq`
for file in $files; do
  gnuplot << PLOT
  set term postscript eps enhanced
  set output "$file.sol.bi.eps"
  set xlabel "seconds"
  set ylabel "best-iteration treewidth"
  set style line 1 lt 1 lw 3
  set style line 2 lt 2 lw 3
  set style line 3 lt 3 lw 3
  set style line 4 lt 4 lw 3
  set style line 5 lt 6 lw 3
  plot "$file.sol.1.plot" using 2:4 title "$file.sol.1.plot" with lines linestyle 1,\
  "$file.sol.2.plot" using 2:4 title "$file.sol.2.plot" with lines linestyle 2,\
  "$file.sol.3.plot" using 2:4 title "$file.sol.3.plot" with lines linestyle 3,\
  "$file.sol.4.plot" using 2:4 title "$file.sol.4.plot" with lines linestyle 4,\
  "$file.sol.5.plot" using 2:4 title "$file.sol.5.plot" with lines linestyle 5
PLOT
done
