#! /usr/local/bin/gnuplot --persist

set terminal svg enhanced font "Helvetica, 16" size 800,800

set xtics 1
set ytics .1
set grid y x
set key inside b r
set xlabel "Threads"
set ylabel "Speedup"

set title "RBSOR (size: 2048)\n{/*0.7 Intel(R) Core(TM) i5-2500 CPU \\\@ 3.30GHz}"
set output "results/rbsor.svg"

f(x)=x
plot f(x) dashtype 30 title 'Ideal', "results/rbsor.dat" using 1:2 title 'OpenMP' with linespoints, "results/rbsor.dat" using 1:3 title 'PThreads' with linespoints

