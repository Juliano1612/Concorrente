#! /usr/local/bin/gnuplot --persist

set terminal svg enhanced font "Helvetica, 16" size 800, 800

set xtics 1
set ytics $1
set grid y x
set key inside b r
set xlabel "Threads"
set ylabel "Speedup"

set title "$0 (size: $3)\n{/*0.7 $2}"
set output "results/$0.svg"

f(x)=x
plot f(x) dashtype 30 title 'Ideal', "results/$0.dat" using 1:2 title 'Pthread' with lines, "results/$0.dat" using 1:3 title 'OpenMP' with lines
