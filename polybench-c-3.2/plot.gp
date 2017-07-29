#! /usr/local/bin/gnuplot --persist

set terminal png enhanced font "Arial, 12" size 800, 800

set xtics 1
set ytics 1
set xrange [0:]
set yrange [0:]
set grid y x
set key inside b r
set xlabel "Threads"
set ylabel "Speedup"

set title "3MM"
set output "results/3mmMPI.png"

f(x)=x
plot f(x) dashtype 30 title 'Ideal', "results/3mm1.dat" using 1:2 title 'Pthread' with lines, "results/3mm1.dat" using 1:3 title 'OpenMP' with lines, "results/3mm1.dat" using 1:4 title 'OpenMPI com Comunicação' with lines, "results/3mm1.dat" using 1:5 title 'OpenMPI sem comunicação' with lines
