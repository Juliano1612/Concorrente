#! /usr/local/bin/gnuplot --persist

set terminal png enhanced font "Arial, 12" size 800, 800

set xtics 1
set ytics 10000000000
set xrange [1:]
set yrange [150000000000:]
set grid y x
set key inside b r
set xlabel "Threads"
set ylabel "Speedup"

set title "Cache miss L1 (size: 1024)"
set output "results/3mmCML1.png"

f(x)=x
plot "results/3mmCacheL1.dat" using 1:2 title 'Sequencial' with lines, "results/3mmCacheL1.dat" using 1:3 title 'Pthread' with lines, "results/3mmCacheL1.dat" using 1:4 title 'OpenMP' with lines

Pthread2	
Pthread4	
Pthread8	
Pthread16	
