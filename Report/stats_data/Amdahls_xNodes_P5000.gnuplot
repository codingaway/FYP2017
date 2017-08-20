# gnuplot script file for ploting data in file "sga_vs_pga_pgab_xPopSize.dat"
# Plot graph for execution time SGA vs. PGA vs. PGA_BCAST
filename

#! /bin/sh
name=$1
gnuplot << EOF
set terminal png
set output "graphs/Amdahls_xNodes_P5000.png"
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Graph title"
set xlabel "Population size"
set ylabel "Execution time"
plot "sga_vs_pga_pgab_xPopSize.dat" using 1:2 title "SGA" with lines, "sga_vs_pga_pgab_xPopSize.dat" using 1:3 title "PGA"  with lines, "sga_vs_pga_pgab_xPopSize.dat" using 1:4 title "PGA\_Bcast" noenhanced  with lines
replot
EOF