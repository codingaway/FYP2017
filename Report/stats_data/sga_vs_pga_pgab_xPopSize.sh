#! /bin/sh
# gnuplot script file for ploting data in file "sga_vs_pga_pgab_xPopSize.dat"
# Plot graph for execution time SGA vs. PGA vs. PGA_BCAST
filename=popsize_vs_execution_new
datafile="${filename}.dat"
# outfile="${filename}.tex"
outfile="graphs/${filename}.png"
# title="Speed up comparison using population size: 5000"
xlabel="Population size"
ylabel="Execution time(sec)"
gnuplot << EOF
# set terminal postscript
# set output "${outfile}"
# set autoscale
# unset log
# unset label
# set xtic auto
# set ytic auto

# plot "${datafile}" using 1:2 title "SGA" with lines, "${datafile}" using 1:3 title "PGA"  with lines, "${datafile}" using 1:4 title "PGA\_Bcast" noenhanced  with lines
reset
# wxt
# set terminal svg enhanced font 'Verdana,9'
# set terminal epslatex color font 'Verdana,9'
set terminal pngcairo enhanced font 'Verdana,10'
set output "${outfile}"
# png
#set terminal pngcairo size 410,250 enhanced font 'Verdana,9'
#set output 'nice_web_plot.png'
# svg
#set terminal svg size 410,250 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '9' rounded dashed
#set output 'nice_web_plot.svg'

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror
# define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

# color definitions
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 3 lc rgb '#365e9c' pt 6 ps 1 lt 1 lw 2 # --- green

set key bottom right
# set title "${title}"
set xlabel "${xlabel}"
set ylabel "${ylabel}"

# set xrange [0:1]
# set yrange [0:1]

# plot 'nice_web_plot.dat' u 1:2 t 'Example line' w lp ls 1, \
#      ''                  u 1:3 t 'Another example' w lp ls 2

plot "${datafile}" using 1:2 title "SGA" w lp ls 1, \
	"${datafile}" using 1:3 title "PGA"  w lp ls 2, \
	"${datafile}" using 1:4 title "PGA\_Bcast" noenhanced w lp ls 3
set output
EOF