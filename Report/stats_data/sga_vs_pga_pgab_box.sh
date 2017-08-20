#! /bin/sh
# gnuplot script file for ploting data in file "sga_vs_pga_pgab_xPopSize.dat"
# Plot graph for execution time SGA vs. PGA vs. PGA_BCAST
filename=popsize_vs_execution_new
datafile="${filename}.dat"
# outfile="${filename}.tex"
outfile="graphs/${filename}_box.png"
# title="Speed up comparison using population size: 5000"
xlabel="Population size"
ylabel="Execution time(sec)"
gnuplot << EOF
# plot "${datafile}" using 1:2 title "SGA" with lines, "${datafile}" using 1:3 title "PGA"  with lines, "${datafile}" using 1:4 title "PGA\_Bcast" noenhanced  with lines
reset
set terminal pngcairo enhanced font 'Verdana,10'
# set terminal qt font 'Verdana,10'
set output "${outfile}"

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
# set border 3 back ls 11
set tics nomirror

# Thinner, filled bars
set boxwidth 0.8
set style fill solid 1.00 

# define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

# Rotate X labels and get rid of the small stripes at the top (nomirror)
set xtics nomirror rotate by -45

# Replace small stripes on the Y-axis with a horizontal gridlines
set tic scale 0
set grid ytics lc rgb "#505050"
# Remove border around chart
unset border

# color definitions
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 3 lc rgb '#365e9c' pt 6 ps 1 lt 1 lw 2 # --- green

set key top left
# set title "${title}"
set xlabel "${xlabel}"
set ylabel "${ylabel}"
set style data histograms

# set xrange [0:1]
# set yrange [0:1]

plot "${datafile}" using 2:xtic(1) title "SGA" lt rgb "#406090",\
     "${datafile}" using 3 title "PGA" lt rgb "#40FF00", \
     "${datafile}" using 4 title "PGA\_Bcast" noenhanced lt rgb "#FF4060"
set output
EOF