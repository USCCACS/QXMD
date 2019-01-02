set term png
set output 'eng.png'
set key box opaque
set ylabel "Energy [eV]"
set xlabel "Step Number"

set title "Total Energy Vs. Step Number"

plot "../../data/md_eng.d" using 1:2 title " Energy" w l lc rgb "#000000", \
