set term png
set output 'eig.png'
set key box opaque
set xrange [0:0.12]
set yrange [-2.0:2.0]
set ylabel "Energy [eV]"
set xlabel "Time [ps]"

set title "Kohn-Sham Energy Eigenvalues Vs. Time"

plot "EIG.dat" using ($1*0.0012):2 title "  Unoccupied" pt 7 ps 0.75 lc rgb "#000000", \
"EIG_occ-one.dat" using ($1*0.0012):2 title "  Singly Occupied" pt 7 ps 0.75 lc rgb "#FF0000", \
"EIG_occ-two.dat"  using ($1*0.0012):2 title "  Doubly Occupied" pt 7 ps 0.75 lc rgb "#00FF00" 
