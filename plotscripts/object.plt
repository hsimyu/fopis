#!/gnuplot
set terminal pdf color enhanced font "Arial,22" size 12,8
set grid
set xlabel "Time [s]"
set ylabel "Potential [V]"
set format x "%2.1e"
set format y "%5.3e"
set tmargin 1.5

set output 'opject_potential_history.pdf'
plot 'data/object_Sp0.txt' u 2:3 w l lw 3 ti "potential",\
'' u 2:5 w l lw 3 ti "electron" axes x1y2 ,\
'' u 2:6 w l lw 3 ti "proton" axes x1y2
unset output
