#!/gnuplot
set terminal pdf color enhanced font "Arial,22" size 12,8
set grid
set xlabel "Time [s]"
set ylabel "Potential [V]"
set format x "%2.1e"
set format y "%5.3e"
set tmargin 1.5

set output 'object_potential_history.pdf'
plot 'data/object_Spacecraft_0.txt' u 2:3 w l lw 3 ti "potential"
unset output

set ylabel "Current [A]"
set output 'object_current_history.pdf'
plot 'data/object_Spacecraft_0.txt' u 2:5 w l lw 3 ti "electron" ,\
'' u 2:6 w l lw 3 ti "proton"
unset output