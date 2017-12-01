#!/gnuplot
set terminal pdf color enhanced font "Arial,30" size 12,8
set grid
set xlabel "Time [s]"
set ylabel "Force [N]"
set format x "%2.1e"
set format y "%5.3e"
set tmargin 1.5

set output 'object_force_history.pdf'
plot 'data/object_Spacecraft.txt' u 2:7 w l lw 3 ti "F_x"
unset output