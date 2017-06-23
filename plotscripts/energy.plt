#!/gnuplot
set terminal pdf color enhanced font "Arial,22" size 12,8
set grid
set xlabel "Time [s]"
set ylabel "Energy [J]"
set format x "%2.1e"
set format y "%5.3e"
set tmargin 1.5

set output 'energy_history.pdf'
plot 'data/energy.txt' u 2:3 w l lw 3 ti "total",\
'data/energy_with_inject.txt' u 2:3 w l lw 3 ti "w/ inject, total"
unset output

system('open ./')
