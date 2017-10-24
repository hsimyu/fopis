#!/gnuplot
set terminal pdf color enhanced font "Arial,22" size 12,8
set grid
set xlabel "Time [s]"
set format x "%2.1e"
set format y "%5.3e"
set tmargin 1.5

set ylabel "Current [A]"
set output 'test_current_history.pdf'
plot 'data/object_Spacecraft.txt' u 2:5 w l lw 3 ti "electron",\
'' u 2:6 w l lw 3 ti "proton"
unset output

set ylabel "Energy [J]"
set output 'test_energy_history.pdf'
plot 'data/energy.txt' u 2:3 w l lw 3 ti "Total",\
'' u 2:4 w l lw 3 ti "Particle",\
'' u 2:5 w l lw 3 ti "Efield"
unset output

set ylabel "Particle Number"
set output 'test_pnum_history.pdf'
plot 'data/total_particle_number.txt' u 2:3 w l lw 3 ti "Total",\
'' u 2:4 w l lw 3 ti "Electron",\
'' u 2:5 w l lw 3 ti "Proton"
unset output