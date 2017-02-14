#!/gnuplot
set terminal pdf color enhanced font "Arial,22" size 12,8
set grid
set xlabel "Energy [eV]"
set ylabel "Ratio"
set format x "%2.0f"
set tmargin 1.5
me = 9.1093818872e-31
e = 1.6021773e-19
velo_to_ev(v) = me*(v**2)*0.5/e
ev_to_velo(ev) = sqrt(2.0*ev*e/me)
ev_to_velo2(ev) = 2.0*ev*e/me
maxwellian(x, kbte) = (4.0 * pi) * (x**2) * ((me/(2.0*pi*kbte*e))**1.5) * exp(-1.0 * (me/(2.0*kbte*e)) * (x**2))

set output 'energy_distribution_electron.pdf'
plot \
'data/energy_distribution_0000_0000.csv' ind 0 u 1:2 w boxes lw 3 ti "electron, t = 0"
unset output

set output 'energy_distribution_proton.pdf'
plot \
'data/energy_distribution_0000_0000.csv' ind 1 u 1:2 w boxes lw 3 ti "proton, t = 0"
unset output

# example plot maxwellian(x*1e3, 1.0)/maxwellian(ev_to_velo(1.0), 1.0) lw 3
set xlabel "Velocity [km/s]"
set output 'velocity_distribution_electron.pdf'
plot \
'data/velocity_distribution_0000_0000.csv' ind 0 u 1:2 w boxes lw 3 ti "electron, t = 0",\
maxwellian(x*1e3, 1.0)/maxwellian(ev_to_velo(1.0), 1.0) lw 3
unset output

set output 'velocity_distribution_proton.pdf'
plot \
'data/velocity_distribution_0000_0000.csv' ind 1 u 1:2 w boxes lw 3 ti "proton, t = 0"
unset output
