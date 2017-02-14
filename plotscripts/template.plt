#!/gnuplot
set terminal pdf color enhanced font "Arial,22" size 12,8
set output 'test.pdf'
set grid
set xlabel ""
set ylabel ""
set format x "%2.0f"
set tmargin 1.5
plot '' u 1:2 w l lw 3 ti ""
unset output
