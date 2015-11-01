set term pos eps enh col
set out 'testbackprop_cubic.eps'
set multiplot layout 2,1
plot 'data/cubic.dat' u 1:2, 'testbackprop_cubic.out' u 2:3
set log y
plot 'testbackprop_cubic.error' u 1:2 w li
unset multiplot
