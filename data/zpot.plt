set term pos eps enh col
set out 'testbackprop_zpot.eps'
set multiplot layout 2,1
plot 'data/zpot.dat' u 1:2, 'testbackprop_zpot.out' u 2:3
set log y
plot 'testbackprop_zpot.error' u 1:2 w li
unset multiplot
