set term pos eps enh col
set out 'testbackprop_xorpoisson.eps'
set multiplot layout 2,1
plot 'testbackprop_xorpoisson.out' u ($4-int($2+$3)%2) w imp ti 'xorpoisson error'
set log y
plot 'testbackprop_xorpoisson.error' u 1:2 w li
unset multiplot
