set term pos eps enh col
set out 'testbackprop_xorsinusoid.eps'
set multiplot layout 2,1
plot 'testbackprop_xorsinusoid.out' u ($4-int($2+$3)%2) w imp ti 'xorsinusoid error'
set log y
plot 'testbackprop_xorsinusoid.error' u 1:2 w li
unset multiplot
