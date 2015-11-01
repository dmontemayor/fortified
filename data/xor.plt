set term pos eps enh col
set out 'testbackprop_xor.eps'
set multiplot layout 2,1
plot 'testbackprop_xor.out' u ($4-int($2+$3)%2) w imp ti 'xor error'
set log y
plot 'testbackprop_xor.error' u 1:2 w li
unset multiplot
