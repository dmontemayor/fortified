set term pos eps enh col
set out 'testbackprop_xorseq.eps'
set multiplot layout 2,1
plot 'data/xorseq.dat' u 2, 'testbackprop_xorseq.out' u 3
set log y
plot 'testbackprop_xorseq.error' u 1:2 w li
unset multiplot
