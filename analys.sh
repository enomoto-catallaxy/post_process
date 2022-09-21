#!/bin/sh

### \lambda = 5 ###
for k in 2 4 8 12;
  do
  echo $k
  sed -i "28c\ \t char \t fname[256]=\"osci_sphere_channel_D40umL40umRe02Ca${k}lam5_freq01\";" ./main.c
  make clean;make;./run
done

#sed -i "28c\ \t char \t fname[256]=\"osci_single_rbc_rectangle_ca01gs4e-06k1_2e-19lam5_freq3_2\";" ./main.c
#make clean;make;./run
