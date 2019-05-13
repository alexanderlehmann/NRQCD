#!/bin/bash

processes=1
exename=./run.exe

mode=4

nheaderlines=1
#filename_input="/home/lehmann/results/nrqcd/free_spectrum/O_m2/2019-04-01/32x32x32_long_dt1/correlator_3s1.txt"
#filename_output="/home/lehmann/results/nrqcd/free_spectrum/O_m2/2019-04-01/32x32x32_long_dt1/spectrum_3s1.txt"
filename_input=mesoncorrelator_3s1.txt
filename_output=spectrum_3s1.txt
tmax=10.


mpiexec -n $processes \
	$exename \
	$mode $nheaderlines $filename_input $filename_output $tmax
