#!/bin/bash

processes=1
exename=./run.exe

mode=4

nheaderlines=0
filename_input="correlator_3s1.txt"
filename_output="spectrum_3s1.txt"

mpiexec -n $processes \
	$exename \
	$mode $nheaderlines $filename_input $filename_output
