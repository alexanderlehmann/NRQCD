#!/bin/bash

processes=4
exename=./run.exe

mode=51
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

StartTime=0.
TimeRange=50.
TimeSpacing=0.1

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
Beta=16.

nMeasurement=10

FileName_WilsonLoops="loops.txt"
#FileName_TimeDerivativeWilsonLoops="tdivloops.txt"
#FileName_Potential="potential.txt"

#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $TimeRange $RandomSeed $Beta $nMeasurement $FileName_WilsonLoops #$FileName_TimeDerivativeWilsonLoops $FileName_Potential
