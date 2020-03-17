#!/bin/bash

processes=4
exename=./run.exe

mode=53
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

StartTime=0.
TimeRange=30.
TimeSpacing=0.1
Rmax=6

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
Beta=16.
nefieldinit=10
nequilibriumtimesteps=300
FileName_WilsonLoops="wilsonloops.txt"
MeasureEnergy=1
FileName_Energy="thermalization_energy.txt"

#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $TimeRange $Rmax $RandomSeed $Beta $nefieldinit $nequilibriumtimesteps $FileName_WilsonLoops $MeasureEnergy $FileName_Energy
