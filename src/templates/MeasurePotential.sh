#!/bin/bash

processes=4
exename=./run.exe

mode=8
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

StartTime=0.
TimeRange=1.
TimeSpacing=0.1

# Initialisation of gluon field
RandomSeed=2
# Given in units of az
InitialGluonSaturationScale=1.2
InitialGluonOccupationAmplitude=1.0
Coupling=0.01

FileName_WilsonLoops="wilsonloops.txt"
FileName_PotentialWilsonLoops="potentialwilsonloops.txt"


#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $TimeRange $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $FileName_WilsonLoops $FileName_PotentialWilsonLoops

