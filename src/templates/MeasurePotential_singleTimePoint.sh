#!/bin/bash

processes=4
exename=./run.exe

mode=9
Nx=24
Ny=24
Nz=24
ax=1.0
ay=1.0
az=1.0

StartTime=0.
TimeSpacing=0.1

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
InitialGluonSaturationScale=1.2
InitialGluonOccupationAmplitude=1.0
Coupling=0.01

FileName_WilsonLoops="wilsonloops_24x24x24_t0.txt"
FileName_PotentialWilsonLoops="potentialwilsonloops_24x24x24_t0.txt"


#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $FileName_WilsonLoops $FileName_PotentialWilsonLoops

