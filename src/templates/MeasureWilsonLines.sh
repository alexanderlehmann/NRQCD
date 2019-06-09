#!/bin/bash

processes=4
exename=./run.exe

mode=5
Nx=8
Ny=8
Nz=8
ax=1.0
ay=1.0
az=1.0

StartTime=0.
TimeRange=10.
TimeSpacing=0.1
Rmax=6

# Initialisation of gluon field
RandomSeed=2
# Given in units of az
InitialGluonSaturationScale=2.0
InitialGluonOccupationAmplitude=1.0
Coupling=0.01

FileGluonicWilsonLoop="gluonicwilsonLoop.txt"

#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $TimeRange $Rmax $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $FileGluonicWilsonLoop

