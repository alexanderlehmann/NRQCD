#!/bin/bash

processes=4
exename=./run.exe

mode=1
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

CoMTime=0.
TimeRange=100.
TimeSpacing=0.1

Ensemblesize=10

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
InitialGluonSaturationScale=1.
InitialGluonOccupationAmplitude=1.
Coupling=0.0
TimeBetweenGluonMeasurement=1.
GaugefixingMaxSteps=10000
GaugefixingTolerance=1.E-8
GaugefixingCoefficient=0.08

mpiexec -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $CoMTime $TimeRange $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $Ensemblesize $TimeBetweenGluonMeasurement $GaugefixingMaxSteps $GaugefixingTolerance $GaugefixingCoefficient