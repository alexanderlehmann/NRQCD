#!/bin/bash

processes=4
exename=./run.exe

Nx=8
Ny=8
Nz=8
ax=1.0
ay=1.0
az=1.0

CoMTime=0.
TimeRange=10.
TimeSpacing=0.1

Ensemblesize=100

# Initialisation of gluon field
RandomSeed=13
# Given in units of az
InitialGluonSaturationScale=2.0
InitialGluonOccupationAmplitude=1.0
Coupling=0.0
TimeBetweenGluonMeasurement=10.
GaugefixingMaxSteps=100000
GaugefixingTolerance=1.E-8
GaugefixingCoefficient=0.08
# Quarks
Mass=5.


mpiexec -n $processes \
	$exename \
	$Nx $Ny $Nz $TimeSpacing $ax $ay $az $CoMTime $TimeRange $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $Ensemblesize $TimeBetweenGluonMeasurement $GaugefixingMaxSteps $GaugefixingTolerance $GaugefixingCoefficient
