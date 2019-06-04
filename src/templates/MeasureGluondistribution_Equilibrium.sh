#!/bin/bash

processes=4
exename=./run.exe

mode=11
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

TimeRange=0.
TimeSpacing=0.1
TimeBetweenGluonMeasurement=1.

FileNameGluondist="gluondist.txt"

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
Beta=16.
nefieldinit=20
nequilibriumtimesteps=200

GaugefixingMaxSteps=100000
GaugefixingTolerance=1.E-10
GaugefixingCoefficient=0.08
MeasureThermalizationEnergy=1
FileName_ThermalizationEnergy="thermalization_energy.txt"

mpiexec -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $TimeRange $RandomSeed $Beta $nefieldinit $nequilibriumtimesteps $TimeBetweenGluonMeasurement $GaugefixingMaxSteps $GaugefixingTolerance $GaugefixingCoefficient $FileNameGluondist $MeasureThermalizationEnergy $FileName_ThermalizationEnergy
