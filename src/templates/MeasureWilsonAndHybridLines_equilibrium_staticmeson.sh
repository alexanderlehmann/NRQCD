#!/bin/bash

processes=4
exename=./run.exe

mode=52
Nx=8
Ny=8
Nz=8
ax=1.
ay=1.
az=1.

TimeRange=30.
TimeSpacing=0.1


MessDir=3
R=4

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
Beta=16.
tolerance_Eprojection=1.E-6
nefieldinit=10
nequilibriumtimesteps=1000
MeasureEnergy=1
FileName_Energy="thermalization_energy.txt"
FileName_ChargeDensity="chargedensity.txt"
FileName_WilsonLoops="wilsonloop_r"$R".txt"


#---- Creating Charge Density ----
xPositionOfPositiveCharge=1
yPositionOfPositiveCharge=1
zPositionOfPositiveCharge=1

xPositionOfNegativeCharge=$xPositionOfPositiveCharge
yPositionOfNegativeCharge=$yPositionOfPositiveCharge
zPositionOfNegativeCharge=$(($zPositionOfPositiveCharge + $R))

echo $xPositionOfPositiveCharge,$yPositionOfPositiveCharge,$zPositionOfPositiveCharge,'+0.000000E+00,+1.000000E+00,+0.000000E+00' > $FileName_ChargeDensity
echo $xPositionOfNegativeCharge,$yPositionOfNegativeCharge,$zPositionOfNegativeCharge,'+0.000000E+00,-1.000000E+00,+0.000000E+00' >> $FileName_ChargeDensity

#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $TimeRange $MessDir $R $RandomSeed $Beta $tolerance_Eprojection $nefieldinit $nequilibriumtimesteps $MeasureEnergy $FileName_Energy $FileName_ChargeDensity $FileName_WilsonLoops
