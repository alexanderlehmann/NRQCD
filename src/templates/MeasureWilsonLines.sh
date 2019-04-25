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

# Initialisation of gluon field
RandomSeed=2
# Given in units of az
InitialGluonSaturationScale=2.0
InitialGluonOccupationAmplitude=1.0
Coupling=0.01
# Quarks
iterativeTol=1.E-7
Mass=100.
c0_re=1. ; c0_im=0.
c1_re=1. ; c1_im=0.
c2_re=1. ; c2_im=0.
c3_re=1. ; c3_im=0.
FileGluonicWilsonLoop="gluonicwilsonLoop.txt"
FileFermionicWilsonLoop="fermionicwilsonLoop.txt"
FileMesonCorrelator="mesoncorrelator_3s1.txt"
FileNorm="norm.txt"


#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $TimeRange $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $iterativeTol $Mass $c0_re $c0_im $c1_re $c1_im $c2_re $c2_im $c3_re $c3_im $FileGluonicWilsonLoop $FileFermionicWilsonLoop $FileMesonCorrelator

