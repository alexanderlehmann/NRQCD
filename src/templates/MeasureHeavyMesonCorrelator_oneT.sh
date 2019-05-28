#!/bin/bash

processes=4
exename=./run.exe

mode=31
Nx=16
Ny=16
Nz=16
ax=1.0
ay=1.0
az=1.0

at=0.1
t=100.
smax=2.
UpdateQuarksEveryNsteps=5

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
InitialGluonSaturationScale=1.
InitialGluonOccupationAmplitude=1.
Coupling=0.001
# Quarks
iterativeTol=1.E-8
Mass=1.570796
c0_re=1. ; c0_im=0.
c1_re=1. ; c1_im=0.
c2_re=1. ; c2_im=0.
c3_re=1. ; c3_im=0.

FileMesonCorrelators="mesoncorrelators_gherm_5.txt"
FileNorm="norm.txt"
FileMesonCorrelators_gconjg="mesoncorrelators_gconjg_5.txt"

#--------------------------------------------------------------

mpiexec -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $at $ax $ay $az $t $smax $UpdateQuarksEveryNsteps $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $iterativeTol $Mass $c0_re $c0_im $c1_re $c1_im $c2_re $c2_im $c3_re $c3_im $FileMesonCorrelators $FileNorm $FileMesonCorrelators_gconjg
