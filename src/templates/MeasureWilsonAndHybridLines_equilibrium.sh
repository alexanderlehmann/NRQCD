#!/bin/bash

processes=4
exename=./run.exe

mode=51
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

StartTime=0.
TimeRange=30.
TimeSpacing=0.01
Rmax=6

# Initialisation of gluon field
RandomSeed=1
# Given in units of az
Beta=16.
nefieldinit=10
nequilibriumtimesteps=300

FileName_WilsonLoops="wilsonloops.txt"

# Quarks
MeasureHybridLoops=0
iterativeTol=1.E-8
Mass=100.
c0_re=1. ; c0_im=0.
c1_re=1. ; c1_im=0.
c2_re=1. ; c2_im=0.
c3_re=1. ; c3_im=0.




FileName_HybridLoops="hybridloops.txt"
FileName_FreeStepHybridLoops="freestephybridloops.txt"
FileName_FreeHybridLoops="freehybridloops.txt"
FileName_QNorm="quarknorm.txt"
FileName_ANorm="antiquarknorm.txt"

#--------------------------------------------------------------

mpirun -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $StartTime $TimeRange $Rmax $RandomSeed $Beta $nefieldinit $nequilibriumtimesteps $MeasureHybridLoops $iterativeTol $Mass $c0_re $c0_im $c1_re $c1_im $c2_re $c2_im $c3_re $c3_im $FileName_WilsonLoops $FileName_HybridLoops $FileName_FreeStepHybridLoops $FileName_FreeHybridLoops $FileName_QNorm $FileName_ANorm
