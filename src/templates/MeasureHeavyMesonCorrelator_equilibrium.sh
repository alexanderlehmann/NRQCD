#!/bin/bash

processes=4
exename=./run.exe

mode=33
Nx=12
Ny=12
Nz=12
ax=1.0
ay=1.0
az=1.0

at=0.1
tmax=10.


# Initialisation of gluon field
RandomSeed=1
Beta=16.
nefieldinit=30
nequilibriumtimesteps=300

# Quarks
iterativeTol=1.E-8
Mass=1.570796
c0_re=1. ; c0_im=0.
c1_re=1. ; c1_im=0.
c2_re=1. ; c2_im=0.
c3_re=1. ; c3_im=0.

FileMesonCorrelators="mesoncorrelators.txt"
FileNorm="norm.txt"

#--------------------------------------------------------------

mpiexec -n $processes \
	$exename \
	$mode $Nx $Ny $Nz $at $ax $ay $az $tmax $RandomSeed $Beta $nefieldinit $nequilibriumtimesteps $iterativeTol $Mass $c0_re $c0_im $c1_re $c1_im $c2_re $c2_im $c3_re $c3_im $FileMesonCorrelators $FileNorm
