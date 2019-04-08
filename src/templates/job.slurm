#!/bin/bash
## Project name as seen in the queue
#SBATCH --account=nn9578k
#SBATCH --job-name=Test
#SBATCH --qos=devel

## Allocating amount of resources:
#SBATCH --nodes=1
## Number of tasks (aka processes) to start on each node: Pure mpi, one task per core
#SBATCH --ntasks-per-node=4

## Application running time. The wall clock time helps the scheduler assess priorities and tasks.
#SBATCH --time=0-00:10:00

## Memory usage per core
##SBATCH --mem-per-cpu=256

#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

# Location of the executable
EXEDIR="/cluster/home/propagator/software/nrqcd"
EXENAME="run.exe"

EXECUTABLE=$EXEDIR/$EXENAME #"/cluster/home/propagator/software/nrqcd/run.exe"

# Parameters
mode=5
Nx=8
Ny=8
Nz=8
ax=1.0
ay=1.0
az=1.0
CoMTime=0.
TimeRange=0.3
TimeSpacing=0.1
# Initialisation of gluon field
RandomSeed=1
# Given in units of az
InitialGluonSaturationScale=0.7
InitialGluonOccupationAmplitude=0.2
Coupling=0.01
# Quarks
iterativeTol=1.E-7
Mass=10.
c0_re=1. ; c0_im=0.
c1_re=1. ; c1_im=0.
c2_re=1. ; c2_im=0.
c3_re=1. ; c3_im=0.
#Output-filenames
FileGluonicWilsonLoop="GluonicWilsonLoop.txt"
FileFermionicWilsonLoop="FermionicWilsonLoop.txt"

#######################################################
## Making sure output is copied back and taken care of
##-----------------------------------------------------
savefile $FileGluonicWilsonLoop $FileFermionicWilsonLoop

cp $EXECUTABLE $SCRATCH/
cd $SCRATCH

## Loading Software modules
## Allways be explicit on loading modules and setting run time environment!!!
module restore system   ## Restore loaded modules to the default

## Intel-compiler
module load intel/2018b
## Intel-MPI
module load iimpi/2018b
## Intel-MKL
module load imkl/2018.3.222-iimpi-2018b
## PETSc
module load PETSc/3.10.1-intel-2018b

# List all loaded modules
module list  

## For IntelMPI (intel toolchain), mpirun is recommended:
##mpirun ./$EXENAME \
srun ./$EXENAME \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $CoMTime $TimeRange $RandomSeed $InitialGluonSaturationScale $InitialGluonOccupationAmplitude $Coupling $iterativeTol $Mass $c0_re $c0_im $c1_re $c1_im $c2_re $c2_im $c3_re $c3_im $FileGluonicWilsonLoop $FileFermionicWilsonLoop

#########################################################
# That was about all this time; lets call it a day...
##-------------------------------------------------------
# Finish the script
exit 0

