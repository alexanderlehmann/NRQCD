#!/bin/bash
## Project name as seen in the queue
##SBATCH --account=nn9578k
#SBATCH --job-name=CorrEq32
#SBATCH --array 1-20
##SBATCH --qos=devel

## Allocating amount of resources:
#SBATCH --nodes=1
## Number of tasks (aka processes) to start on each node: Pure mpi, one task per core
#SBATCH --ntasks-per-node=4

## Application running time. The wall clock time helps the scheduler assess priorities and tasks.
#SBATCH --time=0-08:00:00

## Memory usage per core
##SBATCH --mem-per-cpu=1024

#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

# Location of the executable
EXEDIR="/home/ansatt/lehmann/NRQCD/src"
EXENAME="run.exe"

EXECUTABLE=$EXEDIR/$EXENAME #"/cluster/home/propagator/software/nrqcd/run.exe"

# Parameters
mode=33

Nx=32
Ny=32
Nz=32

ax=1.0
ay=1.0
az=1.0

at=0.1
tmax=100.

# Initialisation of gluon field
RandomSeed=$SLURM_ARRAY_TASK_ID
Beta=16.
nefieldinit=30
nequilibriumtimesteps=300
MeasureEnergy=1
FileName_Energy="thermalization_energy.txt"

# Quarks
iterativeTol=1.E-8
Mass=1.570796
c0_re=1. ; c0_im=0.
c1_re=1. ; c1_im=0.
c2_re=1. ; c2_im=0.
c3_re=1. ; c3_im=0.

# Initialisation of gluon field
FileName_MesonCorrelators="mesoncorrelators_"$SLURM_ARRAY_TASK_ID".txt"
FileName_Norm="norm_"$SLURM_ARRAY_TASK_ID".txt"
#######################################################
## Making sure output is copied back and taken care of
##-----------------------------------------------------
#savefile $FileName_WilsonLoops $FileName_HybridLoops $FileName_FreeStepHybridLoops $FileName_FreeHybridLoops $FileName_QNorm $FileName_ANorm
scratchdir="scratch"$SLURM_ARRAY_TASK_ID
mkdir $scratchdir
cd $scratchdir

cp $EXECUTABLE .

#cd $SCRATCH

## Loading Software modules
## Allways be explicit on loading modules and setting run time environment!!!
#module restore system   ## Restore loaded modules to the default

## Intel-compiler
#module load intel/2019a
## Intel-MPI
#module load iimpi/2019a
## Intel-MKL
#module load imkl/2019.1.144-iimpi-2019a

## PETSc
#module load PETSc/3.10.1-intel-2018b
# Custom PETSc (no complex cluster version available)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. #$SCRATCH   
export PETSC_DIR=/home/ansatt/lehmann/software/petsc

cp $PETSC_DIR/lib/libpetsc.so.3.10 .

# List all loaded modules
#module list  

## For IntelMPI (intel toolchain), mpirun is recommended:
export I_MPI_HYDRA_ENV=$LD_LIBRARY_PATH
mpirun ./$EXENAME \
		$mode $Nx $Ny $Nz $at $ax $ay $az $tmax $RandomSeed $Beta $nefieldinit $nequilibriumtimesteps $iterativeTol $Mass $c0_re $c0_im $c1_re $c1_im $c2_re $c2_im $c3_re $c3_im $FileName_MesonCorrelators $FileName_Norm $MeasureEnergy $FileName_Energy

cp *.txt ..

cd ..
rm -rf $scratchdir


########################################################
# That was about all this time; lets call it a day...
##-------------------------------------------------------
# Finish the script
exit 0

