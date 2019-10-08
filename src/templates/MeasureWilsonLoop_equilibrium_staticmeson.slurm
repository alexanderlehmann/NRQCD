#!/bin/bash
## Project name as seen in the queue
##SBATCH --account=nn9578k
#SBATCH --job-name=sWL8_4
#SBATCH --array 1-2
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

EXECUTABLE=$EXEDIR/$EXENAME

# Parameters
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
RandomSeed=$SLURM_ARRAY_TASK_ID
# Given in units of az
Beta=16.
tolerance_Eprojection=1.E-6
nefieldinit=10
nequilibriumtimesteps=1000
MeasureEnergy=1
FileName_Energy="thermalization_energy_"$SLURM_ARRAY_TASK_ID".txt"
FileName_ChargeDensity="chargedensity.txt"
FileName_WilsonLoops="wilsonloop_r"$R"_"$SLURM_ARRAY_TASK_ID".txt"

#######################################################
## 1. Creating scratch directory
##-----------------------------------------------------
scratchdir="scratch"$SLURM_ARRAY_TASK_ID
mkdir $scratchdir
cd $scratchdir
cp $EXECUTABLE .

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. #$SCRATCH   
export PETSC_DIR=/home/ansatt/lehmann/software/petsc

export I_MPI_HYDRA_ENV=$LD_LIBRARY_PATH

#---- Creating Charge Density ----
xPositionOfPositiveCharge=1
yPositionOfPositiveCharge=1
zPositionOfPositiveCharge=1

xPositionOfNegativeCharge=$xPositionOfPositiveCharge
yPositionOfNegativeCharge=$yPositionOfPositiveCharge
zPositionOfNegativeCharge=$(($zPositionOfPositiveCharge + $R))

echo $xPositionOfPositiveCharge,$yPositionOfPositiveCharge,$zPositionOfPositiveCharge,'+0.000000E+00,+1.000000E+00,+0.000000E+00' > $FileName_ChargeDensity
echo $xPositionOfNegativeCharge,$yPositionOfNegativeCharge,$zPositionOfNegativeCharge,'+0.000000E+00,-1.000000E+00,+0.000000E+00' >> $FileName_ChargeDensity


mpirun  ./$EXENAME \
	$mode $Nx $Ny $Nz $TimeSpacing $ax $ay $az $TimeRange $MessDir $R $RandomSeed $Beta $tolerance_Eprojection $nefieldinit $nequilibriumtimesteps $MeasureEnergy $FileName_Energy $FileName_ChargeDensity $FileName_WilsonLoops

#######################################################
## 2. Making sure output is copied back and taken care of
##-----------------------------------------------------
cp *.txt ..

cd ..
rm -rf $scratchdir


########################################################
# That was about all this time; lets call it a day...
##-------------------------------------------------------
# Finish the script
exit 0
