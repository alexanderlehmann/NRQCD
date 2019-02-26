program simulation
  use, intrinsic :: iso_fortran_env
  use precision
  use mpiinterface, only: ThisProc, NumProcs, MPIstop, intmpi, syncall
  use lattice, only: ndim
  use lattice
  use halocomm
  use xpfft
  !use gaugeconfiguration_su3
  use mpi
  use io
  implicit none

  ! Simulation parameters
  integer(int64) :: LatticeExtensions(ndim)!(ndim)
  real(fp)       :: LatticeSpacings(0:ndim)!(0:ndim)
  integer(int64) :: TimeSteps
  integer(int64) :: RandomNumberSeed

  real(fp)   :: GluonSaturationScale !qs
  real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
  real(fp)   :: GluonCoupling
  integer(int64) :: EnsembleSize
  real(fp)   :: Quarkmass
  real(fp)   :: CoMTime
  real(fp)   :: TimeRange
  !real(fp)  :: Wilsoncoeffs(nWilsonCoeffs)

  real(fp) :: GaugefixingCoefficient
  integer(int64) :: Number_of_Measurements_of_Gluondistribution
  integer(int64) :: TimePoints_between_Measurement_of_Gluondistribution
  real(fp) :: TimeBetweenGluonMeasurements
  
  ! Physical fields
  !type(GaugeConfiguration) :: GaugeConf
  
  complex(fp), allocatable :: testcomm(:)

  ! Monitoring variables
  integer(int64) :: it, idistmeas

  integer(int64) :: MemoryIndex, LatticeIndex, i
  integer(int64), allocatable :: LocalLatticeIndices(:)
  
  ! MPI
  integer(intmpi) :: mpierr, mpistatus(mpi_status_size), proc, tag, source, dest, buffersize

  real(fp), allocatable :: correlator(:)
  ! Ensemble variables
  real(fp), allocatable :: AA_correlator_ensemble(:,:,:), EE_correlator_ensemble(:,:,:)
  real(fp), allocatable :: AA_correlator(:,:), AA_correlator_stderror(:,:)
  real(fp), allocatable :: EE_correlator(:,:), EE_correlator_stderror(:,:)
  real(fp), allocatable :: GluonDistribution(:,:), GluonDistribution_stderror(:,:)
  real(fp), allocatable :: &
       GluonDistribution_AtTimePoint(:), GluonDistribution_AtTimePoint_stderror(:)
  integer :: iensemble
  real(fp), allocatable :: momenta(:)

  
  ! Output
  integer(int8) :: FileID
  character(len=128) :: filename,time_tag
  
  ! Control observable output (Gauss, Energy)
  integer(int8) :: fileID_eg
  real(fp) :: gauss, energy
  real(fp) :: time

  complex(fp), allocatable :: data(:)
  
  call InitSimulation
  
  call EndSimulation

contains
  
  !>@brief Initialisation of the simulation
  !!@details
  !! MPI\n
  !! Lattice-module\n
  !! Random number generator\n
  !! etc.
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 15.02.2019
  !!@version 1.0
  impure subroutine InitSimulation
    use precision, only: fp
    use, intrinsic :: iso_fortran_env
    
    use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
    use lattice,            only: InitModule_Lattice            => InitModule, nDim
    use halocomm,           only: InitModule_HaloComm           => InitModule
    use random,             only: InitModule_Random             => InitModule
    use xpfft,              only: InitModule_xpFFT              => InitModule
    use tolerances,         only: InitModule_tolerances         => InitModule
    implicit none

    integer(int64) :: arg_count
    character(len=80) :: arg
    integer(int8) :: i
    
    !..--** Reading simulation parameters **--..
    arg_count = 0
    
    ! Spatial lattice parameters (extensions, spacings)
    do i=1,ndim
       arg_count = arg_count +1; call get_command_argument(arg_count,arg);
       read(arg,'(I4)') LatticeExtensions(i)
    end do

    do i=0,ndim
       arg_count = arg_count +1; call get_command_argument(arg_count,arg);
       read(arg,'(F10.13)') LatticeSpacings(i)
    end do

    ! Center of mass time T
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(F10.13)') CoMTime
    ! Center of mass time-range smax
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(F10.13)') TimeRange

    TimeSteps=ceiling(TimeRange/LatticeSpacings(0))
    
    ! Seed for random number generator
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(I4)') RandomNumberSeed
    
    ! Initial gluon distribution (box): Saturation scale
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(F10.13)') GluonSaturationScale

    ! Initial gluon distribution (box): Amplitude
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(F10.13)') GluonOccupationAmplitude

    ! Coupling (only relevant in initialisation)
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(F10.13)') GluonCoupling

    ! Ensemble size (statistical average of classical statistical simulation)
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(I6)') EnsembleSize
    
    ! Time between gluon measurements
    arg_count = arg_count +1; call get_command_argument(arg_count,arg);
    read(arg,'(F10.13)') TimeBetweenGluonMeasurements
    Number_of_Measurements_of_Gluondistribution = aint(TimeRange/TimeBetweenGluonMeasurements,int64)
    
    !..--** Module initialisations **--..
    call InitModule_MPIinterface
    call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
    call InitModule_HaloComm
    call InitModule_xpFFT
    call InitModule_Random(RandomNumberSeed + ThisProc())
    call InitModule_tolerances

    call SyncAll
  end subroutine InitSimulation

  !>@brief Ending of the simulation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 15.02.2019
  !!@version 1.0
  subroutine EndSimulation
    use mpiinterface, only: FinalizeModule_MPIinterface => FinalizeModule
    use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
    implicit none

    call FinalizeModule_xpFFT
    call FinalizeModule_MPIinterface

    STOP "Simulation completed"
  end subroutine EndSimulation


  pure subroutine StatisticsCorrelator(&
       correlator_ensemble,correlator,correlator_stderror,lbound_time)
    use precision, only: fp
    use statistics, only: GetMean, GetStdError

    implicit none
    real(fp),intent(in) :: correlator_ensemble(:,lbound_time:,:)
    real(fp),allocatable,intent(out):: correlator(:,:)
    real(fp),allocatable,intent(out):: correlator_stderror(:,:)
    integer, intent(in) :: lbound_time

    integer :: latticeindex, it

    allocate(correlator(&
         lbound(correlator_ensemble,1):ubound(correlator_ensemble,1),&
         lbound(correlator_ensemble,2):ubound(correlator_ensemble,2)))

    allocate(correlator_stderror(&
         lbound(correlator_ensemble,1):ubound(correlator_ensemble,1),&
         lbound(correlator_ensemble,2):ubound(correlator_ensemble,2)))

    forall(&
         it          =lbound(correlator_ensemble,2):ubound(correlator_ensemble,2),&
         latticeindex=lbound(correlator_ensemble,1):ubound(correlator_ensemble,1))&
         correlator(latticeindex,it) = GetMean(correlator_ensemble(latticeindex,it,:))

    if(size(correlator_ensemble,3)>1) then
       forall(&
            it          =lbound(correlator_ensemble,2):ubound(correlator_ensemble,2),&
            latticeindex=lbound(correlator_ensemble,1):ubound(correlator_ensemble,1))&
            correlator_stderror(latticeindex,it) = GetStdError(correlator_ensemble(latticeindex,it,:))
    else
       correlator_stderror = 0
    end if
  end subroutine StatisticsCorrelator
end program simulation
