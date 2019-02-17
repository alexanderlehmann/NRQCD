program simulation
  use, intrinsic :: iso_fortran_env

  use mpiinterface, only: ThisProc, NumProcs, MPIstop
  use lattice, only: ndim
  implicit none

  ! Simulation parameters
  integer(int64) :: LatticeExtensions(ndim)!(ndim)
  real(real64)   :: LatticeSpacings(0:ndim)!(0:ndim)
  integer(int64) :: TimeSteps
  integer(int64) :: RandomNumberSeed

  real(real64)   :: GluonSaturationScale !qs
  real(real64)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
  real(real64)   :: GluonCoupling
  integer(int64) :: EnsembleSize
  real(real64)   :: Quarkmass
  real(real64)   :: CoMTime
  real(real64)   :: TimeRange
  !real(real64)  :: Wilsoncoeffs(nWilsonCoeffs)

  call InitSimulation

  call mpistop
  
  call EndSimulation

contains
  
  !> @brief
  !! Initialisation of the simulation
  !! @details
  !! MPI\n
  !! Lattice-module\n
  !! Random number generator\n
  !! etc.
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 15.02.2019
  !! @version 1.0
  impure subroutine InitSimulation
    use, intrinsic :: iso_fortran_env
    use mpi
    
    use mpiinterface,       only: InitModule_MPIinterface       => InitModule, thisProc
    use lattice,            only: InitModule_Lattice            => InitModule, ndim
    use halocomm,           only: InitModule_HaloComm           => InitModule
    use random,             only: InitModule_Random             => InitModule
    implicit none

    integer(int64) :: arg_count
    character(len=80) :: arg
    integer(int8) :: i

    integer(int32) :: mpierr
    
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
    
    !..--** Module initialisations **--..
    call InitModule_MPIinterface
    call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
    call InitModule_HaloComm
    call InitModule_Random(RandomNumberSeed + ThisProc())

    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end subroutine InitSimulation

  subroutine EndSimulation
    use mpiinterface, only: FinalizeModule_MPIinterface => FinalizeModule
    !use dft,          only: FinalizeModule_DFT          => FinalizeModule
    implicit none

    !call FinalizeModule_DFT
    call FinalizeModule_MPIinterface

    STOP
  end subroutine EndSimulation
end program simulation
