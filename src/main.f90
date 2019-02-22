program simulation
  use, intrinsic :: iso_fortran_env
  use precision
  use mpiinterface, only: ThisProc, NumProcs, MPIstop, intmpi
  use lattice, only: ndim

  use lattice, only: GetLocalLatticeSize_IncludingHalo, GetLocalLatticeSize
  use lattice
  use halocomm
  use xpfft
  use gaugeconfiguration_su3
  implicit none

  ! Simulation parameters
  integer(int64) :: LatticeExtensions(ndim)!(ndim)
  real(fp)   :: LatticeSpacings(0:ndim)!(0:ndim)
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

  complex(fp), allocatable :: testcomm(:)


  integer(int64) :: LocalIndex, LatticeIndex, i
  integer(intmpi) :: proc
  integer(int64), allocatable :: localindices(:)
  complex(fp), allocatable :: data(:)
  
  call InitSimulation
  
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
    
    use mpiinterface,       only: InitModule_MPIinterface       => InitModule, thisProc, intmpi
    use lattice,            only: InitModule_Lattice            => InitModule, ndim
    use halocomm,           only: InitModule_HaloComm           => InitModule
    use random,             only: InitModule_Random             => InitModule
    use xpfft,              only: InitModule_xpFFT              => InitModule
    implicit none

    integer(int64) :: arg_count
    character(len=80) :: arg
    integer(int8) :: i

    integer(intmpi) :: mpierr
    
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
    call InitModule_xpFFT
    call InitModule_Random(RandomNumberSeed + ThisProc())

    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end subroutine InitSimulation

  subroutine EndSimulation
    use mpiinterface, only: FinalizeModule_MPIinterface => FinalizeModule
    use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
    implicit none

    call FinalizeModule_xpFFT
    call FinalizeModule_MPIinterface

    STOP "Simulation completed"
  end subroutine EndSimulation
end program simulation
