!----------------------------------------------------------------------
! Lattice module
!----------------------------------------------------------------------
!
! MODULE: lattice
!> @brief
!! Lattice index, extensions, momenta etc.
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 17.02.2019
!! @version
!! 1.0
! REVISION HISTORY:
! 17 02 2019 - Initial version
!----------------------------------------------------------------------
module lattice
  use, intrinsic :: iso_fortran_env
  implicit none

  PRIVATE
  
  public :: &
       InitModule

  !> number of dimensions
  integer(int8), parameter, public :: nDim = 3_int8

  !> number of halo points
  integer(int8), parameter, public :: nHalo= 1_int8

  !> lattice extensions
  integer(int64) :: GlobalLatticeExtensions(ndim) = -1
  !> lattice size
  integer(int64) :: GlobalLatticeSize = -1
  !> lattice spacings
  real(real64)   :: LatticeSpacings(0:ndim) = -1._real64
  !> volume
  real(real64)   :: GlobalVolume = -1._real64

  !> global maximum 2-norm of momentum
  real(real64) :: maxmomentum=-1._real64

  !..--** MPI: process-dependend parameters**--..
  !> local lattice sizes
  integer(int64), allocatable :: LocalLatticeSizes(:)
  !> local lower lattice boundaries
  integer(int64), allocatable :: LocalLowerLatticeBoundaries(:,:)
  !> local upper lattice boundaries
  integer(int64), allocatable :: LocalUpperLatticeBoundaries(:,:)
  !> local lattice extensions
  integer(int64), allocatable :: LocalLatticeExtensions(:,:)
  !> local lattice indices
  integer(int64), allocatable :: LocalLatticeIndices(:)
  !> local lattice size including halo
  integer(int64), allocatable :: LocalLatticeSizes_IncludingHalo(:)
  !> local lower lattice boundaries including halo
  integer(int64), allocatable :: LocalLowerLatticeBoundaries_IncludingHalo(:,:)
  !> local upper lattice boundaries including halo
  integer(int64), allocatable :: LocalUpperLatticeBoundaries_IncludingHalo(:,:)
  !> local lattice extensions including halo
  integer(int64), allocatable :: LocalLatticeExtensions_includingHalo(:,:)
  !> local lattice indices including halo
  integer(int64), allocatable :: LocalLatticeIndices_includingHalo(:)

  !> @brief
  !! Lattice momentum
  !! @returns
  !! Lattice momentum
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 17.02.2019
  !! @version
  !! 1.0
  !interface GetMomentum
  !   module procedure GetMomentum_BackwardDerivative
     !module procedure GetMomentum_CentralDerivative
     !module procedure GetMomentum_ForwardDerivative
  !end interface GetMomentum

contains

  impure subroutine InitModule(GlobalLatticeExtensions_,LatticeSpacings_)
    use mpi
    use mpiinterface, only: ThisProc, NumProcs, MPISTOP
    implicit none
    !> global lattice extensions
    integer(int64), intent(in) :: GlobalLatticeExtensions_(ndim)
    !> lattice spacings
    real(real64),   intent(in) :: LatticeSpacings_(0:ndim)

    integer(int64) :: partitions(ndim)
    integer :: proc
    integer(int8) :: log2_latticeext(ndim), idivision, numdivisions, ipartition
    integer(int64) :: globallatticeindex

    LatticeSpacings = LatticeSpacings_
    GlobalLatticeExtensions = GlobalLatticeExtensions_
    
    GlobalLatticeSize = product(GlobalLatticeExtensions)
    GlobalVolume = product(GlobalLatticeExtensions*LatticeSpacings(1:ndim))

    ! ..--** START: Distributing lattice points over all the MPI-processes (partitioning) **--..
    ! Dividing the lattice in sub-domains, distributed over all processes
    ! Restrictions:
    ! 1. Lattice size is integer-divisible by number of processes
    ! 2. Number of processes is a power of 2
    if(modulo(GlobalLatticeSize,NumProcs()) /= 0) then
       call MPIstop('Total lattice size has to be integer-divisible by number of MPI-processes')
    end if
    
    numdivisions = int(log(real(NumProcs(),real64))/log(2._real64),int8) !log2 of numprocs
    if(2**numdivisions /= NumProcs()) then
       call MPIstop('Number of MPI-processes has to be a power of 2')
    end if

    partitions = GlobalLatticeExtensions
    do idivision=1,numdivisions
       ! Find biggest extensions ...
       ipartition = MaxLoc(partitions, DIM=1, &
            ! Prefer highest possible dimension (fastest in memory access later on)
            BACK=.true.,&
            ! Number of points has to be divisible into 2 parts (...obviously)
            MASK = modulo(partitions,2)==0&
            )
       
       ! ... and divide it by 2
       partitions(ipartition) = partitions(ipartition)/2
    end do

    ! Check if lattice is too small for given number of processes and number of halo points
    if(any(partitions < nHalo)) &
         call MPISTOP('Lattice extensions after partitioning smaller than number of halo points')

    
  end subroutine InitModule
end module lattice
