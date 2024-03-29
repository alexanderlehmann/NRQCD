!----------------------------------------------------------------------
! Lattice module
!----------------------------------------------------------------------
!
! MODULE: lattice
!>@brief Lattice index, extensions, momenta etc.
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 17.02.2019
!!@version 1.0
! REVISION HISTORY:
! 17 02 2019 - Initial version
!----------------------------------------------------------------------
module lattice
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  implicit none

  PRIVATE

  public :: &
       InitModule, &
       IsModuleInitialised,&
       !GetLocalLatticeIndices,&
       !GetLocalLatticeIndices_allocatable,&
       GetLatticePosition,&
       GetNegativeLatticeIndex,&
       GetProc_G, GetProc_M,&
       GetLatticeIndex, GetLatticeIndex_M, GetMemoryIndex,&
       GetNeib_M, GetNeib_G,&
       GetMomentum_G, GetMomentum_M,&
       GetNorm2Momentum_G, GetNorm2Momentum_M,&
       GetProc_fromGeneralIndex,&
       GetLatticeExtension,&
       GetLatticeSize,&
       InitLatticeIndices,&
       GetMemorySize,&
       GetLocalLatticeSize,&
       GetLocalLowerLatticeBoundary,&
       GetLocalUpperLatticeBoundary,&
       GetIndex_fromPosition,&
       GetLatticeSpacing,&
       GetVolume,&
       GetMaxNorm2Momentum,&
       GetPolarisationVectors,&
       GetTransverseProjector

  !> Module name
  character(len=7), parameter, public ::  modulename='lattice'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.
  
  !> number of dimensions
  integer(int8), parameter, public :: nDim = 3_int8

  !> number of halo points
  integer(int8), parameter, public :: nHalo= 2_int8

  !> lattice extensions
  integer(int64) :: LatticeExtensions(ndim) = -1
  !> lattice size
  integer(int64) :: LatticeSize = -1
  !> local lattice size
  integer(int64) :: LocalLatticeSize = -1
  !> lattice spacings
  real(fp)   :: LatticeSpacings(0:ndim) = -1._fp
  !> volume
  real(fp)   :: Volume = -1._fp
  !> local lattice size including halo
  integer(int64) :: MemorySize=-1

  !> maximum 2-norm of momentum
  real(fp) :: MaxNorm2Momentum=-1._fp

  !..--** MPI: process-dependend parameters**--..
  !> local lower lattice boundaries
  integer(int64), allocatable :: LocalLowerLatticeBoundaries(:,:)
  !> local upper lattice boundaries
  integer(int64), allocatable :: LocalUpperLatticeBoundaries(:,:)
  !> local lower lattice boundaries including halo
  integer(int64), allocatable :: LocalLowerLatticeBoundaries_IncludingHalo(:,:)
  !> local upper lattice boundaries including halo
  integer(int64), allocatable :: LocalUpperLatticeBoundaries_IncludingHalo(:,:)
  !> local lattice indices
  integer(int64), allocatable :: LocalLatticeIndices(:)
  !> memory indices of neighbouring lattice sites
  integer(int64), allocatable :: neib_m(:,:)
  !> local lattice extensions
  integer(int64) :: LocalLatticeExtensions(ndim)
  !> local lattice extensions including halo
  integer(int64) :: LocalLatticeExtensions_includingHalo(ndim)

  !>@brief Lattice momentum
  !!@returns Lattice momentum
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  interface GetMomentum_G
     module procedure GetMomentum_BackwardDerivative
     !module procedure GetMomentum_CentralDerivative
     !module procedure GetMomentum_ForwardDerivative
  end interface GetMomentum_G
contains
  
  !>@brief Returns, if module is initialised
  !!@returns module's initialisation status
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure logical function IsModuleInitialised()
    implicit none
    IsModuleInitialised = IsInitialised
  end function IsModuleInitialised
  
  !>@brief Initialises module
  !!@details Lattice-extensions, lattice-size, volume etc
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  impure subroutine InitModule(LatticeExtensions_,LatticeSpacings_)
    use mpiinterface, only: ThisProc, NumProcs, MPISTOP, intmpi
    use arrayoperations, only: RemoveDuplicates, Sort
    use, intrinsic :: iso_fortran_env
    implicit none
    !> lattice extensions
    integer(int64), intent(in) :: LatticeExtensions_(ndim)
    !> lattice spacings
    real(fp),   intent(in) :: LatticeSpacings_(0:ndim)

    integer(int64) :: partitions(ndim)
    integer(int8) :: idivision, numdivisions, ipartition
    integer(int64) :: latticeindex

    integer(intmpi) :: proc

    if(isInitialised) then
       call MPISTOP('Error in init of '//modulename//': already initialised.')
    else
       call CheckObligatoryInitialisations

       LatticeSpacings = LatticeSpacings_
       LatticeExtensions = LatticeExtensions_

       LatticeSize = product(LatticeExtensions)
       Volume = product(LatticeExtensions*LatticeSpacings(1:ndim))

       ! ..--** START: Distributing lattice points over all the MPI-processes (partitioning) **--..
       ! Dividing the lattice in sub-domains, distributed over all processes
       ! Restrictions:
       ! 1. Lattice size is integer-divisible by number of processes
       ! 2. Number of processes is a power of 2
       if(modulo(LatticeSize,NumProcs()) /= 0) then
          call MPIstop('Total lattice size has to be integer-divisible by number of MPI-processes')
       end if

       numdivisions = int(log(real(NumProcs(),real64))/log(2._real64),int8) !log2 of numprocs

       if(2**numdivisions /= NumProcs()) then
          call MPIstop('Number of MPI-processes has to be a power of 2')
       end if

       partitions = LatticeExtensions
       do idivision=1,numdivisions
          ! Find biggest extensions ...
          ! Prefer highest possible dimension (fastest in memory access later on)
          ! Number of points has to be divisible into 2 parts (...obviously)
          ipartition = MaxLoc(partitions, DIM=1, &
               BACK=.true.,&                    ! Prefering highest possible dimension
               MASK = modulo(partitions,2)==0&  ! Number of points divisible by 2
               )

          ! ... and divide it by 2
          partitions(ipartition) = partitions(ipartition)/2
       end do

       ! Check if lattice is too small for given number of processes and number of halo points
       if(any(partitions < nHalo)) &
            call MPISTOP('Lattice extensions after partitioning smaller than number of halo points')

       ! Local lattice boundaries (including and without halo)
       call InitLocalLatticeBoundaries(&
            LocalLowerLatticeBoundaries,&
            LocalUpperLatticeBoundaries,&
            partitions,nHalo=0_int8)
       call InitLocalLatticeBoundaries(&
            LocalLowerLatticeBoundaries_includingHalo,&
            LocalUpperLatticeBoundaries_includingHalo,&
            partitions,nHalo=nHalo)

       ! Local lattice points
       call InitLatticeIndices(&
            LocalLatticeIndices,&
            LocalLowerLatticeBoundaries_includingHalo(:,ThisProc()+1),&
            LocalUpperLatticeBoundaries_includingHalo(:,ThisProc()+1))

       MemorySize = size(LocalLatticeIndices)

       ! Deleting lattice points which are contained more than once ...
       ! (important in case of small partitions)
       call RemoveDuplicates(&
            LocalLatticeIndices)
       ! ... and sort the arrays
       call Sort(LocalLatticeIndices)

       ! Setting local lattice size (including and without halo)
       MemorySize = size(LocalLatticeIndices)

       LocalLatticeExtensions= 1 + LocalUpperLatticeBoundaries(:,1) -LocalLowerLatticeBoundaries(:,1)
       LocalLatticeExtensions_includingHalo=&
            +1&
            +LocalUpperLatticeBoundaries_includingHalo(:,1)&
            -LocalLowerLatticeBoundaries_includingHalo(:,1)
       
       LocalLatticeSize = product(LocalLatticeExtensions)
       
       call InitNeib_m(Neib_m)
       
       ! Find maximum momentum
       MaxNorm2Momentum = FindMaxNorm2Momentum()

       ! DONE
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule

  !>@brief Initialising the local neighbour indices including the ghost region
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  pure subroutine InitNeib_m(neib_m)
    use mpiinterface, only: ThisProc
    implicit none
    !> Neighbouring index field for memory indices
    integer(int64), allocatable, intent(out) :: neib_m(:,:)
    
    integer(int8)  :: i
    integer(int64) :: MemoryIndex

    allocate(neib_m(-nDim:+nDim,MemorySize))
    neib_m = huge(MemoryIndex)
    forall(MemoryIndex=1:MemorySize, i=-nDim:+nDim&
         !,any(LocalLatticeIndices==GetNeib_G(i,GetLatticeIndex_M(MemoryIndex)))
         )
       neib_m(i,MemoryIndex) = GetMemoryIndex(GetNeib_G(i,GetLatticeIndex_M(MemoryIndex)))
    end forall
  end subroutine InitNeib_M

  !>@brief Neighbour index routine
  !!@returns Memory index of lattice neighbour in i'th direction
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  pure elemental integer(int64) function GetNeib_M(i,MemoryIndex)
    implicit none
    !> Direction
    integer(int8),  intent(in) :: i
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    GetNeib_M = Neib_M(i,MemoryIndex)
  end function GetNeib_M
  
  !>@brief Checks previous necessary initialisations
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
  impure subroutine CheckObligatoryInitialisations
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: IsMPIinterfaceInitialised => IsModuleInitialised, mpiname => modulename,&
         mpistop
    implicit none
    
    if(.not.IsMPIinterfaceInitialised()) then
       call mpistop('Error in init of '//modulename//': '//mpiname//' is not initialised.')
    end if
  end subroutine CheckObligatoryInitialisations

  !>@brief Initialises local lattice boundaries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure recursive subroutine InitLocalLatticeBoundaries(&
       LocalLowerLatticeBoundaries,LocalUpperLatticeBoundaries,partitions,&
       proc,i,ipart,nHalo)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: NumProcs
    implicit none
    !> Lower boundaries of local partition
    integer(int64), intent(inout), allocatable :: LocalLowerLatticeBoundaries(:,:)
    !> Upper boundaries of local partition
    integer(int64), intent(inout), allocatable :: LocalUpperLatticeBoundaries(:,:)
    !> Global sizes of each partition
    integer(int64), intent(in)    :: partitions(ndim)
    !> Process number
    integer(int64), intent(inout), optional :: proc
    !> Direction
    integer(int8),  intent(in),    optional :: i
    !> Direction, used as index in partitions-array
    integer(int64), intent(in),    optional :: ipart(ndim)
    !> Number of halo-points (boundary values)
    integer(int8),  intent(in),    optional :: nHalo

    
    integer(int8)  :: i_next
    integer(int64) :: ipart_next(ndim)

    integer(int8) :: j
    integer(int64) :: proc_init
    integer(int64) :: ip
    integer(int8) :: nHalo_

    if(present(nHalo)) then
       nHalo_=nHalo
    else ! Default value
       nHalo_=0_int8
    end if

    ! If i, ipart or proc not given --> start of recursion
    if(.not.present(i).or..not.present(proc).or..not.present(ipart)) then
       i_next    =ndim
       ipart_next=1
       proc_init =0

       if(allocated(LocalLowerLatticeBoundaries)) deallocate(LocalLowerLatticeBoundaries)
       if(allocated(LocalUpperLatticeBoundaries)) deallocate(LocalUpperLatticeBoundaries)
       allocate(LocalLowerLatticeBoundaries(ndim,NumProcs()))
       allocate(LocalUpperLatticeBoundaries(ndim,NumProcs()))
       call InitLocalLatticeBoundaries(LocalLowerLatticeBoundaries,LocalUpperLatticeBoundaries,&
            partitions,proc_init,i_next,ipart_next,nHalo_)
    else
       ipart_next=ipart
       if(i>0) then
          ! Recursive do-loop continued and handed over to a further do-loop
          do ip=1,LatticeExtensions(i)/partitions(i)
             ipart_next(i)=ip
             i_next       =i-1
             call InitLocalLatticeBoundaries(&
                  LocalLowerLatticeBoundaries,LocalUpperLatticeBoundaries,partitions,&
                  proc,i_next,ipart_next,nHalo_)
          end do
       else
          ! End of recursion
          do j=1,ndim
             LocalLowerLatticeBoundaries(j,proc+1)=&
                  1 + (ipart(j)-1)*partitions(j) -nHalo_
             LocalUpperLatticeBoundaries(j,proc+1)=&
                  ipart(j)*partitions(j) + nHalo_
          end do
          proc = proc+1
       end if
    end if

  end subroutine InitLocalLatticeBoundaries

  !>@brief Initialises local lattice indices
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure recursive subroutine InitLatticeIndices(LatticeIndices,&
       LowerBoundaries,UpperBoundaries,i,x,index)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> To be filled array of lattice indices
    integer(int64), intent(inout), allocatable :: LatticeIndices(:)
    !> Lower boundaries of local partition
    integer(int64), intent(in)                 :: LowerBoundaries(ndim)
    !> Upper boundaries of local partition
    integer(int64), intent(in)                 :: UpperBoundaries(ndim)
    !> Direction
    integer(int8),  intent(in),    optional    :: i
    !> Position vector
    integer(int64), intent(in),    optional    :: x(ndim)
    !> Local array index in LatticeIndices
    integer(int64), intent(inout), optional    :: index

    integer(int64) :: v(ndim), n, index_init, latticesize
    integer(int8) :: j

    ! If i or x not given --> start of recursion
    if(.not.present(i).or..not.present(x).or..not.present(index)) then
       latticesize = product(1+UpperBoundaries-LowerBoundaries)
       if(allocated(LatticeIndices)) deallocate(LatticeIndices)
       allocate(LatticeIndices(latticesize))

       index_init=1
       call InitLatticeIndices(LatticeIndices,LowerBoundaries,UpperBoundaries,&
            ndim,LowerBoundaries,index_init)
    else
       v = x
       if(i>0) then
          ! Recursive do-loop continued and handed over to a further do-loop
          do n=LowerBoundaries(i),&
               UpperBoundaries(i)
             v(i) = n

             call InitLatticeIndices(LatticeIndices,LowerBoundaries,UpperBoundaries,&
                  i-1_int8,v,index)
          end do
       else
          ! End of recursion
          v = GetPeriodicPosition(x,LatticeExtensions)
          LatticeIndices(index) = GetLatticeIndex(v)
          index = index+1
       end if
    end if
  end subroutine InitLatticeIndices

  !>@brief Finds maximum of \f$||\vec{p}||_2\f$
  !!@returns maximum of \f$||\vec{p}||_2\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  impure real(fp) function FindMaxNorm2Momentum
    use, intrinsic :: iso_fortran_env
    use mpi
    use mpiinterface, only: ThisProc, NumProcs, intmpi, MPIstop
    implicit none

    integer(int64) :: MemoryIndex, LatticeIndex
    real(fp)   :: norm2momentum,currentMax
    real(fp), allocatable :: maxmomenta(:)

    integer(intmpi) :: mpierr, sendtype
    
    currentMax = 0
    do MemoryIndex=1,MemorySize
       LatticeIndex = GetLatticeIndex_M(MemoryIndex)
       if(GetProc_G(LatticeIndex)==ThisProc()) then
          norm2momentum = GetNorm2Momentum_G(LatticeIndex)
          if(norm2momentum.gt.currentMax) currentMax = norm2Momentum
       end if
    end do

    ! Communicate all local maxima
    allocate(maxmomenta(NumProcs()))

    select case(fp)
    case(real32)
       sendtype = MPI_REAL4
    case(real64)
       sendtype = MPI_REAL8
    case(real128)
       sendtype = MPI_REAL16
    case default
       call MPISTOP('Error in FindMaxNorm2Momentum of '//modulename&
            //': unsupported floating point precision.')
    end select
    
    call mpi_allgather(&
         currentMax,1,sendtype,& ! Sending
         maxmomenta,1,sendtype,& ! Recieving
         MPI_COMM_WORLD,mpierr)

    FindMaxNorm2Momentum = maxval(maxmomenta)
    deallocate(maxmomenta)
  end function FindMaxNorm2Momentum

  !>@brief Returns lattice index
  !!@returns lattice index
  !!@details
  !! Lattice index is given by
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure integer(int64) function GetLatticeIndex(Position)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> position vector
    integer(int64), intent(in) :: Position(ndim)

    GetLatticeIndex &
         = GetIndex_fromPosition(Position,LatticeExtensions)
  end function GetLatticeIndex

  !>@brief Returns the lattice index corresponding to \f$-j\f$
  !!@returns lattice index
  !!@details
  !! Lattice index is given by
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
 pure integer(int64) function GetNegativeLatticeIndex(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex

    integer(int64) :: Position(ndim), NegativePosition(ndim)
    Position = GetLatticePosition(LatticeIndex)
    NegativePosition = GetNegativeLatticePosition(Position)
    GetNegativeLatticeIndex = GetLatticeIndex(NegativePosition)
  end function GetNegativeLatticeIndex

  !>@brief Returns the lattice position corresponding to \f$-j\f$
  !!@returns Lattice position corresponding to \f$-j\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
 pure function GetNegativeLatticePosition(Position)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Lattice position
    integer(int64), intent(in) :: Position(ndim)
    !> Negative lattice position
    integer(int64) :: GetNegativeLatticePosition(ndim)
    GetNegativeLatticePosition&
         = 1_int64 + modulo(LatticeExtensions - (Position-1_int64),LatticeExtensions)
  end function GetNegativeLatticePosition
    
  !>@brief Returns lattice index of neighbouring point in i'th direction
  !!@returns lattice index of neighbouring point in i'th direction
  !!@details Periodic shift i'th direction
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure elemental integer(int64) function GetNeib_G(i,LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: LatticeIndex

    integer(int64) :: LatticePosition(ndim)

    GetNeib_G &
         = GetNeib_FromIndex(i,LatticeIndex,LatticeExtensions)
  end function GetNeib_G

  !>@brief Returns local lattice indices without halo
  !!@returns local lattice indices without halo
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure subroutine GetLocalLatticeIndices(indices)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> local lattice indices
    integer(int64), intent(out) :: indices(:)
    indices = LocalLatticeIndices
  end subroutine GetLocalLatticeIndices

  !>@brief Returns local lattice indices without halo with allocation
  !!@returns local lattice indices without halo with allocation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure subroutine GetLocalLatticeIndices_Allocatable(indices)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> local lattice indices
    integer(int64), intent(out), allocatable :: indices(:)
    allocate(indices(MemorySize))
    indices = LocalLatticeIndices
  end subroutine GetLocalLatticeIndices_Allocatable

  !>@brief  Returns volume
  !!@returns volume
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure real(fp) function GetVolume()
    use, intrinsic :: iso_fortran_env
    implicit none
    GetVolume = Volume
  end function GetVolume
  
!>@brief  Returns lattice spacing
  !!@returns lattice spacing
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
 pure elemental real(fp) function GetLatticeSpacing(i)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Direction \f$\in\{0,...,ndim\}\f$
    integer(int8), intent(in) :: i
    GetLatticeSpacing = LatticeSpacings(i)
  end function GetLatticeSpacing
  
  !>@brief Returns i'th lattice extension
  !!@returns i'th lattice extension
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
 pure elemental integer(int64) function GetLatticeExtension(i)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), intent(in) :: i
    GetLatticeExtension = LatticeExtensions(i)
  end function GetLatticeExtension

  !>@brief Returns i'th lower lattice boundaries
  !!@returns i'th lower lattice boundaries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date  11.03.2019
  !!@version 1.0
  pure elemental integer(int64) function GetLocalLowerLatticeBoundary(i,proc)
    use mpiinterface, only: intmpi, thisproc
    implicit none
    !> Direction
    integer(int8),   intent(in) :: i
    !> MPI-process rank
    integer(intmpi), intent(in), optional :: proc

    if(present(proc)) then
       GetLocalLowerLatticeBoundary = LocalLowerLatticeBoundaries(i,proc+1_intmpi)
    else
       GetLocalLowerLatticeBoundary = LocalLowerLatticeBoundaries(i,ThisProc()+1_intmpi)
    end if
  end function GetLocalLowerLatticeBoundary

  !>@brief Returns i'th upper lattice boundaries
  !!@returns i'th upper lattice boundaries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date  11.03.2019
  !!@version 1.0
  pure elemental integer(int64) function GetLocalUpperLatticeBoundary(i,proc)
    use mpiinterface, only: intmpi, thisproc
    implicit none
    !> Direction
    integer(int8),   intent(in) :: i
    !> MPI-process rank
    integer(intmpi), intent(in), optional :: proc

    if(present(proc)) then
       GetLocalUpperLatticeBoundary = LocalUpperLatticeBoundaries(i,proc+1_intmpi)
    else
       GetLocalUpperLatticeBoundary = LocalUpperLatticeBoundaries(i,ThisProc()+1_intmpi)
    end if
  end function GetLocalUpperLatticeBoundary
  
  !>@brief Returns lattice size
  !!@returns Number of points on the whole lattice
  !!@details Lattice size is \f$|\Lambda|=\prod\limits_{i=1}^{d}N_i\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure integer(int64) function GetLatticeSize()
    use, intrinsic :: iso_fortran_env
    implicit none
    GetLatticeSize = LatticeSize
  end function GetLatticeSize

  !>@brief Returns local lattice size including halo
  !!@returns local lattice size including halo
  !!@details Lattice size is \f$|\Lambda|=\prod\limits_{i=1}^{d}N_i\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure integer(int64) function GetMemorySize()
    use, intrinsic :: iso_fortran_env
    implicit none
    GetMemorySize = MemorySize
  end function GetMemorySize

  !>@brief Returns local lattice size
  !!@returns local lattice size
  !!@details Lattice size is \f$|\Lambda|=\prod\limits_{i=1}^{d}N_i\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure integer(int64) function GetLocalLatticeSize()
    use, intrinsic :: iso_fortran_env
    implicit none
    GetLocalLatticeSize = LocalLatticeSize
  end function GetLocalLatticeSize
  
  !>@brief Returns MPI-process-rank corresponding to lattice index
  !!@returns MPI-process-rank corresponding to lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.1
 pure elemental integer(intmpi) function GetProc_G(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: intmpi
    implicit none
    !> lattice index
    integer(int64), intent(in) :: LatticeIndex
    GetProc_G = GetProc_fromGeneralIndex(LatticeIndex,&
         LocalLowerLatticeBoundaries,LocalUpperLatticeBoundaries,LatticeExtensions)
  end function GetProc_G

  !>@brief Returns MPI-process-rank corresponding to memory index
  !!@returns MPI-process-rank corresponding to memory index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  pure elemental integer(intmpi) function GetProc_M(MemoryIndex)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: intmpi
    implicit none
    !> lattice index
    integer(int64), intent(in) :: MemoryIndex
    GetProc_M = GetProc_G(GetLatticeIndex_M(MemoryIndex))
  end function GetProc_M

  !>@brief General function, returning MPI-process-rank corresponding to given index
  !!@returns MPI-process-rank corresponding to index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
 pure integer(intmpi) function GetProc_fromGeneralIndex(&
       Index,LowerBoundaries,UpperBoundaries,Extensions)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: intmpi
    implicit none
    !> index
    integer(int64), intent(in) :: Index
    !> Lower boundaries
    integer(int64), intent(in) :: LowerBoundaries(:,:)
    !> Upper boundaries
    integer(int64), intent(in) :: UpperBoundaries(:,:)
    !> Extensions
    integer(int64), intent(in) :: Extensions(ndim)
    
    integer(int64) :: Position(ndim)
    integer(intmpi) :: proc
    
    Position = GetPosition_fromIndex(Index,Extensions)

    GetProc_fromGeneralIndex = GetProc_fromGeneralPosition(Position,&
         LowerBoundaries,UpperBoundaries)
  end function GetProc_fromGeneralIndex

  !>@brief General function, returning MPI-process-rank corresponding to given position
  !!@returns MPI-process-rank corresponding to index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
 pure integer(intmpi) function GetProc_fromGeneralPosition(&
       Position,LowerBoundaries,UpperBoundaries)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: intmpi
    implicit none
    !> Position
    integer(int64), intent(in) :: Position(ndim)
    !> Lower boundaries
    integer(int64), intent(in) :: LowerBoundaries(:,:)
    !> Upper boundaries
    integer(int64), intent(in) :: UpperBoundaries(:,:)
    
    integer(intmpi) :: proc
    
    do concurrent(proc=1:size(LowerBoundaries,2))
       if(all(LowerBoundaries(:,proc)<=Position)&
            .and.all(UpperBoundaries(:,proc)>=Position))then
          GetProc_fromGeneralPosition=proc-1
       end if
    end do
  end function GetProc_fromGeneralPosition

  !>@brief Deconstructs lattice index into position indices
  !!@returns position indices
  !!@details
  !! Lattice index is given by \n
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !! This relation is inverted by\n
  !! \f$ v_i = (index-1)/N_i+1\f$
  !! Then \f$index - v_i\f$ and that as input for \f$v_{i-1}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure function GetLatticePosition(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> lattice position
    integer(int64), dimension(ndim) :: GetLatticePosition

    GetLatticePosition = &
         GetPosition_FromIndex(LatticeIndex,LatticeExtensions)
  end function GetLatticePosition

  !>@brief Returns local lattice index based on global lattice index
  !!@returns local lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure elemental integer(int64) function GetMemoryIndex(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int64), intent(in) :: LatticeIndex

    GetMemoryIndex = FindLoc(&
                                ! Where to look
         LocalLatticeIndices,dim=1,&
                                ! What to look for
         value = LatticeIndex)
  end function GetMemoryIndex
  
  !>@brief Returns lattice index based on memory index
  !!@returns lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure elemental integer(int64) function GetLatticeIndex_M(MemoryIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Local index
    integer(int64), intent(in) :: MemoryIndex
    GetLatticeIndex_M = LocalLatticeIndices(MemoryIndex)
  end function GetLatticeIndex_M

  ! ..--** Generic Index Functions **--..
  !>@brief Returns  index of neighbouring point in i'th direction
  !!@returns lattice index of neighbouring point in i'th direction
  !!@details Periodic shift i'th direction
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure integer(int64) function GetNeib_FromIndex(i,index,extensions)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Direction
    integer(int8),  intent(in) :: i
    !> Index
    integer(int64), intent(in) :: index
    !> Extensions
    integer(int64), intent(in) :: extensions(ndim)

    integer(int64) :: position(ndim)

    if(i == 0_int8) then
       GetNeib_FromIndex = index
    else
       ! Getting position
       position = GetPosition_FromIndex(index,extensions)

       ! Periodicially shifting i'th index
       position(abs(i))&
            = GetPeriodicPosition(&
            position(abs(i)) + sign(1,i), & ! Shift
            Extensions(abs(i)))! Lattice extension, used for periodicity

       GetNeib_FromIndex&
            = GetIndex_FromPosition(position,extensions)
    end if

  end function GetNeib_FromIndex

  !>@brief Generic index function
  !!@returns index
  !!@details index is given by \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure integer(int64) function GetIndex_FromPosition(position,extensions)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> position vector
    integer(int64), intent(in) :: position(ndim)
    !> extensions
    integer(int64), intent(in) :: extensions(ndim)

    integer(int8) :: i

    GetIndex_FromPosition = 1
    do i=1,ndim
       GetIndex_FromPosition = &
            GetIndex_FromPosition + product(extensions(1:i-1))*(position(i)-1)
    end do
  end function GetIndex_FromPosition

  !>@brief
  !! Generic deconstruction function for index into position for given extensions
  !!@returns position indices
  !!@details Index is given by \n
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !! This relation is inverted by\n
  !! \f$ v_i = (index-1)/N_i+1\f$
  !! Then \f$index - v_i\f$ and that as input for \f$v_{i-1}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure function GetPosition_FromIndex(index,extensions)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> index
    integer(int64), intent(in) :: index
    !> extensions
    integer(int64), intent(in) :: extensions(ndim)
    !> position
    integer(int64), dimension(ndim) :: GetPosition_FromIndex

    integer(int64) :: ind
    integer(int8)  :: i

    ind = index
    do i=ndim,1,-1
       GetPosition_FromIndex(i) = &
            (ind-1_int64)/product(extensions(1:i-1))+1_int64
       ind = ind - (GetPosition_FromIndex(i)-1_int64)*product(extensions(1:i-1))
    end do
  end function GetPosition_FromIndex

  !>@brief Implements periodic boundary conditions
  !!@returns periodically projected position
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure elemental integer(int64) function GetPeriodicPosition(index,extension)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> component of position
    integer(int64), intent(in) :: index
    !> extension
    integer(int64), intent(in) :: extension
    GetPeriodicPosition = modulo(index-1_int64,extension)+1_int64
  end function GetPeriodicPosition

  ! ..--** Momenta **--..
  !>@brief Eigenvalue of momentum operator, defined via backwards derivative
  !!@returns Eigenvalue of momentum operator, defined via backwards derivative
  !!@details
  !! \f$\left(\hat{p}_k^{\text{F}}\right)_{m,n}
  !! =\frac{-\mathrm{i}}{a_k}\left(\delta_{m,n}-\delta_{m-1,n}\right)\f$
  !! with eigenvaues
  !! \f$p_{j,k}^{\text{F}}=\frac{2\sin\left(\frac{\pi j}{N_k}\right)}{a_k}
  !! e^{-\mathrm{i}\frac{\pi j}{N_k}}
  !! =\frac{2\sin\left(\frac{a_kp_{j,k}}{2}\right)}{a_k}e^{-\mathrm{i}\frac{a_kp_{j,k}}{2}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure function GetMomentum_BackwardDerivative(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    use mathconstants, only: pi
    implicit none
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Momentum
    complex(fp) :: GetMomentum_BackwardDerivative(ndim)

    integer(int64) :: LatticePosition(ndim)

    LatticePosition = GetLatticePosition(LatticeIndex)

    GetMomentum_BackwardDerivative = 2*sin(pi*(LatticePosition-1)/LatticeExtensions)&
         /LatticeSpacings(1:ndim)&
         *Exp(cmplx(0,-pi*(LatticePosition-1)/LatticeExtensions,fp))
  end function GetMomentum_BackwardDerivative

  !>@brief Eigenvalue of momentum operator, defined via forward derivative
  !!@returns Eigenvalue of momentum operator, defined via backwards derivative
  !!@details
  !! \f$\left(\hat{p}_k^{\text{F}}\right)_{m,n}
  !! =\frac{-\mathrm{i}}{a_k}\left(\delta_{m,n}-\delta_{m-1,n}\right)\f$
  !! with eigenvaues
  !! \f$p_{j,k}^{\text{F}}=\frac{2\sin\left(\frac{\pi j}{N_k}\right)}{a_k}
  !! e^{-\mathrm{i}\frac{\pi j}{N_k}}
  !! =\frac{2\sin\left(\frac{a_kp_{j,k}}{2}\right)}{a_k}e^{+\mathrm{i}\frac{a_kp_{j,k}}{2}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure function GetMomentum_ForwardDerivative(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    use mathconstants, only: pi
    implicit none
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Momentum
    complex(fp) :: GetMomentum_ForwardDerivative(ndim)

    integer(int64) :: LatticePosition(ndim)

    LatticePosition = GetLatticePosition(LatticeIndex)

    GetMomentum_ForwardDerivative = 2*sin(pi*(LatticePosition-1)/LatticeExtensions)&
         /LatticeSpacings(1:ndim)&
         *Exp(cmplx(0,+pi*(LatticePosition-1)/LatticeExtensions,fp))
  end function GetMomentum_ForwardDerivative

  !>@brief Eigenvalue of momentum operator, defined via central derivative
  !!@returns Eigenvalue of momentum operator, defined via central derivative
  !!@details
  !! \f$\left(\hat{p}_k^{\text{F}}\right)_{m,n}
  !! =\frac{-\mathrm{i}}{a_k}\left(\delta_{m+1,n}-\delta_{m-1,n}\right)\f$
  !! with eigenvaues
  !! \f$p_{j,k}^{\text{F}}=\frac{\sin\left(\frac{2\pi j}{N_k}\right)}{a_k}
  !! =\frac{\sin\left(a_kp_{j,k}\right)}{a_k}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure function GetMomentum_CentralDerivative(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use mathconstants, only: pi
    implicit none
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Momentum
    complex(fp) :: GetMomentum_CentralDerivative(ndim)

    integer(int64) :: LatticePosition(ndim)

    LatticePosition = GetLatticePosition(LatticeIndex)

    GetMomentum_CentralDerivative = sin(2*pi*(LatticePosition-1)/LatticeExtensions)&
         /LatticeSpacings(1:ndim)
  end function GetMomentum_CentralDerivative

  !>@brief Momentum from memory index
  !!@returns Momentum from memory index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  pure function GetMomentum_M(MemoryIndex)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    implicit none
    !> Memory Index
    integer(int64), intent(in)  :: MemoryIndex
    !> Momentum
    complex(fp) :: GetMomentum_M(nDim)

    GetMomentum_M = GetMomentum_G(GetLatticeIndex_M(MemoryIndex))
  end function GetMomentum_M

   !>@brief 2-Norm of lattice momentum from memory index
  !!@returns 2-Norm of lattice momentum
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure real(fp) function GetNorm2Momentum_M(MemoryIndex)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    implicit none
    !> Memory Index
    integer(int64), intent(in) :: MemoryIndex
    GetNorm2Momentum_M = GetNorm2Momentum_G(GetLatticeIndex_M(MemoryIndex))
  end function GetNorm2Momentum_M
  
  !>@brief 2-Norm of lattice momentum
  !!@returns 2-Norm of lattice momentum
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure real(fp) function GetNorm2Momentum_G(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> lattice index
    integer(int64), intent(in) :: LatticeIndex
    complex(fp) :: Momentum(ndim)
    Momentum = GetMomentum_G(LatticeIndex)
    GetNorm2Momentum_G = norm2(abs(Momentum))
  end function GetNorm2Momentum_G

  !>@brief
  !! Returns biggest mometum norm on the lattice
  !!@returns
  !! Biggest mometum norm on the lattice
  !!@author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
 pure real(fp) function GetMaxNorm2Momentum()
    use, intrinsic :: iso_fortran_env
    implicit none
    GetMaxNorm2Momentum = MaxNorm2Momentum
  end function GetMaxNorm2Momentum

  ! ..--** Polarisation Vectors **--..
  !>@brief Returns polarisation vectors
  !!@details The polarisation vectors are constructed based on the
  !! lattice momentum index \f$\vec{p}\f$ with minimum momentum position \f$\vec{0}\f$.
  !! The choice of polarisation vectors is based on the paper from Kaspar and Hebenstreit in
  !! <a href="https://arxiv.org/abs/1403.4849">arxiv:1403.4849</a> 
  !!@author Alexander Lehmann, !UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  !!@todo Finding a dimension-independend way of programming this -- or at least one which compiles
  !! for different numbers of dimensions. For now one would have to overwrite this code when
  !! one wants to change the numbers of dimensions
 pure subroutine GetPolarisationVectors(p,PolarisationVectors)
    use precision, only: fp
    use tolerances, only: GetZeroTol
    implicit none
    !> Lattice momentum vector
    complex(fp), intent(in) :: p(nDim)
    !> Polarisations
    complex(fp), intent(out) :: PolarisationVectors(nDim,nDim)

    real(fp) :: abs_p
    
    select case(nDim)
    case(1)
       PolarisationVectors(1,1) = 1 !(only one component)
    case(3)
       abs_p = norm2(abs(p))
       if(abs_p < GetZeroTol()) then
          PolarisationVectors = 0
       elseif( abs(p(1)) > GetZeroTol()) then
          ! First transversal polarisation vector
          PolarisationVectors(:,1) = [&
               -p(2)/norm2(abs(p(1:2))),& ! First component
               +p(1)/norm2(abs(p(1:2))),& ! Second component
               cmplx(0,0,fp)&             ! Third component
               ]

          ! Second transversal polarisation vector
          PolarisationVectors(:,2) = [&
               conjg(p(1))*p(3)/(abs_p*norm2(abs(p(1:2)))),& ! First component
               conjg(p(2))*p(3)/(abs_p*norm2(abs(p(1:2)))),& ! Second component
               cmplx(-norm2(abs(p(1:2)))/abs_p,0,fp)&        ! Third component
               ]

          ! Third == longitudinal polarisation vector
          PolarisationVectors(:,3) = p/abs_p
       else
          ! First transversal polarisation vector
          PolarisationVectors(:,1) = [&
               cmplx(0,0,fp),&            ! First component
               +p(3)/norm2(abs(p(2:3))),& ! Second component
               -p(2)/norm2(abs(p(2:3)))&  ! Third component
               ]
          
          ! Second transversal polarisation vector
          PolarisationVectors(:,2) = [&
               cmplx(1,0,fp),& ! First component
               cmplx(0,0,fp),& ! Second comonent
               cmplx(0,0,fp) & ! Third componetn
               ]

          ! Third == longitudinal polarisation vector
          PolarisationVectors(:,3) = p/abs_p
       end if
    case default
       PolarisationVectors = 0
    end select
  end subroutine GetPolarisationVectors

  !>@brief Transverse projection operator on the lattice
  !!@returns Transverse projection operator on the lattice
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
 pure real(fp) function GetTransverseProjector(i,j,momentum)
    use tolerances, only: GetZeroTol
    implicit none
    !> Direction
    integer(int8), intent(in) :: i,j
    !> Lattice momentum
    complex(fp),   intent(in) :: momentum(ndim)

    if(norm2(abs(momentum)) < GetZeroTol()) then
       GetTransverseProjector = 0
    else
       GetTransverseProjector = &
            (GetKroneckerDelta(i,j) - conjg(momentum(i))*momentum(j)/norm2(abs(momentum))**2)/2
    end if
  end function GetTransverseProjector

  ! ..--** Auxiliary Mathematical Routines **--..
  
  !>@brief Kronecker-Delta
  !!@returns Kronecker-Delta \f$\delta_{i,j}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 26.11.2018
  !!@version 1.0
 pure elemental integer(int8) function GetKroneckerDelta(i,j)
    implicit none
    !> Index \f$i\f$
    integer(int8), intent(in) :: i
    !> Index \f$j\f$
    integer(int8), intent(in) :: j

    if(i==j) then
       GetKroneckerDelta = 1
    else
       GetKroneckerDelta = 0
    end if
  end function GetKroneckerDelta
end module lattice
