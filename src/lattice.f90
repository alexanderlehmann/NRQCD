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
       GetGlobalLatticeIndex,&
       GetLocalIndex, &
       GetLocalLatticeIndices,&
       GetLocalLatticeIndices_includingHalo,&
       GetLocalLatticeIndices_allocatable,&
       GetLocalLatticeIndices_includingHalo_allocatable,&
       GetLatticePosition,&
       GetLatticeIndex,&
       GetProc,&
       GetProc_fromGeneralIndex,&
       GetLatticeExtension,&
       GetLatticeSize,&
       GetLocalLatticeSize,&
       InitLatticeIndices,&
       GetLocalLatticeSize_includingHalo,&
       GetLatticeSpacing,&
       GetVolume,&
       GetNeib,&
       GetMomentum,&
       GetNorm2Momentum,&
       GetMaxNorm2Momentum

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
  !> lattice spacings
  real(fp)   :: LatticeSpacings(0:ndim) = -1._fp
  !> volume
  real(fp)   :: Volume = -1._fp
  !> local lattice sizes
  integer(int64) :: LocalLatticeSize=-1
  !> local lattice sizes including halo
  integer(int64) :: LocalLatticeSize_includingHalo=-1

  !> maximum 2-norm of momentum
  real(fp) :: MaxNorm2Momentum=-1._fp

  !..--** MPI: process-dependend parameters**--..
  !> local lower lattice boundaries
  integer(int64), allocatable :: LocalLowerLatticeBoundaries(:,:)
  !> local upper lattice boundaries
  integer(int64), allocatable :: LocalUpperLatticeBoundaries(:,:)
  !> local lattice indices
  integer(int64), allocatable :: LocalLatticeIndices(:)
  !> local lattice size including halo
  integer(int64), allocatable :: LocalLatticeSizes_IncludingHalo(:)
  !> local lower lattice boundaries including halo
  integer(int64), allocatable :: LocalLowerLatticeBoundaries_IncludingHalo(:,:)
  !> local upper lattice boundaries including halo
  integer(int64), allocatable :: LocalUpperLatticeBoundaries_IncludingHalo(:,:)
  !> local lattice indices including halo
  integer(int64), allocatable :: LocalLatticeIndices_includingHalo(:)

  !>@brief
  !! Lattice momentum
  !!@returns
  !! Lattice momentum
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date
  !! 17.02.2019
  !!@version
  !! 1.0
  interface GetMomentum
     module procedure GetMomentum_BackwardDerivative
     !module procedure GetMomentum_CentralDerivative
     !module procedure GetMomentum_ForwardDerivative
  end interface GetMomentum

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
    
    character(len=100) :: errormessage

    if(isInitialised) then
       errormessage = 'Error in init of '//modulename//': already initialised.'
       call MPISTOP(errormessage)
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
            LocalLatticeIndices,                        &
            LocalLowerLatticeBoundaries(:,ThisProc()+1),&
            LocalUpperLatticeBoundaries(:,ThisProc()+1))
       call InitLatticeIndices(&
            LocalLatticeIndices_includingHalo,&
            LocalLowerLatticeBoundaries_includingHalo(:,ThisProc()+1),&
            LocalUpperLatticeBoundaries_includingHalo(:,ThisProc()+1))

       LocalLatticeSize_includingHalo = size(LocalLatticeIndices_includingHalo)

       ! Deleting lattice points which are contained more than once ...
       ! (important in case of small partitions)
       call RemoveDuplicates(&
            LocalLatticeIndices)
       call RemoveDuplicates(&
            LocalLatticeIndices_includingHalo)
       ! ... and sort the arrays
       call Sort(LocalLatticeIndices)
       call Sort(LocalLatticeIndices_includingHalo)

       ! Setting local lattice size (including and without halo)
       LocalLatticeSize = size(LocalLatticeIndices)
       LocalLatticeSize_includingHalo = size(LocalLatticeIndices_includingHalo)

       ! Find maximum momentum
       MaxNorm2Momentum = FindMaxNorm2Momentum()

       ! DONE
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule
  
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
    character(len=70) :: errormessage
    
    if(.not.IsMPIinterfaceInitialised()) then
       errormessage = 'Error in init of '//modulename//': '//mpiname//' is not initialised.'
       call mpistop(errormessage)
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
    use mpiinterface, only: NumProcs, intmpi, MPIstop
    implicit none

    integer(int64) :: LocalLatticeindex, LatticeIndex
    real(fp)   :: norm2momentum,currentMax
    real(fp), allocatable :: maxmomenta(:)

    integer(intmpi) :: mpierr, sendtype
    
    character(len=100) :: errormessage
    
    currentMax = 0
    do LocalLatticeIndex=1,LocalLatticeSize
       LatticeIndex = LocalLatticeIndices(LocalLatticeIndex)
       norm2momentum = GetNorm2Momentum(LatticeIndex)
       if(norm2momentum.gt.currentMax) currentMax = norm2Momentum
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
       errormessage = 'Error in FindMaxNorm2Momentum of '//modulename&
            //': unsupported floating point precision.'
       call MPISTOP(errormessage)
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
         = GetIndex_FromPosition(Position,LatticeExtensions)
  end function GetLatticeIndex

  !>@brief Returns lattice index of neighbouring point in i'th direction
  !!@returns lattice index of neighbouring point in i'th direction
  !!@details Periodic shift i'th direction
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure integer(int64) function GetNeib(i,LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: LatticeIndex

    integer(int64) :: LatticePosition(ndim)

    GetNeib &
         = GetNeib_FromIndex(i,LatticeIndex,LatticeExtensions)
  end function GetNeib

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
    allocate(indices(LocalLatticeSize))
    indices = LocalLatticeIndices
  end subroutine GetLocalLatticeIndices_Allocatable

  !>@brief Returns local lattice indices including halo
  !!@returns local lattice indices including halo
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure subroutine GetLocalLatticeIndices_includingHalo(indices)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> local lattice indices
    integer(int64), intent(out) :: indices(:)
    indices = LocalLatticeIndices_includingHalo
  end subroutine GetLocalLatticeIndices_IncludingHalo

  !>@brief Returns local lattice indices including halo with allocation
  !!@returns local lattice indices including halo
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure subroutine GetLocalLatticeIndices_includingHalo_Allocatable(indices)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> local lattice indices
    integer(int64), intent(out), allocatable :: indices(:)
    allocate(indices(LocalLatticeSize_includingHalo))
    indices = LocalLatticeIndices_includingHalo
  end subroutine GetLocalLatticeIndices_IncludingHalo_Allocatable

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

  !>@brief Returns local lattice size without halo
  !!@returns local lattice size without halo
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

  !>@brief Returns local lattice size including halo
  !!@returns local lattice size including halo
  !!@details Lattice size is \f$|\Lambda|=\prod\limits_{i=1}^{d}N_i\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure integer(int64) function GetLocalLatticeSize_includingHalo()
    use, intrinsic :: iso_fortran_env
    implicit none
    GetLocalLatticeSize_includingHalo = LocalLatticeSize_includingHalo
  end function GetLocalLatticeSize_includingHalo

  !>@brief Returns MPI-process-rank corresponding to lattice index
  !!@returns MPI-process-rank corresponding to lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.1
  pure elemental integer(intmpi) function GetProc(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: intmpi
    implicit none
    !> lattice index
    integer(int64), intent(in) :: LatticeIndex
    GetProc = GetProc_fromGeneralIndex(LatticeIndex,&
         LocalLowerLatticeBoundaries,LocalUpperLatticeBoundaries,LatticeExtensions)
  end function GetProc

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

    do concurrent(proc=1:size(LowerBoundaries,2))
       if(all(LowerBoundaries(:,proc)<=Position)&
            .and.all(UpperBoundaries(:,proc)>=Position))then
          GetProc_fromGeneralIndex=proc-1
       end if
    end do
  end function GetProc_fromGeneralIndex

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
  pure elemental integer(int64) function GetLocalIndex(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int64), intent(in) :: LatticeIndex

    GetLocalIndex = FindLoc(&
                                ! Where to look
         LocalLatticeIndices_includingHalo,dim=1,&
                                ! What to look for
         value = LatticeIndex)
  end function GetLocalIndex
  
  !>@brief Returns lattice index based on local index
  !!@returns lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure elemental integer(int64) function GetGlobalLatticeIndex(LocalIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int64), intent(in) :: LocalIndex
    GetGlobalLatticeIndex = LocalLatticeIndices_IncludingHalo(LocalIndex)
  end function GetGlobalLatticeIndex

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

  !>@brief 2-Norm of lattice momentum
  !!@returns 2-Norm of lattice momentum
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure real(fp) function GetNorm2Momentum(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    implicit none
    !> lattice index
    integer(int64), intent(in) :: LatticeIndex
    complex(fp) :: Momentum(ndim)
    Momentum = GetMomentum(LatticeIndex)
    GetNorm2Momentum = norm2(abs(Momentum))
  end function GetNorm2Momentum

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
end module lattice
