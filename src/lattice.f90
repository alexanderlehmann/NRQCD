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
  integer(int64) :: LatticeExtensions(ndim) = -1
  !> lattice size
  integer(int64) :: LatticeSize = -1
  !> lattice spacings
  real(real64)   :: LatticeSpacings(0:ndim) = -1._real64
  !> volume
  real(real64)   :: Volume = -1._real64
  !> local lattice sizes
  integer(int64) :: LocalLatticeSize=-1
  !> local lattice sizes including halo
  integer(int64) :: LocalLatticeSize_includingHalo=-1
  !> local volume
  real(real64)   :: LocalVolume  = -1._real64

  !> global maximum 2-norm of momentum
  real(real64) :: maxmomentum=-1._real64

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

  impure subroutine InitModule(LatticeExtensions_,LatticeSpacings_)
    use mpi
    use mpiinterface, only: ThisProc, NumProcs, MPISTOP

    use arrayoperations, only: RemoveDuplicates
    implicit none
    !> lattice extensions
    integer(int64), intent(in) :: LatticeExtensions_(ndim)
    !> lattice spacings
    real(real64),   intent(in) :: LatticeSpacings_(0:ndim)

    integer(int64) :: partitions(ndim)
    integer :: proc
    integer(int8) :: log2_latticeext(ndim), idivision, numdivisions, ipartition
    integer(int64) :: latticeindex

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
    
    ! Deleting lattice points which are contained more than once
    ! (important in case of small partitions)
    call RemoveDuplicates(&
         LocalLatticeIndices)
    call RemoveDuplicates(&
         LocalLatticeIndices_includingHalo)
    ! Setting local lattice size (including and without halo)
    LocalLatticeSize = size(LocalLatticeIndices)
    LocalLatticeSize_includingHalo = size(LocalLatticeIndices_includingHalo)

    !print*,LocalLatticeSize
    
    !if(ThisProc()==0) then
    !  do proc=0,NumProcs()-1
    !      print*,'Process',proc
    !      print*,LocalLowerLatticeBoundaries_includingHalo(:,proc+1)
    !      print*,LocalLowerLatticeBoundaries(:,proc+1)
    !      print*,LocalUpperLatticeBoundaries(:,proc+1)
    !      print*,LocalUpperLatticeBoundaries_includingHalo(:,proc+1)
    !      call flush(6)
    !   end do
    !end if

    
  contains

    pure recursive subroutine InitLocalLatticeBoundaries(&
         LocalLowerLatticeBoundaries,LocalUpperLatticeBoundaries,partitions,&
         proc,i,ipart,nHalo)
      implicit none
      integer(int64), intent(inout), allocatable :: LocalLowerLatticeBoundaries(:,:)
      integer(int64), intent(inout), allocatable :: LocalUpperLatticeBoundaries(:,:)
      integer(int64), intent(in)    :: partitions(ndim)
      integer(int64), intent(inout), optional :: proc
      integer(int8),  intent(in),    optional :: i
      integer(int64), intent(in),    optional :: ipart(ndim)
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


    pure recursive subroutine InitLatticeIndices(LatticeIndices,&
         LowerBoundaries,UpperBoundaries,i,x,index)
      implicit none
      integer(int64), intent(inout), allocatable :: LatticeIndices(:)
      integer(int64), intent(in)                 :: LowerBoundaries(ndim),UpperBoundaries(ndim)
      integer(int8),  intent(in),    optional :: i
      integer(int64), intent(in),    optional :: x(ndim)
      integer(int64), intent(inout), optional :: index

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
    
  end subroutine InitModule

  !> @brief
  !! Returns lattice index
  !! @returns
  !! lattice index
  !! @details
  !! Lattice index is given by
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 17.02.2019
  !! @version
  !! 1.0
  pure integer(int64) function GetLatticeIndex(Position)
    implicit none
    !> position vector
    integer(int64), intent(in) :: Position(ndim)
    
    GetLatticeIndex &
         = GetIndex_FromPosition(Position,LatticeExtensions)
  end function GetLatticeIndex






  ! ..--** Generic Index Functions **--..
  !> @brief
  !! Returns  index of neighbouring point in i'th direction
  !! @returns
  !! lattice index of neighbouring point in i'th direction
  !! @details
  !! Periodic shift i'th direction
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 16.02.2019
  !! @version
  !! 1.0
  pure integer(int64) function GetNeib_FromIndex(i,index,extensions)
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

  !> @brief
  !! Generic index function
  !! @returns
  !! index
  !! @details
  !! Lattice index is given by
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 16.02.2019
  !! @version
  !! 1.0
  pure integer(int64) function GetIndex_FromPosition(position,extensions)
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
  
  !> @brief
  !! Generic deconstruction function for index into position for given extensions
  !! @returns
  !! position indices
  !! @details
  !! Index is given by \n
  !! \f$index(\vec{v}) = v_i+N_i\cdot index(\vec{v}_{i+1..ndim})-1\f$\n
  !! This relation is inverted by\n
  !! \f$ v_i = (index-1)/N_i+1\f$
  !! Then \f$index - v_i\f$ and that as input for \f$v_{i-1}\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 16.02.2019
  !! @version
  !! 1.0
  pure function GetPosition_FromIndex(index,extensions)
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

  !> @brief
  !! Implements periodic boundary conditions
  !! @returns
  !! periodically projected position
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 17.02.2019
  !! @version
  !! 1.0
  pure elemental integer(int64) function GetPeriodicPosition(index,extension)
    implicit none
    !> component of position
    integer(int64), intent(in) :: index
    !> extension
    integer(int64), intent(in) :: extension
    GetPeriodicPosition = modulo(index-1_int64,extension)+1_int64
  end function GetPeriodicPosition

end module lattice
