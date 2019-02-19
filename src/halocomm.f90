!----------------------------------------------------------------------
! Halo communcation module
!----------------------------------------------------------------------
!
! MODULE: halocomm
!> @brief Communication of halo values
!! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date 17.02.2019
!! @version 1.0
! REVISION HISTORY:
! 17 02 2019 - Initial version
!----------------------------------------------------------------------
module halocomm
  use, intrinsic :: iso_fortran_env
  implicit none

  PRIVATE

  public :: &
       InitModule,&
       IsModuleInitialised,&
       CommunicateBoundary

  !> Module name
  character(len=8), parameter, public ::  modulename='halocomm'
  
  !> Contains information, if module is initialised
  logical :: IsInitialised = .false.

  !> List of MPI-ranks of neighbours
  integer,        allocatable :: HaloProcs(:)
  !> Neighbouring points for each neighbour
  integer,        allocatable :: NeibPoints(:)
  !> List of points and from this process recieving processes
  integer(int64), allocatable :: SendList(:,:)
  !> List of points and to this process sending processes
  integer(int64), allocatable :: RecvList(:,:)
  !> Number of neighbours
  integer :: neibs

  !> @brief Communication of boundary values
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  interface CommunicateBoundary
     ! ..--** REAL **--..
     ! real(real32) aka single precision
     module procedure CommunicateBoundary_real_single_rank1
     module procedure CommunicateBoundary_real_single_rank2
     module procedure CommunicateBoundary_real_single_rank3
     ! real(real64) aka double precision
     module procedure CommunicateBoundary_real_double_rank1
     module procedure CommunicateBoundary_real_double_rank2
     module procedure CommunicateBoundary_real_double_rank3
     ! real(real128) aka quad precision
     module procedure CommunicateBoundary_real_quad_rank1
     module procedure CommunicateBoundary_real_quad_rank2
     module procedure CommunicateBoundary_real_quad_rank3

     ! ..--** COMPLEX **--..
     ! complex(real32) aka complex single precision
     module procedure CommunicateBoundary_complex_single_rank1
     module procedure CommunicateBoundary_complex_single_rank2
     module procedure CommunicateBoundary_complex_single_rank3
     module procedure CommunicateBoundary_complex_single_rank4
     ! complex(real64) aka complex double precision
     module procedure CommunicateBoundary_complex_double_rank1
     module procedure CommunicateBoundary_complex_double_rank2
     module procedure CommunicateBoundary_complex_double_rank3
     module procedure CommunicateBoundary_complex_double_rank4
     ! complex(real128) aka complex quad precision
     module procedure CommunicateBoundary_complex_quad_rank1
     module procedure CommunicateBoundary_complex_quad_rank2
     module procedure CommunicateBoundary_complex_quad_rank3
     module procedure CommunicateBoundary_complex_quad_rank4
  end interface CommunicateBoundary
contains
  !> @brief Initialises module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    implicit none

    call CheckDependencies
    
    ! Initialise list of lattice points which are to be recieved from which other process
    call InitSendRecvLists(Neibs,HaloProcs,NeibPoints,SendList,RecvList)

    ! DONE
    IsInitialised = .TRUE.
    
  end subroutine InitModule
  
  !> @brief Checks previous necessary initialisations
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine CheckDependencies
    use, intrinsic :: iso_fortran_env
    use lattice, only: IsLatticeInitialised => IsModuleInitialised
    use mpi
    implicit none

    integer :: proc, mpierr
    if(.not.IsLatticeInitialised()) then
       call mpi_comm_rank(MPI_COMM_WORLD, proc, mpierr)
       if(proc==0) then
          call flush(ERROR_UNIT)
          write(ERROR_UNIT,*) 'Error in init of ',modulename&
               ,': Lattice-module is not initialised.'
          call flush(ERROR_UNIT)
       end if
       STOP
    end if
  end subroutine CheckDependencies
  
  !>@brief Returns, if module is initialised
  !! @returns module's initialisation status
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  pure logical function IsModuleInitialised()
    use, intrinsic :: iso_fortran_env
    implicit none
    IsModuleInitialised = isInitialised
  end function IsModuleInitialised

  !>@brief Initialises list of processes and points which send to this process
  !! @returns list of processes and points which send to this process
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine InitSendRecvLists(Neibs,HaloProcs,NeibPoints,SendList,RecvList)
    use, intrinsic :: iso_fortran_env
    use lattice, only: GetLocalLatticeIndices_includingHalo_Allocatable, GetProc
    use mpiinterface, only: NumProcs, ThisProc, MPISTOP
    use arrayoperations, only: RemoveDuplicates, Sort
    use mpi
    implicit none
    !> Number of neighbours
    integer,                     intent(out) :: Neibs
    !> List of MPI-ranks of neighbours
    integer,        allocatable, intent(out) :: HaloProcs(:)
    !> Neighbouring points for each neighbour
    integer,        allocatable, intent(out) :: NeibPoints(:)
    !> List of points and from this process recieving processes
    integer(int64), allocatable, intent(out) :: SendList(:,:)
    !> List of points and to this process sending processes
    integer(int64), allocatable, intent(out) :: RecvList(:,:)

    integer(int64), allocatable :: LocalLatticeIndices(:)
    integer(int64), allocatable :: PointsPerProc_includingThisProc(:)
    integer, allocatable :: HaloProcs_includingThisProc(:)
    integer, allocatable :: Procs(:)
    integer(int64) :: LocalIndex, neibpoint
    integer :: proc, neib
    integer(int64) :: MaxHaloPoints

    ! MPI
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr
    
    ! 1. Compute halo points and the associated process numbers
    call GetLocalLatticeIndices_includingHalo_Allocatable(LocalLatticeIndices)

    allocate(Procs(NumProcs()))
    allocate(HaloProcs_includingThisProc(NumProcs()))
    Procs     = GetProc(LocalLatticeIndices)
    HaloProcs_includingThisProc = Procs
    call RemoveDuplicates(HaloProcs_includingThisProc,PointsPerProc_includingThisProc)
    call Sort(HaloProcs_includingThisProc)
    
    neibs = size(HaloProcs_includingThisProc)-1
    allocate(HaloProcs(neibs))
    allocate(NeibPoints(neibs))
    MaxHaloPoints = MaxVal(&
         PointsPerProc_includingThisProc,DIM=1,& !Where to look
         MASK=HaloProcs_includingThisProc/=ThisProc() ) !What to ignore

    allocate(RecvList(MaxHaloPoints,neibs))
    RecvList = -1
    
    neib = 0
    do proc=1,size(HaloProcs_includingThisProc)
       if(HaloProcs_includingThisProc(proc) /= ThisProc()) then
          neib = neib + 1
          HaloProcs(neib)  = HaloProcs_includingThisProc(proc)
          NeibPoints(neib) = PointsPerProc_includingThisProc(proc)

          ! Assign global lattice indices of halo to the process which sends them
          neibpoint=0
          do LocalIndex=1,size(LocalLatticeIndices)
             if(GetProc(LocalLatticeIndices(LocalIndex))==HaloProcs(neib)) then
                neibpoint = neibpoint + 1_int64
                RecvList(neibpoint,neib) = LocalLatticeIndices(LocalIndex)
             end if
          end do
       end if
    end do

    ! 2. Get from all the neighbours the lists of points they wish to recieve
    !    using the symmetry isSender <=> isReciever

    allocate(SendList(MaxHaloPoints,neibs))
    do neib=1,neibs
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            RecvList(&          ! What to send ...
            1,neib),&           ! ... and it's first index
            NeibPoints(neib),&  ! How many points
            MPI_INT64_T,&       ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            SendList(&          ! What to recieve ...
            1,neib),&           ! ... and it's first index
            NeibPoints(neib),&  ! How many points
            MPI_INT64_T,&       ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code
    end do
  end subroutine InitSendRecvLists

  !>@brief Communication of 1D-array with real entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_single_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    real(kind), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(BufferIndex) = data(LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1),&                ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL4,&         ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1),&                ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL4,&         ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(LocalIndex) = buffer(BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_single_rank1

  !>@brief Communication of 2D-array with real entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_single_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    real(kind), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,BufferIndex) = data(:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1,1),&              ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL4,&         ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1,1),&              ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL4,&         ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,LocalIndex) = buffer(:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_single_rank2

  !>@brief Communication of 3D-array with real entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_single_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    real(kind), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,BufferIndex) = data(:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1,1,1),&            ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL4,&         ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1,1,1),&            ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL4,&         ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,LocalIndex) = buffer(:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_single_rank3
  
  !>@brief Communication of 1D-array with real entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_double_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    real(kind), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(BufferIndex) = data(LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1),&                ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL8,&         ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1),&                ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL8,&         ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(LocalIndex) = buffer(BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_double_rank1

  !>@brief Communication of 2D-array with real entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_double_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    real(kind), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,BufferIndex) = data(:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1,1),&              ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL8,&         ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1,1),&              ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL8,&         ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,LocalIndex) = buffer(:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_double_rank2

  !>@brief Communication of 3D-array with real entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_double_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    real(kind), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,BufferIndex) = data(:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1,1,1),&            ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL8,&         ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1,1,1),&            ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL8,&         ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,LocalIndex) = buffer(:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_double_rank3

  !>@brief Communication of 1D-array with real entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_quad_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    real(kind), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(BufferIndex) = data(LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1),&                ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL16,&        ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1),&                ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL16,&        ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(LocalIndex) = buffer(BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_quad_rank1

  !>@brief Communication of 2D-array with real entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_quad_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    real(kind), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,BufferIndex) = data(:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1,1),&              ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL16,&        ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1,1),&              ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL16,&        ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,LocalIndex) = buffer(:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_quad_rank2

  !>@brief Communication of 3D-array with real entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_real_quad_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    real(kind), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(kind), allocatable :: buffer(:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,BufferIndex) = data(:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&            ! What to send ...
            1,1,1),&            ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL16,&        ! What type to send
            dest, sendtag,&     ! Destination and sendtag
            buffer(&            ! What to recieve ...
            1,1,1),&            ! ... and it's first index
            buffersize,&        ! How many points
            MPI_REAL16,&        ! What type to recieve
            src,  recvtag,&     ! Source and recvtag
            mpi_comm_world,&    ! Communicator
            status, mpierr)     ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,LocalIndex) = buffer(:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_real_quad_rank3

  




  ! START COMPLEX

  !>@brief Communication of 1D-array with complex entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_single_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    complex(kind), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(BufferIndex) = data(LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1),&             ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1),&             ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(LocalIndex) = buffer(BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_single_rank1

  !>@brief Communication of 2D-array with complex entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_single_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    complex(kind), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,BufferIndex) = data(:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1),&           ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1),&           ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,LocalIndex) = buffer(:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_single_rank2

  !>@brief Communication of 3D-array with complex entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_single_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    complex(kind), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,BufferIndex) = data(:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1,1),&         ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1,1),&         ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,LocalIndex) = buffer(:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_single_rank3

  !>@brief Communication of 4D-array with complex entries of single precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_single_rank4(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real32
    complex(kind), intent(inout) :: data(:,:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),size(data,3),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,:,BufferIndex) = data(:,:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1,1,1),&       ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1,1,1),&       ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX8,&   ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,:,LocalIndex) = buffer(:,:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_single_rank4
  
  !>@brief Communication of 1D-array with complex entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_double_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    complex(kind), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(BufferIndex) = data(LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1),&             ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1),&             ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(LocalIndex) = buffer(BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_double_rank1

  !>@brief Communication of 2D-array with complex entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_double_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    complex(kind), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,BufferIndex) = data(:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1),&           ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1),&           ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,LocalIndex) = buffer(:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_double_rank2

  !>@brief Communication of 3D-array with complex entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_double_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    complex(kind), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,BufferIndex) = data(:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1,1),&         ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1,1),&         ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,LocalIndex) = buffer(:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_double_rank3

  !>@brief Communication of 4D-array with complex entries of double precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_double_rank4(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real64
    complex(kind), intent(inout) :: data(:,:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),size(data,3),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,:,BufferIndex) = data(:,:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1,1,1),&       ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1,1,1),&       ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,:,LocalIndex) = buffer(:,:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_double_rank4

  !>@brief Communication of 1D-array with complex entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_quad_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    complex(kind), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(BufferIndex) = data(LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1),&             ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1),&             ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(LocalIndex) = buffer(BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_quad_rank1

  !>@brief Communication of 2D-array with complex entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_quad_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    complex(kind), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,BufferIndex) = data(:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1),&           ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1),&           ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,LocalIndex) = buffer(:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_quad_rank2

  !>@brief Communication of 3D-array with complex entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_quad_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    complex(kind), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,BufferIndex) = data(:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1,1),&         ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1,1),&         ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,LocalIndex) = buffer(:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_quad_rank3

  !>@brief Communication of 4D-array with complex entries of quad precision
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 18.02.2019
  !! @version 1.0
  impure subroutine CommunicateBoundary_complex_quad_rank4(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetLocalIndex
    use mpiinterface, only: ThisProc
    implicit none
    integer(int8), parameter :: kind = real128
    complex(kind), intent(inout) :: data(:,:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(kind), allocatable :: buffer(:,:,:,:)
    integer :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr

    
    integer :: neib, proc
    integer :: ValuesPerPoint, MaxHaloPoints
    integer :: BufferIndex
    integer(int64) :: LocalIndex
    integer(int64) :: LatticeIndex

    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(buffer(size(data,1),size(data,2),size(data,3),MaxHaloPoints))

    ValuesPerPoint = size(data)/size(data,r)
    
    AllNeighbours: do neib=1,Neibs
       buffersize = ValuesPerPoint*NeibPoints(neib)

       ! Prepare buffer with send-data
       packing: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          buffer(:,:,:,BufferIndex) = data(:,:,:,LocalIndex)
       end do packing

       ! Send and recieve
       dest = HaloProcs(neib)
       src  = dest

       sendtag = ThisProc()
       recvtag = dest

       call MPI_SendRecv(&
            buffer(&         ! What to send ...
            1,1,1,1),&       ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            buffer(&         ! What to recieve ...
            1,1,1,1),&       ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX32,&  ! What type to recieve
            src,  recvtag,&  ! Source and recvtag
            mpi_comm_world,& ! Communicator
            status, mpierr)  ! Status and error-code

       ! Unpack recieved data from buffer
       unpacking: do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)
          
          ! Get local index in data
          LocalIndex = GetLocalIndex(LatticeIndex)

          ! Assign data to buffer
          data(:,:,:,LocalIndex) = buffer(:,:,:,BufferIndex)
       end do unpacking
    end do AllNeighbours
    deallocate(buffer)
  end subroutine CommunicateBoundary_complex_quad_rank4
end module halocomm
