!----------------------------------------------------------------------
! Halo communcation module
!----------------------------------------------------------------------
!
! MODULE: halocomm
!>@brief Communication of halo values
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 17.02.2019
!!@version 1.0
! REVISION HISTORY:
! 17 02 2019 - Initial version
!----------------------------------------------------------------------
module halocomm
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  use mpiinterface, only: intmpi
  implicit none

  PRIVATE

  public :: &
       InitModule,&
       IsModuleInitialised,&
       CommunicateBoundary

  !> Module name
  character(len=8), parameter, public ::  modulename='halocomm'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.

  !> List of MPI-ranks of neighbours
  integer(intmpi), allocatable :: HaloProcs(:)
  !> Number of neighbouring points for each neighbour
  integer(intmpi), allocatable :: NeibPoints(:)
  !> List of points and from this process recieving processes
  integer(int64), allocatable :: SendList(:,:)
  !> List of points and to this process sending processes
  integer(int64), allocatable :: RecvList(:,:)
  !> Number of neighbours
  integer(intmpi) :: neibs=-1
  
  !>@brief Communication of boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
  interface CommunicateBoundary
     ! ..--** REAL **--..
     module procedure CommunicateBoundary_real_rank1
     module procedure CommunicateBoundary_real_rank2
     module procedure CommunicateBoundary_real_rank3
     module procedure CommunicateBoundary_real_rank4

     ! ..--** COMPLEX **--..
     module procedure CommunicateBoundary_complex_rank1
     module procedure CommunicateBoundary_complex_rank2
     module procedure CommunicateBoundary_complex_rank3
     module procedure CommunicateBoundary_complex_rank4
  end interface CommunicateBoundary
  
contains ! Module procedures
  
  !>@brief Initialises module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: mpistop
    use mpi, only: &
         MPI_REAL4,    MPI_REAL8,     MPI_REAL16,&
         MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32
    implicit none

    if(isInitialised) then
       call MPISTOP('Error in init of '//modulename//': already initialised.')
    else
       
       call CheckObligatoryInitialisations

       ! Initialise list of lattice points which are to be recieved from which other process
       call InitSendRecvLists(Neibs,HaloProcs,NeibPoints,SendList,RecvList)
       
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
    use lattice, only: IsLatticeInitialised => IsModuleInitialised, latticename => modulename
    use mpiinterface, only: mpistop
    implicit none
    
    if(.not.IsLatticeInitialised()) then
       call mpistop('Error in init of '//modulename//': '//latticename//' is not initialised.')
    end if
  end subroutine CheckObligatoryInitialisations
  
  !>@brief Returns, if module is initialised
  !!@returns module's initialisation status
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure logical function IsModuleInitialised()
    use, intrinsic :: iso_fortran_env
    implicit none
    IsModuleInitialised = isInitialised
  end function IsModuleInitialised

  !>@brief Initialises list of processes and points which send to this process
  !!@returns list of processes and points which send to this process
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine InitSendRecvLists(Neibs,HaloProcs,NeibPoints,SendList,RecvList)
    use, intrinsic :: iso_fortran_env
    use lattice, only: GetLatticeIndex_M, GetMemorySize, GetProc_G
    use mpiinterface, only: NumProcs, ThisProc, MPISTOP
    use arrayoperations, only: RemoveDuplicates, Sort
    use mpi
    implicit none
    !> Number of neighbours
    integer(intmpi),              intent(out) :: Neibs
    !> List of MPI-ranks of neighbours
    integer(intmpi), allocatable, intent(out) :: HaloProcs(:)
    !> Neighbouring points for each neighbour
    integer(intmpi), allocatable, intent(out) :: NeibPoints(:)
    !> List of points and from this process recieving processes
    integer(int64),    allocatable, intent(out) :: SendList(:,:)
    !> List of points and to this process sending processes
    integer(int64),    allocatable, intent(out) :: RecvList(:,:)

    integer(int64), allocatable :: PointsPerProc_includingThisProc(:)
    integer(intmpi), allocatable :: HaloProcs_includingThisProc(:)
    integer(int64) :: LatticeIndex, MemoryIndex, neibpoint
    integer(intmpi) :: proc, neib
    integer(int64) :: MaxHaloPoints

    ! MPI
    integer(intmpi) :: dest, src, sendtag, recvtag, buffersize, status(mpi_status_size), mpierr
    
    ! 1. Compute halo points and the associated process numbers
    !call GetLocalLatticeIndices_includingHalo_Allocatable(LocalLatticeIndices)

    !allocate(HaloProcs_includingThisProc(size(LocalLatticeIndices)))
    !HaloProcs_includingThisProc = GetProc(LocalLatticeIndices)
    allocate(HaloProcs_includingThisProc(GetMemorySize()))
    do MemoryIndex=1,GetMemorySize()
       HaloProcs_includingThisProc(MemoryIndex) = GetProc_G(GetLatticeIndex_M(MemoryIndex))
    end do
    
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
          do MemoryIndex=1,GetMemorySize()
             LatticeIndex = GetLatticeIndex_M(MemoryIndex)
             if(GetProc_G(LatticeIndex)==HaloProcs(neib)) then
                neibpoint = neibpoint + 1_int64
                RecvList(neibpoint,neib) = LatticeIndex
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

  !>@brief Communication of 1D-array with real entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_real_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetRealSendType
    implicit none
    real(fp), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(fp), asynchronous, allocatable :: sendbuffer(:,:),recvbuffer(:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(BufferIndex,neib) = data(MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,neib),              & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,neib),              & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(MemoryIndex) = recvbuffer(BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_real_rank1

  !>@brief Communication of 2D-array with real entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_real_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetRealSendType
    implicit none
    real(fp), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(fp), asynchronous, allocatable :: sendbuffer(:,:,:),recvbuffer(:,:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(size(data,1),MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(:,BufferIndex,neib) = data(:,MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,1,neib),            & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,1,neib),            & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(:,MemoryIndex) = recvbuffer(:,BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_real_rank2

  !>@brief Communication of 3D-array with real entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_real_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetRealSendType
    implicit none
    real(fp), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(fp), asynchronous, allocatable :: sendbuffer(:,:,:,:),recvbuffer(:,:,:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(size(data,1),size(data,2),MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(:,:,BufferIndex,neib) = data(:,:,MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,1,1,neib),          & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,1,1,neib),          & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(:,:,MemoryIndex) = recvbuffer(:,:,BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_real_rank3

  !>@brief Communication of 4D-array with real entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_real_rank4(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetRealSendType
    implicit none
    real(fp), intent(inout) :: data(:,:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    real(fp), asynchronous, allocatable :: sendbuffer(:,:,:,:,:),recvbuffer(:,:,:,:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(size(data,1),size(data,2),size(data,3),MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(:,:,:,BufferIndex,neib) = data(:,:,:,MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,1,1,1,neib),        & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,1,1,1,neib),        & ! ... and it's first index
            buffersize,           & ! How many points
            GetRealSendType(),    & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(:,:,:,MemoryIndex) = recvbuffer(:,:,:,BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_real_rank4

  !>@brief Communication of 1D-array with complex entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_complex_rank1(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetComplexSendType
    implicit none
    complex(fp), intent(inout) :: data(:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(fp), asynchronous, allocatable :: sendbuffer(:,:),recvbuffer(:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(BufferIndex,neib) = data(MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,neib),              & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,neib),              & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(MemoryIndex) = recvbuffer(BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_complex_rank1

  !>@brief Communication of 2D-array with complex entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_complex_rank2(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetComplexSendType
    implicit none
    complex(fp), intent(inout) :: data(:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(fp), asynchronous, allocatable :: sendbuffer(:,:,:),recvbuffer(:,:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(size(data,1),MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(:,BufferIndex,neib) = data(:,MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,1,neib),            & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,1,neib),            & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(:,MemoryIndex) = recvbuffer(:,BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_complex_rank2

  !>@brief Communication of 3D-array with complex entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_complex_rank3(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetComplexSendType
    implicit none
    complex(fp), intent(inout) :: data(:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(fp), asynchronous, allocatable :: sendbuffer(:,:,:,:),recvbuffer(:,:,:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(size(data,1),size(data,2),MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(:,:,BufferIndex,neib) = data(:,:,MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,1,1,neib),          & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,1,1,neib),          & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(:,:,MemoryIndex) = recvbuffer(:,:,BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_complex_rank3

  !>@brief Communication of 4D-array with complex entries
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_complex_rank4(data)
    use, intrinsic :: iso_fortran_env
    use mpi
    use lattice, only: GetMemoryIndex
    use mpiinterface, only: ThisProc,GetComplexSendType
    implicit none
    complex(fp), intent(inout) :: data(:,:,:,:)

    integer(int8), parameter :: r=rank(data)
    
    ! MPI
    complex(fp), asynchronous, allocatable :: sendbuffer(:,:,:,:,:),recvbuffer(:,:,:,:,:)
    integer(intmpi) :: dest, src, tag, buffersize, mpierr, status(mpi_status_size)

    integer(intmpi) :: neib
    integer(intmpi) :: ValuesPerPoint, MaxHaloPoints
    integer(intmpi) :: BufferIndex
    integer(intmpi), allocatable :: sendrequests(:)
    integer(int64) :: MemoryIndex
    integer(int64) :: LatticeIndex
    
    MaxHaloPoints = MaxVal(NeibPoints,DIM=1)
    
    allocate(sendbuffer(size(data,1),size(data,2),size(data,3),MaxHaloPoints,neibs))
    allocate(recvbuffer,source=sendbuffer)

    ValuesPerPoint = size(data)/size(data,r)
    allocate(sendrequests(neibs))

    packing: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = SendList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(:,:,:,BufferIndex,neib) = data(:,:,:,MemoryIndex)
       end do
    end do packing

    sending: do neib=1,neibs
       dest = HaloProcs(neib)
       tag  = ThisProc()
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_ISEND(&
            sendbuffer(           & ! What to communicate...
            1,1,1,1,neib),        & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            dest, tag,            & ! Destination and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            sendrequests(neib),   & ! Request handle
            mpierr)                 ! Error-code
    end do sending
    
    recieving: do neib=1,neibs
       src = HaloProcs(neib)
       tag = src
       buffersize = ValuesPerPoint*NeibPoints(neib)
       call MPI_RECV(&
            recvbuffer(           & ! What to communicate...
            1,1,1,1,neib),        & ! ... and it's first index
            buffersize,           & ! How many points
            GetComplexSendType(), & ! What type
            src,  tag,            & ! Source and tag
            MPI_COMM_WORLD,       & ! MPI-Communicator
            status,               & ! Status
            mpierr)                 ! Error code
    end do recieving

    unpacking: do neib=1,neibs
       do BufferIndex=1,NeibPoints(neib)
          ! Get global lattice index
          LatticeIndex = RecvList(BufferIndex,neib)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign buffer to data
          data(:,:,:,MemoryIndex) = recvbuffer(:,:,:,BufferIndex,neib)
       end do
    end do unpacking
    
    ! Wait for all sends to finish
    call MPI_WaitAll(neibs,sendrequests,MPI_STATUSES_IGNORE,mpierr)
  end subroutine CommunicateBoundary_complex_rank4
end module halocomm
