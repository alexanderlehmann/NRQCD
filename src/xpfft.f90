!------------------------------------------------------------------------------
! Fast Fourier Transforms between real and momentum space
!------------------------------------------------------------------------------
!
! MODULE: fft
!>@brief
!! Providing interfaces for fast fourier transforms between real and momentum space
!!@author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date
!! 19.02.2019
!!@version
!! 1.0
  !------------------------------------------------------------------------------
module xpfft
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  use mkl_cdft
  use lattice, only: ndim
  use mpiinterface, only: intmpi
  use mpi

  implicit none

  PRIVATE

  public :: &
       InitModule, FinalizeModule,&
       x2p, p2x
       
  !> Module name
  character(len=5), parameter, public ::  modulename='xpfft'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.

  ! Variables  for (x<->p)-FFT
  !> MKL communicator for (x<->p)-FFT
  integer(intmpi) :: mkl_comm
  !>@brief Colour of this process for (x<->p)-FFT.
  !!@details
  !! 0: This process does not participate in FFT\n
  !! 1: This process does participate in FFT
  integer(intmpi), allocatable :: mkl_colors(:)
  !> Pointer to distributed array for (x<->p)-FFT
  type(dfti_descriptor_dm), pointer :: desc
  !> In-place (x<->p)-FFT field
  complex(real64), allocatable :: xp_data(:)
  !> Number of dimensions
  integer(int64), parameter :: ranks = int(nDim,intmpi)

  !> Number of MKL/FFT-processes
  integer(intmpi) :: mklprocs = -1
  !> Size of the local xp-lattice ✓
  integer(int64) :: xp_datasize
  !> Offset in distributed indexing in last dimension
  integer(int64) :: xp_offset
  !> List of lattice indices corresponding to the data in xp_data ✓
  integer(int64), allocatable :: xp_MKL_Indices(:)
  !> List of MPI-ranks, this MPI_COMM_WORLD process has to communicate with
  integer(intmpi), allocatable :: xp_MPIWorld_Procs(:)
  !> List of MPI-ranks, this MKL-process has to communicate with
  integer(intmpi), allocatable :: xp_MKL_Procs(:)
  !> Number of points, this MPI_COMM_WORLD process has to communicate
  integer(intmpi), allocatable :: xp_MPIWorld_CommPoints(:)
  !> Number of points, this MKL-process has to communicate
  integer(intmpi), allocatable :: xp_MKL_CommPoints(:)
  !> Lattice indices which are to be communicated to and from the MKL-processes
  integer(int64), allocatable :: xp_MPIWorld_SendRecvList(:,:)
  !> Lattice indices which are to be communicated this MKL-process has to send and recv
  integer(int64), allocatable :: xp_MKL_SendRecvList(:,:)
  !> Lower lattice boundaries
  integer(int64), allocatable :: xp_LowerLatticeBoundaries(:,:)
  !> Upper lattice boundaries
  integer(int64), allocatable :: xp_UpperLatticeBoundaries(:,:)
contains

  !>@brief Initialises module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use lattice, only: nDim, GetLatticeExtension,&
         GetProc_G, GetProc_fromGeneralIndex, InitLatticeIndices,&
         GetLatticeSpacing, GetVolume, GetMemorySize, GetLocalLatticeSize, GetLatticeIndex_M
    use mpiinterface, only: ThisProc, NumProcs, mpistop
    use mpi
    use arrayoperations, only: RemoveDuplicates, Sort
    implicit none

    ! MPI
    integer(intmpi) :: mpierr, buffersize
    
    ! MKL
    integer(int64) :: lengths(ndim), status, extension

    integer(intmpi) :: proc,xp_MPIWORLD_ncomms,xp_MKL_ncomms
    
    character(len=100) :: errormessage

    integer(int64),  allocatable :: firstindices(:), extensions(:)
    integer(int64) :: LatticeExtensions(ndim), MKLMemoryIndex, MemoryIndex, commpoint, LatticeIndex,is

    integer(int64), allocatable :: indices_includingThisProc(:)
    integer(int64), allocatable :: PointsPerProc_includingThisProc(:)
    integer(int64) :: MaxCommPoints
    integer(int64) :: localfftsize
    
    if(isInitialised) then
       errormessage = 'Error in init of '//modulename//': already initialised.'
       call MPISTOP(errormessage)
    else
       call CheckObligatoryInitialisations

       ! Assigning a colour to each process
       mklprocs = minval([GetLatticeExtension([nDim-1_int8:nDim]),int(NumProcs(),int64)])
       if(allocated(mkl_colors)) deallocate(mkl_colors)
       allocate(mkl_colors(NumProcs()))
       do concurrent (proc=1:NumProcs())
          if(proc .le. mklprocs) then
             mkl_colors(proc) = 1
          else
             mkl_colors(proc) = 0
          end if
       end do
       ! Split communicator
       call mpi_comm_split(MPI_COMM_WORLD,GetMKLColor(),ThisProc(),mkl_comm,mpierr)
       
       ! 1. Partitioning (done by mkl and then communicated)
       if(GetMKLColor()==1) then
          lengths = int(GetLatticeExtension([1_int8:ndim]))

          status = DftiCreateDescriptorDM(int(mkl_comm,int64),desc,&
               DFTI_DOUBLE, DFTI_COMPLEX, ranks, lengths)

          status = DftiGetValueDM(desc,CDFT_LOCAL_NX,extension)
          status = DftiGetValueDM(desc,CDFT_LOCAL_X_START,xp_offset)
          status = DftiSetValueDM(desc,DFTI_FORWARD_SCALE,&
               real(product(GetLatticeSpacing([1_int8:ndim])),real64))
          status = DftiSetValueDM(desc,DFTI_BACKWARD_SCALE,real(1/GetVolume(),real64))
          status = DftiCommitDescriptorDM(desc)

          status = DftiGetValueDM(desc,CDFT_LOCAL_SIZE,localfftsize)
          if(allocated(xp_data)) deallocate(xp_data)
          allocate(xp_data(localfftsize))

          call mpi_barrier(mkl_comm,mpierr)
       end if

       ! Communication of lattice boundaries
       if(allocated(xp_LowerLatticeBoundaries)) deallocate(xp_LowerLatticeBoundaries)
       allocate(xp_LowerLatticeBoundaries(ndim,mklprocs))
       if(allocated(xp_UpperLatticeBoundaries)) deallocate(xp_UpperLatticeBoundaries)
       allocate(xp_UpperLatticeBoundaries(ndim,mklprocs))

       ! 0th MKL-process collecting the boundaries from the other MKL-processes
       if(GetMKLColor()==1) then
          if(ThisProc()==0) then
             allocate(Extensions(mklprocs))  !buffer
             allocate(Firstindices(mklprocs))!buffer
          end if
          
          call mpi_gather(&
            xp_offset,&         ! What to send
            1_intmpi,&          ! How many points
            MPI_INT64_T,&       ! What type to send
            firstindices(&      ! Where to gather ...
            1),&                ! ... and it's first index
            1_intmpi,&          ! How many points
            MPI_INT64_T,&       ! What type to recieve
            0_intmpi,&          ! Gathering process
            mkl_comm,&          ! Communicator
            mpierr)             ! Error-code

          call mpi_gather(&
            Extension,&         ! What to send
            1_intmpi,&          ! How many points
            MPI_INT64_T,&       ! What type to send
            extensions(&        ! Where to gather ...
            1),&                ! ... and it's first index
            1_intmpi,&          ! How many points
            MPI_INT64_T,&       ! What type to recieve
            0_intmpi,&          ! Gathering process
            mkl_comm,&          ! Communicator
            mpierr)             ! Error-code
          
          if(ThisProc()==0) then
             do proc=1,mklprocs
                xp_LowerLatticeBoundaries(1:ndim-1,proc) = 1
                xp_UpperLatticeBoundaries(1:ndim-1,proc) = GetLatticeExtension([1_int8:ndim-1_int8])
                xp_LowerLatticeBoundaries(ndim,proc) = firstindices(proc)
                xp_UpperLatticeBoundaries(ndim,proc) = firstindices(proc) + extensions(proc) -1
             end do
             deallocate(firstindices,extensions) !buffers not needed anymore
          end if
       end if
       
       ! Let 0th (MKL- as well as WORLD)-process brodcast the MKL-lattice-partitions
       buffersize = int(size(xp_LowerLatticeBoundaries),intmpi)
       call mpi_bcast(&
            xp_LowerLatticeBoundaries ,& ! What to broadcast
            buffersize                ,& ! Number of entries in buffer
            MPI_INT64_T               ,& ! Data type of buffer
            0_intmpi                  ,& ! Broadcasting process (root)
            MPI_COMM_WORLD            ,& ! Communicator
            mpierr)                      ! Error-code
       
       call mpi_bcast(&
            xp_UpperLatticeBoundaries ,& ! What to broadcast
            buffersize                ,& ! Number of entries in buffer
            MPI_INT64_T               ,& ! Data type of buffer
            0_intmpi                  ,& ! Broadcasting process (root)
            MPI_COMM_WORLD            ,& ! Communicator
            mpierr)                      ! Error-code
       
       ! 2. Initialise communication lists
       ! i.   Compute the indices of the local lattice indices
       ! ii.  Compute where on which world (a) and on which mkl (b) process this index lies
       if(allocated(xp_MPIWorld_Procs)) deallocate(xp_MPIWorld_Procs)
       allocate(xp_MPIWorld_Procs(GetLocalLatticeSize()))
       LatticeExtensions = GetLatticeExtension([1_int8:ndim])

       !forall(LocalIndex=1:size(LocalLatticeIndices))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex_M(MemoryIndex)
          if(ThisProc()==GetProc_G(LatticeIndex)) then
             is = is+1
             xp_MPIWorld_Procs(is) = GetProc_fromGeneralIndex(&
                  LatticeIndex,&
                  xp_LowerLatticeBoundaries,&
                  xp_UpperLatticeBoundaries,&
                  LatticeExtensions)
          end if
       end do
       !end forall

       ! Initialising send-recv-list of WORLD-processes
       call RemoveDuplicates(xp_MPIWorld_Procs,PointsPerProc_includingThisProc)
       call Sort(xp_MPIWorld_Procs)
       MaxCommPoints = MaxVal(&
            PointsPerProc_includingThisProc, DIM=1)
       xp_MPIWORLD_ncomms = size(xp_MPIWORLD_Procs)
       
       if(allocated(xp_MPIWorld_CommPoints)) deallocate(xp_MPIWorld_CommPoints)
       allocate(xp_MPIWorld_CommPoints(xp_MPIWORLD_ncomms))
       
       if(allocated(xp_MPIWorld_SendRecvList)) deallocate(xp_MPIWorld_SendRecvList)
       allocate(xp_MPIWorld_SendRecvList(MaxCommPoints,xp_MPIWorld_ncomms))

       do proc=1,size(xp_MPIWorld_Procs)
          xp_MPIWorld_Commpoints(proc) = PointsPerProc_includingThisProc(proc)

          commpoint=0
          is=0
          do MemoryIndex=1,GetMemorySize()
             LatticeIndex = GetLatticeIndex_M(MemoryIndex)
             if(GetProc_G(LatticeIndex)==ThisProc()) then
                is = is + 1
                if(GetProc_fromGeneralIndex(&
                     LatticeIndex,&
                     xp_LowerLatticeBoundaries,&
                     xp_UpperLatticeBoundaries,&
                     LatticeExtensions)&
                     ==xp_MPIWorld_Procs(proc)) then
                   commpoint = commpoint + 1
                   xp_MPIWorld_SendRecvList(commpoint,proc) = LatticeIndex
                end if
             end if
          end do
       end do
       ! Now the World-processes know to which MKL-processes they will send data
       !-------------------------------------------------------------------------
       ! Now letting all MKL-processes know from whom they recieve what
       if(GetMKLColor()==1) then
          xp_datasize = product(&
               1 &
               + xp_UpperLatticeBoundaries(:,ThisProc()+1) &
               - xp_LowerLatticeBoundaries(:,ThisProc()+1) )

          if(allocated(xp_MKL_indices)) deallocate(xp_MKL_indices)
          allocate(xp_MKL_indices(xp_datasize))

          call InitLatticeIndices(xp_MKL_indices,&
               xp_LowerLatticeBoundaries(:,ThisProc()+1),&
               xp_UpperLatticeBoundaries(:,ThisProc()+1))

          if(allocated(xp_MKL_Procs)) deallocate(xp_MKL_Procs)
          allocate(xp_MKL_Procs(size(xp_MKL_indices)))
          xp_MKL_Procs = GetProc_G(xp_MKL_indices)

          call RemoveDuplicates(xp_MKL_Procs,PointsPerProc_includingThisProc)
          call Sort(xp_MKL_Procs)

          xp_MKL_ncomms = size(xp_MKL_Procs)
          
          if(allocated(xp_MKL_CommPoints)) deallocate(xp_MKL_CommPoints)
          allocate(xp_MKL_CommPoints(xp_MKL_ncomms))
          
          MaxCommPoints = MaxVal(&
               PointsPerProc_includingThisProc,DIM=1) ! Where to look

          if(allocated(xp_MKL_SendRecvList)) deallocate(xp_MKL_SendRecvList)
          allocate(xp_MKL_SendRecvList(MaxCommPoints,xp_MKL_ncomms))

          do proc=1,size(xp_MKL_Procs)
             xp_MKL_CommPoints(proc) = PointsPerProc_includingThisProc(proc)

             commpoint=0
             do MKLMemoryIndex=1,size(xp_MKL_indices)
                if(GetProc_G(xp_MKL_indices(MKLMemoryIndex))==xp_MKL_Procs(proc)) then
                   commpoint = commpoint + 1_int64
                   xp_MKL_SendRecvList(commpoint,proc) = xp_MKL_indices(MKLMemoryIndex)
                end if
             end do
          end do
          ! Now the MKL-processes know from whom they will get the data
          !-------------------------------------------------------------
       end if ! is MKL-process
       
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule

  !>@brief Returns MPI-colour regarding MKL's distributed cluster FFT
  !!@returns MPI-colour regarding MKL's distributed cluster FFT
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
  pure integer(intmpi) function GetMKLColor(proc)
    use mpiinterface, only: intmpi, ThisProc
    implicit none
    integer(intmpi), intent(in), optional :: proc
    if(present(proc)) then
       GetMKLColor = mkl_colors(proc+1_intmpi)
    else
       GetMKLColor = mkl_colors(ThisProc()+1_intmpi)
    end if
  end function GetMKLColor

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
    character(len=70) :: errormessage
    
    if(.not.IsLatticeInitialised()) then
       errormessage = 'Error in init of '//modulename//': '//latticename//' is not initialised.'
       call mpistop(errormessage)
    end if
  end subroutine CheckObligatoryInitialisations
  
  !>@brief Finalizes module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.02.2019
  !!@version 1.0
  impure subroutine FinalizeModule
    use mpiinterface, only: mpistop,thisproc
    implicit none
    character(len=70) :: errormessage
    
    integer(int64) :: status
    
    if(isInitialised) then
       FFTprocess: if(mkl_colors(ThisProc()+1)==1) then
          deallocate(&
               xp_data&
               ,xp_MKL_Indices&
               ,xp_MKL_Procs&
               ,xp_MKL_CommPoints&
               ,xp_MKL_SendRecvList&
               )
          status = DftiFreeDescriptorDM(desc)
       end if FFTprocess
       
       deallocate(&
            mkl_colors&
            ,xp_MPIWorld_Procs&
            ,xp_MPIWorld_CommPoints&
            ,xp_MPIWorld_SendRecvList&
            ,xp_LowerLatticeBoundaries&
            ,xp_UpperLatticeBoundaries&
            )
       
       IsInitialised = .FALSE.
    else
       errormessage = 'Error in finalization of '//modulename//': is not initialised.'
       call MPISTOP(errormessage)
    end if
  end subroutine FinalizeModule
  
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

  !>@brief Loads data from the original processes into the distributed FFT-array
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
  impure subroutine LoadDataToFFT(data)
    use precision, only: fp
    use lattice
    use mpi
    use mpiinterface, only: ThisProc, intmpi, SyncAll
    implicit none
    !> Fourier-transform input
    complex(fp), intent(in) :: data(:)

    complex(real64), allocatable, asynchronous :: sendbuffer(:,:), recvbuffer(:,:)
    integer(intmpi) :: dest, src, sendtag, recvtag, buffersize, mpierr
    integer(intmpi) :: recvstatus(MPI_STATUS_SIZE)
    integer(intmpi), allocatable :: sendrequest(:),sendstatus(:,:)
    integer(intmpi) :: proc
    integer(intmpi) :: MaxPoints
    integer(intmpi) :: BufferIndex
    integer(int64) :: MKLMemoryIndex,MemoryIndex, LatticeIndex

    MaxPoints = MaxVal(xp_MPIWorld_Commpoints,DIM=1)

    allocate(sendbuffer(MaxPoints,size(xp_MPIWorld_Procs)))
    allocate(sendrequest(size(xp_MPIWorld_Procs)))

    Destinations: do proc=1,size(xp_MPIWorld_Procs)
       buffersize = xp_MPIWorld_Commpoints(proc)

       ! Prepare buffer with send-data for FFT-process
       packing: do BufferIndex=1,xp_MPIWorld_Commpoints(proc)
          ! Get global lattice index
          LatticeIndex = xp_MPIWorld_SendRecvList(BufferIndex,proc)

          ! Get memory index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          sendbuffer(BufferIndex,proc) = data(MemoryIndex)
       end do packing

       ! Sending the data to the FFT-process
       dest    = xp_MPIWorld_Procs(proc)
       sendtag = ThisProc() + mklprocs*xp_MPIWorld_Procs(proc)
       
       call MPI_ISend(&       
            sendbuffer(&     ! What to send...
            1,proc),&        ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to send
            dest, sendtag,&  ! Destination and sendtag
            mpi_comm_world,& ! Communicator
            sendrequest(proc),& ! Request handle
            mpierr)          ! Error-code
    end do Destinations   

    if(GetMKLColor()==1) then
       MaxPoints = MaxVal(xp_MKL_Commpoints,DIM=1)
       allocate(recvbuffer(MaxPoints,size(xp_MKL_Procs)))
       
       Sources: do proc=1,size(xp_MKL_Procs)
          buffersize = xp_MKL_Commpoints(proc)

          ! Recieving the data
          src     = xp_MKL_Procs(proc)
          recvtag = xp_MKL_Procs(proc) + mklprocs*ThisProc()
          
          call MPI_Recv(&       
               recvbuffer(&     ! Where to recieve...
               1,proc),&        ! ... and it's first index
               buffersize,&     ! How many points
               MPI_COMPLEX16,&  ! What type to recieve
               src, recvtag,&   ! Source and recvtag
               mpi_comm_world,& ! Communicator
               recvstatus, mpierr)  ! Status and error-code

          ! Unpacking buffer and saving into distributed FFT-array
          unpacking: do BufferIndex=1,xp_MKL_Commpoints(proc)
             ! Get global lattice index
             LatticeIndex = xp_MKL_SendRecvList(BufferIndex,proc)

             ! Get local index in data
             MKLMemoryIndex = GetMKLMemoryIndex(LatticeIndex)

             ! Assign buffer to data
             xp_data(MKLMemoryIndex) = Recvbuffer(BufferIndex,proc)
          end do unpacking
       end do Sources

       deallocate(RecvBuffer)
    end if

    ! Checking all status
    !allocate(sendstatus(MPI_STATUS_SIZE,size(sendrequest)))
    !call MPI_WAITALL(int(size(xp_MPIWorld_Procs),intmpi),sendrequest,sendstatus,mpierr)
    !deallocate(sendrequest,sendstatus,SendBuffer)
    !call SyncAll
    deallocate(sendrequest,SendBuffer)
  end subroutine LoadDataToFFT

  !>@brief Distributes data from distributed FFT-array back to the original processes
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
  impure subroutine DistributeDataFromFFT(data)
    use precision, only: fp
    use lattice
    use mpi
    use mpiinterface, only: ThisProc
    implicit none
    !> Fourier-transform output
    complex(fp), intent(out) :: data(:)

    complex(real64), allocatable, asynchronous :: sendbuffer(:,:), recvbuffer(:,:)
    integer(intmpi) :: dest, src, sendtag, recvtag, buffersize, mpierr
    integer(intmpi) :: recvstatus(MPI_STATUS_SIZE)
    integer(intmpi), allocatable :: sendrequest(:),sendstatus(:,:)
    integer(intmpi) :: proc
    integer(intmpi) :: MaxPoints
    integer(intmpi) :: BufferIndex
    integer(int64) :: MKLMemoryIndex, MemoryIndex, LatticeIndex

    if(GetMKLColor()==1) then
       MaxPoints = MaxVal(xp_MKL_Commpoints,DIM=1)
       allocate(sendbuffer(MaxPoints,size(xp_MKL_Procs)))
       allocate(sendrequest(size(xp_MKL_Procs)))
       
       Destinations: do proc=1,size(xp_MKL_Procs)
          buffersize = xp_MKL_Commpoints(proc)

          ! Packing buffer
          packing: do BufferIndex=1,xp_MKL_Commpoints(proc)
             ! Get global lattice index
             LatticeIndex = xp_MKL_SendRecvList(BufferIndex,proc)

             ! Get local index in data
             MKLMemoryIndex = GetMKLMemoryIndex(LatticeIndex)

             ! Assign data to buffer
             SendBuffer(BufferIndex,proc) = xp_data(MKLMemoryIndex)
          end do packing
          
          ! Sending the data
          dest    = xp_MKL_Procs(proc)
          sendtag = xp_MKL_Procs(proc) + mklprocs*ThisProc()
          call MPI_ISend(&       
               sendbuffer(&     ! Where to recieve...
               1,proc),&        ! ... and it's first index
               buffersize,&     ! How many points
               MPI_COMPLEX16,&  ! What type to send
               dest, sendtag,&  ! Destination and sendtag
               mpi_comm_world,& ! Communicator
               sendrequest(proc), mpierr)  ! Status and error-code
       end do Destinations
    end if

    MaxPoints = MaxVal(xp_MPIWorld_Commpoints,DIM=1)
    allocate(RecvBuffer(MaxPoints,size(xp_MPIWorld_Procs)))
    Sources: do proc=1,size(xp_MPIWorld_Procs)
       buffersize = xp_MPIWorld_Commpoints(proc)

       ! Sending the data to the FFT-process
       src     = xp_MPIWorld_Procs(proc)
       recvtag = ThisProc() + mklprocs*xp_MPIWorld_Procs(proc)
       call MPI_Recv(&       
            recvbuffer(&     ! What to recieve...
            1,proc),&        ! ... and it's first index
            buffersize,&     ! How many points
            MPI_COMPLEX16,&  ! What type to recieve
            src, recvtag,&   ! Source and recvtag
            mpi_comm_world,& ! Communicator
            recvstatus,&     ! Status
            mpierr)          ! Error-code

       ! Unpack buffer into FFT-output
       unpacking: do BufferIndex=1,xp_MPIWorld_Commpoints(proc)
          ! Get global lattice index
          LatticeIndex = xp_MPIWorld_SendRecvList(BufferIndex,proc)

          ! Get local index in data
          MemoryIndex = GetMemoryIndex(LatticeIndex)

          ! Assign data to buffer
          data(MemoryIndex) = recvbuffer(BufferIndex,proc)
       end do unpacking
    end do Sources

    ! Checking status
    if(GetMKLColor()==1) then
       allocate(sendstatus(MPI_STATUS_SIZE,size(sendrequest)))
       call MPI_WAITALL(int(size(sendrequest),intmpi),sendrequest,sendstatus,mpierr)
       deallocate(sendbuffer,sendstatus,sendrequest)
    end if
    deallocate(recvbuffer)
  end subroutine DistributeDataFromFFT
  
  !>@brief DFT for ndim-dimensional signal from real space to momentum space using MKL-CDFT
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
  impure subroutine x2p(data)
    use precision, only: fp
    implicit none
    !> Fourier-transform input and output
    complex(fp), intent(inout) :: data(:)

    integer(int64) :: status

    call LoadDataToFFT(data)
    
    ! Perform fourier transform
    status = DftiComputeForwardDM(desc,xp_data)

    call DistributeDataFromFFT(data)
  end subroutine x2p
  
  !>@brief DFT for ndim-dimensional signal from momentum space to real space using MKL-CDFT
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
  impure subroutine p2x(data)
    use precision, only: fp

    use mpiinterface, only: ThisProc
    implicit none
    !> Fourier-transform input and output
    complex(fp), intent(inout) :: data(:)

    integer(int64) :: status

    call LoadDataToFFT(data)
    
    ! Perform fourier transform
    status = DftiComputeBackwardDM(desc,xp_data)

    call DistributeDataFromFFT(data)
    
  end subroutine p2x

  !>@brief Returns local FFT-index based on global lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 21.02.2019
  !!@version 1.0
  pure elemental integer(int64) function GetMKLMemoryIndex(LatticeIndex)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: ThisProc
    implicit none
    integer(int64), intent(in) :: LatticeIndex

    GetMKLMemoryIndex = FindLoc(&
                                ! Where to look
         xp_MKL_Indices,dim=1,&
                                ! What to look for
         value = LatticeIndex) &
         ! Offset
         + (xp_offset-xp_LowerLatticeBoundaries(ndim,ThisProc()+1))&
         *product(&
         1&
         +xp_UpperLatticeBoundaries(1:ndim-1,ThisProc()+1)&
         -xp_LowerLatticeBoundaries(1:ndim-1,ThisProc()+1))
  end function GetMKLMemoryIndex
end module xpfft
