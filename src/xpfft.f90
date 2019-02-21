!------------------------------------------------------------------------------
! Fast Fourier Transforms between real and momentum space
!------------------------------------------------------------------------------
!
! MODULE: fft
!> @brief
!! Providing interfaces for fast fourier transforms between real and momentum space
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 19.02.2019
!! @version
!! 1.0
  !------------------------------------------------------------------------------
module xpfft
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  use mkl_cdft
  use lattice, only: ndim
  use mpiinterface, only: intmpi

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
  !> @brief Colour of this process for (x<->p)-FFT.
  !! @details
  !! 0: This process does not participate in FFT\n
  !! 1: This process does participate in FFT
  integer(intmpi) :: mkl_color=-1
  !> Pointer to distributed array for (x<->p)-FFT
  type(dfti_descriptor_dm), pointer :: desc
  !> In-place (x<->p)-FFT field
  complex(fp), allocatable :: data(:)
  !> Number of dimensions
  integer(int64), parameter :: ranks = int(nDim,intmpi)

  !> Size of the local xp-lattice ✓
  integer(int64) :: xp_datasize
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

  !> @brief Initialises module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 21.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use lattice, only: nDim, GetLatticeExtension, GetLocalLatticeIndices_allocatable,&
         GetProc,GetProc_fromGeneralIndex, InitLatticeIndices,&
         GetLatticeSpacing, GetVolume
    use mpiinterface, only: ThisProc, NumProcs, mpistop
    use mpi
    use arrayoperations, only: RemoveDuplicates, Sort
    implicit none

    ! MPI
    integer(intmpi) :: mpierr, buffersize
    integer(intmpi) :: maxmklprocs
    
    ! MKL
    integer(int64) :: lengths(ndim), status, extension, firstindex

    integer(intmpi) :: proc,xp_MPIWORLD_ncomms,xp_MKL_ncomms
    
    character(len=100) :: errormessage

    integer(int64),  allocatable :: firstindices(:), extensions(:)
    integer(int64),  allocatable :: locallatticeindices(:)
    integer(int64) :: LatticeExtensions(ndim), LocalIndex, commpoint

    integer(int64), allocatable :: indices_includingThisProc(:)
    integer(int64), allocatable :: PointsPerProc_includingThisProc(:)
    integer(int64) :: MaxCommPoints
    
    if(isInitialised) then
       errormessage = 'Error in init of '//modulename//': already initialised.'
       call MPISTOP(errormessage)
    else
       call CheckObligatoryInitialisations

       ! Assigning a colour to each process
       maxmklprocs = minval([GetLatticeExtension([nDim-1_int8:nDim]),int(NumProcs(),int64)])
       if(ThisProc() .le. maxmklprocs-1) then
          mkl_color = 1
       else
          mkl_color = 0
       end if
       ! Split communicator
       call mpi_comm_split(MPI_COMM_WORLD,mkl_color,ThisProc(),mkl_comm,mpierr)
       
       ! 1. Partitioning (done by mkl and then communicated)
       if(mkl_color==1) then
          lengths = int(GetLatticeExtension([1_int8:ndim]))

          status = DftiCreateDescriptorDM(int(mkl_comm,int64),desc,&
               DFTI_DOUBLE, DFTI_COMPLEX, ranks, lengths)

          status = DftiGetValueDM(desc,CDFT_LOCAL_NX,extension)
          status = DftiGetValueDM(desc,CDFT_LOCAL_X_START,firstindex)
          status = DftiSetValueDM(desc,DFTI_FORWARD_SCALE,product(GetLatticeSpacing([1_int8:ndim])))
          status = DftiSetValueDM(desc,DFTI_BACKWARD_SCALE,1._real64/GetVolume())
          status = DftiCommitDescriptorDM(desc)
       end if

       ! Communication of lattice boundaries
       if(allocated(xp_LowerLatticeBoundaries)) deallocate(xp_LowerLatticeBoundaries)
       allocate(xp_LowerLatticeBoundaries(ndim,maxmklprocs))
       if(allocated(xp_UpperLatticeBoundaries)) deallocate(xp_UpperLatticeBoundaries)
       allocate(xp_UpperLatticeBoundaries(ndim,maxmklprocs))

       ! 0th MKL-process collecting the boundaries from the other MKL-processes
       if(mkl_color==1) then
          if(ThisProc()==0) then
             allocate(Extensions(maxmklprocs))  !buffer
             allocate(Firstindices(maxmklprocs))!buffer
          end if
          
          call mpi_gather(&
            firstindex,&        ! What to send
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
             do proc=1,maxmklprocs
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
       call GetLocalLatticeIndices_allocatable(LocalLatticeIndices)
       if(allocated(xp_MPIWorld_Procs)) deallocate(xp_MPIWorld_Procs)
       allocate(xp_MPIWorld_Procs(size(LocalLatticeIndices)))
       LatticeExtensions = GetLatticeExtension([1_int8:ndim])

       forall(LocalIndex=1:size(LocalLatticeIndices))
          xp_MPIWorld_Procs(LocalIndex) = GetProc_fromGeneralIndex(&
               LocalLatticeIndices(LocalIndex),&
               xp_LowerLatticeBoundaries,&
               xp_UpperLatticeBoundaries,&
               LatticeExtensions)
       end forall

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
          do LocalIndex=1,size(LocalLatticeIndices)
             if(GetProc_fromGeneralIndex(LocalLatticeIndices(LocalIndex),&
                  xp_LowerLatticeBoundaries,&
                  xp_UpperLatticeBoundaries,&
                  LatticeExtensions)&
                  ==xp_MPIWorld_Procs(proc)) then
                commpoint = commpoint + 1_int64
                xp_MPIWorld_SendRecvList(commpoint,proc) = LocalLatticeIndices(LocalIndex)
             end if
          end do
       end do
       ! Now the World-processes know to which MKL-processes they will send data
       !-------------------------------------------------------------------------
       ! Now letting all MKL-processes know from whom they recieve what
       if(mkl_color==1) then
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
          xp_MKL_Procs = GetProc(xp_MKL_indices)

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
             do LocalIndex=1,size(xp_MKL_indices)
                if(GetProc(xp_MKL_indices(LocalIndex))==xp_MKL_Procs(proc)) then
                   commpoint = commpoint + 1_int64
                   xp_MKL_SendRecvList(commpoint,proc) = xp_MKL_indices(LocalIndex)
                end if
             end do
          end do
          ! Now the MKL-processes know from whom they will get the data
          !-------------------------------------------------------------
       end if ! is MKL-process
       
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule

  !> @brief Returns MPI-colour regarding MKL's distributed cluster FFT
  !! @returns MPI-colour regarding MKL's distributed cluster FFT
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  pure integer(intmpi) function GetColor()
    implicit none
    GetColor = mkl_color
  end function GetColor

  !> @brief Checks previous necessary initialisations
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
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
  
  !> @brief Finalizes module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine FinalizeModule
    use mpiinterface, only: mpistop
    implicit none
    character(len=70) :: errormessage

    if(isInitialised) then
       

       IsInitialised = .FALSE.
    else
       errormessage = 'Error in finalization of '//modulename//': is not initialised.'
       call MPISTOP(errormessage)
    end if
  end subroutine FinalizeModule
  
  !>@brief Returns, if module is initialised
  !! @returns module's initialisation status
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  pure logical function IsModuleInitialised()
    implicit none
    IsModuleInitialised = IsInitialised
  end function IsModuleInitialised

  !> @brief DFT for ndim-dimensional signal from real space to momentum space using MKL-CDFT
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 21.02.2019
  !! @version 1.0
  impure subroutine x2p(data_)
    implicit none
    !> Fourier-transform input and output
    complex(real64), intent(inout) :: data_(:)

    integer status
    
    !status = DftiComputeForwardDM(xp_desc,data_)
  end subroutine x2p

  !> @brief DFT for ndim-dimensional signal from momentum space to real space using MKL-CDFT
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 21.02.2019
  !! @version 1.0
  impure subroutine p2x(data_)
    implicit none
    !> Fourier-transform input and output
    complex(real64), intent(inout) :: data_(:)

    integer status

    !status = DftiComputeBackwardDM(xp_desc,data_)
    
  end subroutine p2x
end module xpfft
