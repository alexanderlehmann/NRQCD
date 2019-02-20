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
       InitModule, FinalizeModule
       
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

  !> List of lattice indices corresponding to the data in xp_data
  integer(int64), allocatable :: LocalLatticeIndices(:)
  !> List of MPI-ranks to communicate with in (x<->p)-FFT
  integer(intmpi), allocatable :: CommProcs(:)
  !> Lower lattice boundaries
  integer(int64), allocatable :: xp_LowerLatticeBoundaries(:,:)
  !> Upper lattice boundaries
  integer(int64), allocatable :: xp_UpperLatticeBoundaries(:,:)
contains

  !> @brief Initialises module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use lattice, only: nDim, GetLatticeExtension
    use mpiinterface, only: ThisProc, NumProcs, mpistop
    use mpi
    implicit none

    ! MPI
    integer(intmpi) :: mpierr, buffersize
    integer(intmpi) :: maxmklprocs
    
    ! MKL
    integer(int64) :: lengths(ndim), status, extension, firstindex

    integer(intmpi) :: proc
    
    character(len=100) :: errormessage

    integer(int64), allocatable :: firstindices(:), extensions(:)

    if(isInitialised) then
       errormessage = 'Error in init of '//modulename//': already initialised.'
       call MPISTOP(errormessage)
    else
       call CheckObligatoryInitialisations

       ! Assigning a colour to each process
       maxmklprocs = minval([GetLatticeExtension(nDim),int(NumProcs(),int64)])
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
            firstindices(&      ! What to recieve ...
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
            extensions(&        ! What to recieve ...
            1),&                ! ... and it's first index
            1_intmpi,&          ! How many points
            MPI_INT64_T,&       ! What type to recieve
            0_intmpi,&          ! Gathering process
            mkl_comm,&          ! Communicator
            mpierr)             ! Error-code

          if(ThisProc()==0) then
             forall(proc=1:maxmklprocs)
                xp_LowerLatticeBoundaries(1:ndim-1,proc) = 1
                xp_UpperLatticeBoundaries(1:ndim-1,proc) = GetLatticeExtension([1_int8:ndim-1_int8])
                xp_LowerLatticeBoundaries(ndim,proc) = firstindices(proc)
                xp_UpperLatticeBoundaries(ndim,proc) = firstindices(proc) + extensions(proc) -1
             end forall
             deallocate(firstindices,extensions) !buffers not needed anymore
          end if
       end if
       ! Let 0th (MKL- as well as WORLD)-process brodcast the LatticeExtensions
       buffersize = int(size(xp_LowerLatticeBoundaries),intmpi)
       call mpi_bcast(&
            xp_LowerLatticeBoundaries ,& ! What to send
            buffersize                ,& ! Number of entries in buffer
            MPI_INT64_T               ,& ! Data type of buffer
            0_intmpi                  ,& ! Broadcasting process (root)
            MPI_COMM_WORLD            ,& ! Communicator
            mpierr)                      ! Error-code
       
       call mpi_bcast(&
            xp_UpperLatticeBoundaries ,& ! What to send
            buffersize                ,& ! Number of entries in buffer
            MPI_INT64_T               ,& ! Data type of buffer
            0_intmpi                  ,& ! Broadcasting process (root)
            MPI_COMM_WORLD            ,& ! Communicator
            mpierr)                      ! Error-code
       
       !do proc=0,NumProcs()-1
       !  print*,'Process',proc
       !   do status=1,maxmklprocs
       !      print*,'MKL-process:',status-1
       !      print*,xp_LowerLatticeBoundaries(:,status)
       !      print*,xp_UpperLatticeBoundaries(:,status)
       !   end do
       !   call flush(6)
       !   call mpi_barrier(mpi_comm_world,mpierr)
       !end do

       
       ! 2. Initialise list where the which lattice points have sent
       ! 3. Initialise list from where lattice points are recieved
       
       
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule

  !> @brief Returns MPI-colour regarding MKL's distributed cluster FFT
  !! @returns MPI-colour regarding MKL's distributed cluster FFT
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  pure integer function GetColor()
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
end module xpfft
