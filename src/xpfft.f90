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
    integer(intmpi) :: mpierr
    integer(intmpi) :: maxmklprocs
    
    ! MKL
    integer(int64) :: lengths(ndim), status, local_extension, local_firstindex
    
    character(len=100) :: errormessage

    if(isInitialised) then
       errormessage = 'Error in init of '//modulename//': already initialised.'
       call MPISTOP(errormessage)
    else
       call CheckObligatoryInitialisations

       ! Assigning a colour to each process
       maxmklprocs = GetLatticeExtension(nDim)
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
          
          !call mpi_barrier(mkl_comm,mpierr)

          status = DftiCreateDescriptorDM(int(mkl_comm,int64),desc,&
               DFTI_DOUBLE, DFTI_COMPLEX, ranks, lengths)

          print*,thisproc(),mpierr
          call flush(6)
          call mpi_barrier(mkl_comm,mpierr)

          status = DftiGetValueDM(desc,CDFT_LOCAL_NX,local_extension)
          status = DftiGetValueDM(desc,CDFT_LOCAL_X_START,local_firstindex)

          !print*,ThisProc(),local_firstindex, status

       else
          !print*,ThisProc()

       end if


       
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
