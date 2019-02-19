!------------------------------------------------------------------------------
! Fast Fourier Transform interfaces
!------------------------------------------------------------------------------
!
! MODULE: fft
!> @brief
!! Providing interfaces for fast fourier transforms
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 10.01.2019
!! @version
!! 1.0
  !------------------------------------------------------------------------------
module fft
  use, intrinsic :: iso_fortran_env
  use mkl_cdft
  use lattice, only: ndim

  implicit none

  PRIVATE

  public :: &
       InitModule, FinalizeModule
       
  !> Module name
  character(len=3), parameter, public ::  modulename='fft'
  
  !> Contains information, if module is initialised
  logical :: IsInitialised = .false.

  ! Variables  for (x<->p)-FFT
  !> MKL communicator for (x<->p)-FFT
  integer :: mkl_comm
  !> Pointer to distributed array for (x<->p)-FFT
  type(dfti_descriptor_dm), pointer :: xp_desc
  !> @brief Colour of this process for (x<->p) -FFT.
  !! @details
  !! 0: This process does not participate in FFT\n
  !! 1: This process does participate in FFT
  integer :: mkl_color=-1
  
contains

  !> @brief Initialises module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    implicit none

    call CheckObligatoryInitialisations

    ! Assigning a colour to each process
    

    
    IsInitialised = .TRUE.
  end subroutine InitModule

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
    implicit none


    IsInitialised = .FALSE.
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
end module fft
