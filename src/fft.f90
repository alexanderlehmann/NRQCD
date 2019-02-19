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

    IsInitialised = .TRUE.
  end subroutine InitModule

  !> @brief Checks previous necessary initialisations
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 19.02.2019
  !! @version 1.0
  impure subroutine CheckObligatoryInitialisations
    use, intrinsic :: iso_fortran_env
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
