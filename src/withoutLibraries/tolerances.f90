!>@brief Providing numerical tolerances
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 24.02.2019
!!@version 1.1
!!@todo Adding quad precision zero-tolerance
module tolerances
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  implicit none

  PRIVATE

  public :: GetZeroTol,InitModule, IsModuleInitialised
  
  !> Module name
  character(len=10), parameter, public ::  modulename='tolerances'
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.
  !> Zero-tolerance for single precision
  real(real32), parameter :: tol_zero_real32=1.E-6
  !> Zero-tolerance for double precision
  real(real64), parameter :: tol_zero_real64=1.E-12

  !> Zero-tolerance
  real(fp) :: tol_zero=tol_zero_real64

contains
  !>@brief Initialization of module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  !!@todo Adding quad precision zero-tolerance
  impure subroutine InitModule
    use mpiinterface, only: mpistop
    implicit none

    select case(fp)
    case(real32)
       tol_zero = tol_zero_real32
    case(real64)
       tol_zero = tol_zero_real64
    case default
       call MPISTOP('Error in initialisation of '//modulename&
               //': unsupported floating point precision.')
    end select

    isInitialised = .TRUE.
  end subroutine InitModule
  
  !>@brief Returns, if module is initialised
  !!@returns module's initialisation status
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  pure logical function IsModuleInitialised()
    implicit none
    IsModuleInitialised = IsInitialised
  end function IsModuleInitialised

  !>@brief Zero-tolerance
  !!@returns Zero-tolerance
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  pure real(fp) function GetZeroTol()
    implicit none
    GetZeroTol = tol_zero
  end function GetZeroTol
end module tolerances
