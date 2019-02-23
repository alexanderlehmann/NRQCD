!----------------------------------------------------------------------
! Precision module
!----------------------------------------------------------------------
!
! MODULE: precision
!>@brief Defines numerical precision
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 19.02.2019
!!@version 1.0
! REVISION HISTORY:
! 19 02 2019 - Initial version
!----------------------------------------------------------------------
module precision
  use, intrinsic :: iso_fortran_env
  implicit none

  PUBLIC
  !> Floating number precision
  integer, parameter :: fp = real64
end module precision
