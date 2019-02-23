!>@brief Providing numerical tolerances
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 17.01.2019
!!@version 1.0
module tolerances
  use, intrinsic :: iso_fortran_env
  implicit none

  PUBLIC

  !> Tolerance as zero for single precision
  real(real32), parameter :: tol_zero_real32=1.E-6
  !> Tolerance as zero for double precision
  real(real64), parameter :: tol_zero_real64=1.E-14

end module tolerances
