!> @brief
!! Providing numerical tolerances
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 17.01.2019
!! @version
!! 1.0
module tolerances
  use, intrinsic :: iso_fortran_env
  implicit none

  PUBLIC

  real(real64), parameter :: tol_zero_real64=1.E-12

end module tolerances
