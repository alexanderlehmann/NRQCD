!------------------------------------------------------------------------------
! MATHCONSTANTS, Mathematial constants
!------------------------------------------------------------------------------
!
! MODULE: mathcontants
!> @brief Mathematical constants like \f$\pi, \text{i}\f$ etc.
!! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date 03.09.2018
!! @version 1.0
! REVISION HISTORY:
! 03 09 2018
module mathconstants
  use, intrinsic :: iso_fortran_env
  implicit none

  PUBLIC
  
  !> \f$\pi\f$
  real(real64), parameter :: pi=acos(-1._real64)
end module mathconstants
