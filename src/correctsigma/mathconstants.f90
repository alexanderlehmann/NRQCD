!------------------------------------------------------------------------------
! MATHCONSTANTS, Mathematial constants
!------------------------------------------------------------------------------
!
! MODULE: mathcontants
!>@brief Mathematical constants like \f$\pi, \text{i}\f$ etc.
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 03.09.2018
!!@version 1.0
! REVISION HISTORY:
! 03 09 2018
module mathconstants
  use precision
  implicit none

  PUBLIC
  
  !>\f$\pi\f$
  real(fp), parameter :: pi=acos(-1._fp)
end module mathconstants
