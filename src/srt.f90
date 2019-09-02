module srt
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  IMPLICIT NONE

  PRIVATE

  PUBLIC metric
contains

  pure real(fp) function metric(mu,nu)
    IMPLICIT NONE
    integer(int8), intent(in) :: mu,nu
    
    if(mu==nu) then
       if(mu==0) then
          metric = -1
       else
          metric = +1
       end if
    else
       metric = 0
    end if

  end function metric

end module srt
