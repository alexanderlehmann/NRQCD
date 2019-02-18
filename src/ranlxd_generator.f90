!------------------------------------------------------------------------------
! RANLXD_GENRATOR, Random number generator
!------------------------------------------------------------------------------
!
! MODULE: ranlxd_generator
!> @brief Ranlux-module from Lüscher
!! @details
!! See the notes\n
!! User's guide for ranlxs and ranlxd [F90 programs] (December 1997)\n
!! Double precision implementation of the random number generator ranlux (December 1997)
!! @author
!! Martin Lüscher <luscher@mail.desy.de>\n
!! Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date 15.12.1997 (03.09.2018)
!! @version 2.1a
! REVISION HISTORY:
! 15 12 1997 - Version from Martin Lüscher
! 03 09 2018 - Small adjustments for modern fortran
module ranlxd_generator
  use, intrinsic :: iso_fortran_env
  implicit none
  public :: ranlxd,ranlxd_init,ranlxd_get,ranlxd_reset

  private :: define_constants,error,update

  integer(int64),save,private :: pr,ir,jr,ir_old,init=0
  integer(int64),dimension(0:11),save,private :: next 

  real(real64),save,private :: zero,one,sbase,base,one_bit,carry
  real(real64),dimension(0:11),save,private :: xdbl

contains
  !> @brief Error messages
  !! @author Martin Lüscher <luscher@mail.desy.de>\n
  !! Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 15.12.1997 (03.09.2018)
  !! @version 2.1a
  impure subroutine error(no)
    use iso_fortran_env
    !> Number of the error message
    integer(int64),intent(in) :: no

    select case(no)
    case(0)
       write(ERROR_UNIT,*) "Error in ranlxd_init"
       write(ERROR_UNIT,*) "Arithmetic on this machine is not suitable for ranlxd"
    case(1)
       write(ERROR_UNIT,*) "Error in ranlxd_init"
       write(ERROR_UNIT,*) "Bad choice of luxury level (should be 1 or 2)"
    case(2)
       write(ERROR_UNIT,*) "Error in ranlxd_init"
       write(ERROR_UNIT,*) "Bad choice of seed (should be between 1 and 2^31-1)"
    case(3)
       write(ERROR_UNIT,*) "Error in ranlxd_get"
       write(ERROR_UNIT,*) "Undefined state or improper argument"
    case(4)
       write(ERROR_UNIT,*) "Error in ranlxd_reset"
       write(ERROR_UNIT,*) "Arithmetic on this machine is not suitable for ranlxd"
    case(5)
       write(ERROR_UNIT,*) "Error in ranlxd_reset"
       write(ERROR_UNIT,*) "Unexpected input data"
    end select
    stop "Program aborted"
  end subroutine error

  !> @brief Updates the pseudo-random number generator
  !! @author Martin Lüscher <luscher@mail.desy.de>
  !! @date 15.12.1997
  !! @version 2.1
  impure subroutine update()

    integer(int64) :: k,kmax
    real(real64) :: y1,y2,y3

    k=0

    do 

       if (ir==0) then
          exit
       end if
       k=k+1

       y1=xdbl(jr)-xdbl(ir)
       y2=y1-carry
       if (y2<zero) then
          carry=one_bit
          y2=y2+one
       else
          carry=zero
       end if
       xdbl(ir)=y2

       ir=next(ir)
       jr=next(jr)

    end do

    kmax=pr-12

    do

       if (k>kmax) then
          exit
       end if
       k=k+12

       y1=xdbl(7)-xdbl(0)
       y1=y1-carry

       y2=xdbl(8)-xdbl(1)
       if (y1<zero) then
          y2=y2-one_bit
          y1=y1+one
       end if
       xdbl(0)=y1
       y3=xdbl(9)-xdbl(2)
       if (y2<zero) then
          y3=y3-one_bit
          y2=y2+one
       end if
       xdbl(1)=y2
       y1=xdbl(10)-xdbl(3)
       if (y3<zero) then
          y1=y1-one_bit
          y3=y3+one
       end if
       xdbl(2)=y3
       y2=xdbl(11)-xdbl(4)
       if (y1<zero) then
          y2=y2-one_bit
          y1=y1+one
       end if
       xdbl(3)=y1
       y3=xdbl(0)-xdbl(5)
       if (y2<zero) then
          y3=y3-one_bit
          y2=y2+one
       end if
       xdbl(4)=y2
       y1=xdbl(1)-xdbl(6)
       if (y3<zero) then
          y1=y1-one_bit
          y3=y3+one
       end if
       xdbl(5)=y3
       y2=xdbl(2)-xdbl(7)
       if (y1<zero) then
          y2=y2-one_bit
          y1=y1+one
       end if
       xdbl(6)=y1
       y3=xdbl(3)-xdbl(8)
       if (y2<zero) then
          y3=y3-one_bit
          y2=y2+one
       end if
       xdbl(7)=y2
       y1=xdbl(4)-xdbl(9)
       if (y3<zero) then
          y1=y1-one_bit
          y3=y3+one
       end if
       xdbl(8)=y3
       y2=xdbl(5)-xdbl(10)
       if (y1<zero) then
          y2=y2-one_bit
          y1=y1+one
       end if
       xdbl(9)=y1
       y3=xdbl(6)-xdbl(11)
       if (y2<zero) then
          y3=y3-one_bit
          y2=y2+one
       end if
       xdbl(10)=y2
       if (y3<zero) then
          carry=one_bit
          y3=y3+one
       else
          carry=zero
       end if
       xdbl(11)=y3

    end do

    kmax=pr-k

    do k=1,kmax

       y1=xdbl(jr)-xdbl(ir)
       y2=y1-carry
       if (y2<zero) then
          carry=one_bit
          y2=y2+one
       else
          carry=zero
       end if
       xdbl(ir)=y2

       ir=next(ir)
       jr=next(jr)

    end do

    ir_old=ir

  end subroutine update

  !> @brief Defines constants in initialisation step
  !! @author Martin Lüscher <luscher@mail.desy.de>
  !! @date 15.12.1997
  !! @version 2.1
  impure subroutine define_constants()

    integer(int64) :: k

    init=1
    zero=0.0_real64
    one=1.0_real64

!!!!!!!!!!!!!!!!!!!!!!!!
    ! this does not work with the pgi-compiler on a
    ! pc running linux. 
    ! pgi- "scale" only works for single precsision
    ! the NAGware f95 compiler works ok
    sbase=scale(one,24)
    base=scale(one,48)
    one_bit=scale(one,-48)

    ! if you have to use pgi use this instead:
    ! sbase=one*2**24
    ! base=one*2**48
    ! one_bit=one*2**(-48)
    ! this is done by redefining the scale function

    do k=0,11

       next(k)=modulo(k+1,12)

    end do

  end subroutine define_constants

  !> @brief Initialises Ranlux with given level and seed
  !! @author Martin Lüscher <luscher@mail.desy.de>
  !! @date 15.12.1997
  !! @version 2.1
  impure subroutine ranlxd_init(level,seed)
    !> Level for Ranlux
    integer(int64),intent(in) :: level
    !> Seed for Ranlux
    integer(int64),intent(in) :: seed

    integer(int64),dimension(0:30) :: xbit
    integer(int64) :: ibit,jbit,i,k,l
    real(real64) :: x,y 

    if ((huge(1)<2147483647).or.(radix(1.0_real64)/=2).or. &
         (digits(1.0_real64)<48)) then
       call error(0)
    end if

    select case(level)
    case(1)
       pr=202
    case(2)
       pr=397
    case default
       call error(1)
    end select

    call define_constants()

    i=seed

    do k=0,30

       xbit(k)=modulo(i,2)
       i=i/2

    end do

    if ((seed<=0).or.(i/=0)) then
       call error(2)
    end if

    ibit=0
    jbit=18

    do k=0,11

       x=zero

       do l=1,48

          y=real(modulo(xbit(ibit)+1,2),real64)
          x=x+x+y

          xbit(ibit)=modulo(xbit(ibit)+xbit(jbit),2)
          ibit=modulo(ibit+1,31)
          jbit=modulo(jbit+1,31)

       end do

       xdbl(k)=one_bit*x

    end do

    carry=zero
    ir=11
    jr=7
    ir_old=0

  end subroutine ranlxd_init

  !> @brief Draws a uniformly distributed real pseudo-random number with Ranlux
  !! @author Martin Lüscher <luscher@mail.desy.de>
  !! @date 15.12.1997
  !! @version 2.1
  impure subroutine ranlxd(r)

    real(real64),dimension(:),intent(out) :: r
    integer(int64) :: k,kmin,kmax

    if (init==0) then
       call ranlxd_init(1,1)
    end if

    kmin=lbound(r,1)
    kmax=ubound(r,1)

    do k=kmin,kmax

       ir=next(ir)

       if (ir==ir_old) then 

          call update()

       end if

       r(k)=xdbl(ir)

    end do

  end subroutine ranlxd

  !> @brief Gets the state of the pseudo-random number generator (Ranlux)
  !! @author Martin Lüscher <luscher@mail.desy.de>
  !! @date 15.12.1997
  !! @version 2.1
  pure subroutine ranlxd_get(state)

    integer(int64),dimension(:),intent(out) :: state

    integer(int64) :: k
    real(real64) :: x,y1,y2

    if ((init==0).or. &
         ((lbound(state,1)/=1).or.(ubound(state,1)/=25))) then
       !call error(3)
    end if

    do k=0,11

       x=sbase*xdbl(k)
       y2=aint(x,real64)
       y1=sbase*(x-y2)
       state(2*k+1)=int(y1)
       state(2*k+2)=int(y2)

    end do

    k=12*pr+ir
    k=12*k+jr
    k=12*k+ir_old
    state(25)=2*k+int(carry*base)

  end subroutine ranlxd_get

  !> @brief Resets the state of the pseudo-random number generator (Ranlux)
  !! @author Martin Lüscher <luscher@mail.desy.de>
  !! @date 15.12.1997
  !! @version 2.1
  impure subroutine ranlxd_reset(state)

    integer(int64),dimension(:),intent(in) :: state
    integer(int64) :: k
    real(real64) :: y1,y2

    if ((huge(1)<2147483647).or.(radix(1.0_real64)/=2).or. &
         (digits(1.0_real64)<48)) then
       call error(4)
    end if

    if ((lbound(state,1)/=1).or.(ubound(state,1)/=25)) then
       call error(5)
    end if

    call define_constants()

    do k=1,24

       if ((state(k)>=int(sbase)).or.(state(k)<0)) then
          call error(5)
       end if

    end do

    k=state(25)
    if (k<0) then
       call error(5)
    end if
    carry=one_bit*real(modulo(k,2),real64)
    k=k/2
    ir_old=modulo(k,12)
    k=k/12
    jr=modulo(k,12)
    k=k/12
    ir=modulo(k,12)
    pr=k/12

    if (((pr/=202).and.(pr/=397)).or.(jr/=modulo((ir_old+7),12))) then
       call error(5)
    end if

    do k=0,11

       y1=real(state(2*k+1),real64)
       y2=real(state(2*k+2),real64)
       xdbl(k)=one_bit*(y1+y2*sbase)

    end do

  end subroutine ranlxd_reset

end module ranlxd_generator



