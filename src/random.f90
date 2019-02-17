!------------------------------------------------------------------------------
! RANDOM, Random number generation interfaces and routines to Ranlux
!------------------------------------------------------------------------------
!
! MODULE: random
!> @brief Pseudo-random number generation interface using Ranlux
!! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date 17.02.2019
!! @version 1.0
! REVISION HISTORY:
! 03 09 2018 - Initial version
! 14 01 2019 - Added random number drawing from specific process
! 17 02 2019 - Added initialisation status
module random
  use, intrinsic :: iso_fortran_env
  use ranlxd_generator
  implicit none

  private
  public :: InitModule,&
       IsModuleInitialised,&
       GetState,&
       ResetState,&
       GetRandomUniformReal,&
       GetRandomUniformCmplx,&
       GetRandomNormalCmplx,&
       GetRandomNormalCmplx_specificProcess
  
  !> Contains information, if module is initialised
  logical :: IsInitialised = .false.
  
  !> @brief Getting an uniformly distributed real pseudo-random number
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  interface GetRandomUniformReal
     module procedure GetRandomUniformRealScalar
     module procedure GetRandomUniformRealArray
  end interface GetRandomUniformReal

  !> @brief Getting an uniformly distributed complex pseudo-random number
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  interface GetRandomUniformCmplx
     module procedure GetRandomUniformCmplxScalar
     module procedure GetRandomUniformCmplxArray
  end interface GetRandomUniformCmplx

  !> @brief Getting a normal distributed complex pseudo-random number
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  interface GetRandomNormalCmplx
     module procedure GetRandomNormalCmplxScalar
     module procedure GetRandomNormalCmplxArray
  end interface GetRandomNormalCmplx

  !> @brief Getting a normal distributed complex pseudo-random number from specific process
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  interface GetRandomNormalCmplx_specificProcess
     module procedure GetRandomNormalCmplxScalar_specificProcess
     module procedure GetRandomNormalCmplxArray_specificProcess
  end interface GetRandomNormalCmplx_specificProcess
contains
  
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
  
  !> @brief Current state of pseudo-random number generator
  !! @details
  !! The current state is necessary for, e.g., continuing a sequence of pseudo-random numbers (Ranlux)
  !! after closing the program
  !! @returns Current state of pseudo-random number generator (Ranlux)
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function GetState() result(state)
    implicit none
    !> State of pseudo-random number generator
    integer(int64), dimension(25) :: state
    call ranlxd_get(state)
  end function GetState

  !> @brief Sets the random number generator (Ranlux) to given state
  !! @details
  !! The current state is necessary for, e.g., continuing a sequence of pseudo-random numbers (Ranlux)
  !! after closing the program
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine ResetState(state)
    implicit none
    !> To-be set state of pseudo-random number generator
    integer(int64), dimension(25) :: state
    call ranlxd_reset(state)
    IsInitialised = .TRUE.
  end subroutine ResetState
  
  !> @brief Initialises pseudo-random number generator (Ranlux)
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.1
  impure subroutine InitModule(CustomSeed)
    implicit none
    !> seed for random number generator
    integer(int64), intent(in), optional :: CustomSeed
    integer(int64), parameter :: DefaultSeed=1_int64

    if(present(CustomSeed)) then
       call ranlxd_init(2,CustomSeed)
    else
       call ranlxd_init(2,DefaultSeed)
    end if
    
    ! DONE
    IsInitialised = .TRUE.
  end subroutine InitModule
  
  !> @brief Draws a real pseudo-random number (Ranlux) from uniform distribution in [0,1[
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function GetRandomUniformRealScalar() result(res)
    implicit none
    !> Pseudo-random number
    real(real64) :: res
    real(real64) :: r(1)
    call ranlxd(r)
    res = r(1)
  end function GetRandomUniformRealScalar
  
  !> @brief Draws real pseudo-random numbers (Ranlux) from uniform distribution in [0,1[
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function GetRandomUniformRealArray(n) result(res)
    implicit none
    !> Number of random numbers
    integer(int64), intent(in) :: n
    !> Pseudo-random numbers
    real(real64)               :: res(n)
    call ranlxd(res)
  end function GetRandomUniformRealArray
  
  !> @brief Draws a complex pseudo-random number (Ranlux) from uniform distribution in [0,1[
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function GetRandomUniformCmplxScalar() result(res)
    implicit none
    !> Pseudo-random number
    complex(real64) :: res
    real(real64)    :: r(2)
    call ranlxd(r)
    res = cmplx(r(1),r(2),real64)
  end function GetRandomUniformCmplxScalar
  
  !> @brief Draws complex pseudo-random numbers (Ranlux) from uniform distribution in ]0,1]
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure function GetRandomUniformCmplxArray(n) result(res)
    implicit none
    !> Number of random numbers
    integer(int64), intent(in) :: n
    !> Pseudo-random numbers    
    complex(real64) :: res(n)
    
    real(real64)    :: r(2*n)

    integer(int64) :: i

    call ranlxd(r)
    do concurrent(i=1:n)
       res(i) = cmplx(r(2*i-1),r(2*i),real64)
    end do
  end function GetRandomUniformCmplxArray
    
  !> @brief Draws complex pseudo-random number (Ranlux) from normal distribution
  !! @details
  !! Uses Box-M\"uller transform on normal distributed random numbers (Ranlux)
  !! \f$r_{\text{u}}\in]0,1]^2\f$:
  !! \f$r_{\text{normal}}=\sqrt{-2\cdot\log(r_1)}\cdot(\cos(2\pi r_2)+\text{i}\sin(2\pi r_2))\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function GetRandomNormalCmplxScalar() result(res)
    use mathconstants, only: pi
    implicit none
    !> Pseudo-random number
    complex(real64) :: res
    real(real64)    :: r(2)
    real(real64), parameter :: twopi = 2*pi
    call ranlxd(r)
    res = sqrt(-2._real64 * log(r(1)))&
         * cmplx(&
         cos(twopi*r(2)),&
         sin(twopi*r(2)),real64)
  end function GetRandomNormalCmplxScalar
  
  !> @brief Draws n complex pseudo-random numbers (Ranlux) from normal distribution
  !! @details
  !! Uses Box-M\"uller transform on normal distributed random numbers (Ranlux)
  !! \f$r_{\text{u}}\in]0,1]^2\f$:
  !! \f$r_{\text{normal}}=\sqrt{-2\cdot\log(r_1)}\cdot(\cos(2\pi r_2)+\text{i}\sin(2\pi r_2))\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function GetRandomNormalCmplxArray(n) result(res)
    use mathconstants, only: pi
    implicit none
    !> Number of pseudo-random numbers
    integer(int64), intent(in) :: n
    !> Pseudo-random numbers
    complex(real64) :: res(n)
    real(real64)    :: r(2*n)
    integer(int64) :: i
    real(real64), parameter :: twopi = 2*pi
    call ranlxd(r)
    do concurrent(i=1:n)
       res(i) = sqrt(-2._real64 * log(r(2*i-1)))&
            * cmplx(&
            cos(twopi*r(2*i)),&
            sin(twopi*r(2*i)),real64)
    end do
  end function GetRandomNormalCmplxArray

  !> @brief
  !! Draws n complex pseudo-random numbers (Ranlux) from normal distribution on one specific process\n
  !! and sends it to another recieving one
  !! @details
  !! Uses Box-M\"uller transform on normal distributed random numbers (Ranlux)
  !! \f$r_{\text{u}}\in]0,1]^2\f$:
  !! \f$r_{\text{normal}}=\sqrt{-2\cdot\log(r_1)}\cdot(\cos(2\pi r_2)+\text{i}\sin(2\pi r_2))\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 16.02.2019
  !! @version 1.0
  impure function GetRandomNormalCmplxScalar_specificProcess(send_proc,recv_proc) result(res)
    use mpiinterface, only: ThisProc
    use mpi
    implicit none
    !> Sending MPI-process
    integer(int64), intent(in) :: send_proc
    !> Recieving MPI-process
    integer(int64), intent(in) :: recv_proc
    !> Pseudo-random numbers
    complex(real64) :: res

    ! MPI
    integer :: status(mpi_status_size)
    integer :: buffersize
    integer :: dest, source
    integer :: tag
    integer :: mpierr
    
    if(send_proc==recv_proc) then
       if(ThisProc()==send_proc) res = GetRandomNormalCmplxScalar()
    else
       buffersize = 1
       dest       = recv_proc
       source     = send_proc
       tag = recv_proc
       if(ThisProc()==send_proc) then
          res = GetRandomNormalCmplxScalar()
       
          call MPI_SEND(res,buffersize,mpi_double_complex,dest,tag,&
               MPI_COMM_WORLD,mpierr)
          
       elseif(ThisProc()==recv_proc) then
          call MPI_RECV(res,buffersize,mpi_double_complex,source,tag,&
               MPI_COMM_WORLD,status,mpierr)
       end if
    end if
  end function GetRandomNormalCmplxScalar_specificProcess

  !> @brief
  !! Draws n complex pseudo-random numbers (Ranlux) from normal distribution on one specific process\n
  !! and sends it to another recieving one
  !! @details
  !! Uses Box-M\"uller transform on normal distributed random numbers (Ranlux)
  !! \f$r_{\text{u}}\in]0,1]^2\f$:
  !! \f$r_{\text{normal}}=\sqrt{-2\cdot\log(r_1)}\cdot(\cos(2\pi r_2)+\text{i}\sin(2\pi r_2))\f$
  !! @authorAlexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 16.02.2019
  !! @version 1.1
  impure function GetRandomNormalCmplxArray_specificProcess(n,send_proc,recv_proc) result(res)
    use mpiinterface, only: ThisProc
    use mpi
    implicit none
    !> Number of pseudo-random numbers
    integer(int64), intent(in) :: n
    !> Sending MPI-process
    integer(int64), intent(in) :: send_proc
    !> Recieving MPI-process
    integer(int64), intent(in) :: recv_proc
    !> Pseudo-random numbers
    complex(real64) :: res(n)

    ! MPI
    integer :: status(mpi_status_size)
    integer :: buffersize
    integer :: dest, source
    integer :: tag
    integer :: mpierr

    if(ThisProc()==recv_proc .or. ThisProc()==send_proc) then
       if(send_proc==recv_proc) then
          res = GetRandomNormalCmplxArray(n)
       else
          buffersize = n
          dest       = recv_proc
          source     = send_proc
          tag = recv_proc
          if(ThisProc()==send_proc) then
             ! This process is sending
             res = GetRandomNormalCmplxArray(n)

             call MPI_SEND(res,buffersize,mpi_double_complex,dest,tag,&
                  MPI_COMM_WORLD,mpierr)

          else
             ! This process is receiving
             call MPI_RECV(res,buffersize,mpi_double_complex,source,tag,&
                  MPI_COMM_WORLD,status,mpierr)
          end if
       end if
    end if
  end function GetRandomNormalCmplxArray_specificProcess
  
end module random
