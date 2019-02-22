!----------------------------------------------------------------------
! mpiInterface
!----------------------------------------------------------------------
!
! MODULE: mpiInterface
!>@brief Interface to some of MPI's variables like the rank and number of processes
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 17.02.2019
!!@version 1.0
! REVISION HISTORY:
! 15 02 2019 - Initial version
module mpiInterface
  use, intrinsic :: iso_fortran_env
  use mpi
  implicit none

  PRIVATE

  public &
       InitModule,&
       FinalizeModule,&
       IsModuleInitialised,&
       ThisProc,&
       NumProcs,&
       MPIstop,&
       GetComplexSendType,&
       GetRealSendType

  !> Module name
  character(len=12), parameter, public ::  modulename='mpiinterface'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.

  !> Integer-kind of mpi-variables like MPI_COMM_WORLD, process ranks and so on
  integer(int8), parameter, public :: intmpi = kind(MPI_COMM_WORLD)
  
  !> Number of processes
  integer(intmpi), private :: num_procs=-1_intmpi
  !> Process number of this process
  integer(intmpi), private :: this_proc=-1_intmpi
  !> Message type of real variables with precision in precision-module
  integer(intmpi), private :: realsendtype_fp=-1_intmpi
  !> Message type of complex variables with precision in precision-module
  integer(intmpi), private :: complexsendtype_fp=-1_intmpi
contains
  
  !>@brief Returns, if module is initialised
  !!@returns module's initialisation status
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure logical function IsModuleInitialised()
    implicit none
    IsModuleInitialised = IsInitialised
  end function IsModuleInitialised
  
  !>@brief Initialization of module
  !!@details Initialises MPI (call to MPI_INIT), getting of process' number
  !! as well as total number of processes
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  impure subroutine InitModule
    use precision, only: fp
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(intmpi) :: mpierr

    character(len=100) :: errormessage

    if(isInitialised) then
       errormessage = 'Error in init of '//modulename//': already initialised.'
       call MPISTOP(errormessage)
    else
       call MPI_INIT(mpierr)
       if(mpierr /= MPI_SUCCESS) then
          write(ERROR_UNIT,*) 'Error while MPI-initialization. Error code:',mpierr
          STOP
       end if

       call mpi_comm_rank(MPI_COMM_WORLD, this_proc, mpierr)
       if(mpierr /= MPI_SUCCESS) then
          write(ERROR_UNIT,*) 'Error while getting MPI-rank. Error code:',mpierr
          STOP
       end if

       call mpi_comm_size(MPI_COMM_WORLD, num_procs, mpierr)
       if(mpierr /= MPI_SUCCESS) then
          write(ERROR_UNIT,*) 'Error while getting number of MPI-processes. Error code:',mpierr
          STOP
       end if
       
       select case(fp)
       case(real32)
          realsendtype_fp    = MPI_REAL4
          complexsendtype_fp = MPI_COMPLEX8
       case(real64)
          realsendtype_fp    = MPI_REAL8
          complexsendtype_fp = MPI_COMPLEX16
       case(real128)
          realsendtype_fp    = MPI_REAL16
          complexsendtype_fp = MPI_COMPLEX32
       case default
          errormessage = 'Error in GetDeviationFromGaussLaw of '//modulename&
               //': unsupported floating point precision.'
          call MPISTOP(errormessage)
       end select
    
       ! DONE
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule

  !>@brief Finalization of module
  !!@details Finalises MPI (call to MPI_FINALIZE)
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  impure subroutine FinalizeModule
    implicit none
    integer(intmpi) :: mpierr
    character(len=70) :: errormessage
    
    if(IsInitialised) then
       call MPI_FINALIZE(mpierr)
       IsInitialised = .FALSE.
    else
       errormessage = 'Error in finalization of '//modulename//': is not initialised.'
       call MPISTOP(errormessage)
    end if
  end subroutine FinalizeModule

  !>@brief MPI stop with optional message and code
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  impure subroutine MPIstop(errormessage,errorcode)
    use, intrinsic :: iso_fortran_env
    implicit none
    character(len=*),  intent(in), optional :: errormessage
    integer(intmpi), intent(in), optional :: errorcode
    integer(intmpi) :: mpierr

    call flush(6)
    call flush(ERROR_UNIT)
    if(.not.isInitialised) call InitModule
    call mpi_barrier(mpi_comm_world,mpierr)
    call mpi_comm_rank(MPI_COMM_WORLD, this_proc, mpierr)
    if(this_proc==0) then
       if(present(errormessage)) write(ERROR_UNIT,*) trim(errormessage)
       if(present(errorcode))    write(ERROR_UNIT,*)'Error code:',errorcode
    end if
    call flush(ERROR_UNIT)
    call mpi_barrier(mpi_comm_world,mpierr)
    call MPI_Finalize(mpierr)

    if(this_proc==0) then
       STOP "Program aborted"
    else
       STOP
    end if
  end subroutine MPIstop


  !>@brief Getting total number of processes
  !!@returns Total number of processes
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure integer(intmpi) function NumProcs()
    implicit none
    NumProcs = num_procs
  end function NumProcs

  !>@brief Getting ID of this process
  !!@returns ID of this process
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  pure integer(intmpi) function ThisProc()
    implicit none
    ThisProc = this_proc
  end function ThisProc
  
  !>@brief Complex send type
  !!@returns Complex send type
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure integer(intmpi) function GetComplexSendType(kind)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, intent(in), optional :: kind

    if(.not.present(kind)) then
       GetComplexSendType = complexsendtype_fp
    else
       select case(kind)
       case(real32)
          GetComplexSendType = MPI_COMPLEX8
       case(real64)
          GetComplexSendType = MPI_COMPLEX16
       case(real128)
          GetComplexSendType = MPI_COMPLEX32
       end select
    end if
  end function GetComplexSendType

  !>@brief Real send type
  !!@returns Real send type
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure integer(intmpi) function GetRealSendType(kind)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, intent(in), optional :: kind

    if(.not.present(kind)) then
       GetRealSendType = realsendtype_fp
    else
       select case(kind)
       case(real32)
          GetRealSendType = MPI_REAL4
       case(real64)
          GetRealSendType = MPI_REAL8
       case(real128)
          GetRealSendType = MPI_REAL16
       end select
    end if
  end function GetRealSendType
end module mpiInterface
