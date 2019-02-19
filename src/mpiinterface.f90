!----------------------------------------------------------------------
! mpiInterface
!----------------------------------------------------------------------
!
! MODULE: mpiInterface
!> @brief Convenience-interface to most needed MPI-routines
!! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date 17.02.2019
!! @version 1.0
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
       MPIstop

  !> Module name
  character(len=12), parameter, public ::  modulename='mpiinterface'
  
  !> Contains information, if module is initialised
  logical :: IsInitialised = .false.
  
  !> Number of processes
  integer, private :: num_procs=-1
  !> Process number of this process
  integer, private :: this_proc=-1

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
  
  !>@brief Initialization of module
  !! @details Initialises MPI (call to MPI_INIT), getting of process' number
  !! as well as total number of processes
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    implicit none
    integer :: mpierr

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
       
    ! DONE
    IsInitialised = .TRUE.
  end subroutine InitModule

  !>@brief Finalization of module
  !! @details Finalises MPI (call to MPI_FINALIZE)
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine FinalizeModule
    implicit none
    integer :: mpierr
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
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine MPIstop(errormessage,errorcode)
    use, intrinsic :: iso_fortran_env
    implicit none
    character(len=*), intent(in), optional :: errormessage
    integer,          intent(in), optional :: errorcode
    integer :: mpierr

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
  !! @returns Total number of processes
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  pure integer function NumProcs()
    implicit none
    NumProcs = num_procs
  end function NumProcs

  !>@brief Getting ID of this process
  !! @returns ID of this process
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  pure integer function ThisProc()
    implicit none
    ThisProc = this_proc
  end function ThisProc
end module mpiInterface
