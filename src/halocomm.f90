!----------------------------------------------------------------------
! Halo communcation module
!----------------------------------------------------------------------
!
! MODULE: halocomm
!> @brief Communication of halo values
!! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date 17.02.2019
!! @version 1.0
! REVISION HISTORY:
! 17 02 2019 - Initial version
!----------------------------------------------------------------------
module halocomm
  use, intrinsic :: iso_fortran_env
  implicit none

  PRIVATE

  public :: &
       InitModule,&
       IsModuleInitialised

  !> Module name
  character(len=8), parameter, public ::  modulename='halocomm'
  
  !> Contains information, if module is initialised
  logical :: IsInitialised = .false.
contains
  !> @brief Initialises module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use lattice, only: IsModuleInitialised_Lattice => IsModuleInitialised
    use mpi
    implicit none

    integer :: proc, mpierr

    ! ..--** START: Previous necessary initialisations **--..
    if(.not.IsModuleInitialised_Lattice()) then
       call mpi_comm_rank(MPI_COMM_WORLD, proc, mpierr)
       if(proc==0) then
          call flush(ERROR_UNIT)
          write(ERROR_UNIT,*) 'Error in init of ',modulename&
               ,': Lattice-module is not initialised.'
          call flush(ERROR_UNIT)
       end if
       STOP
    end if
    ! ..--**  END : Previous necessary initialisations **--..



    
    ! DONE
    IsInitialised = .TRUE.
  end subroutine InitModule

  !>@brief Returns, if module is initialised
  !! @returns module's initialisation status
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  pure logical function IsModuleInitialised()
    use, intrinsic :: iso_fortran_env
    implicit none
    IsModuleInitialised = isInitialised
  end function IsModuleInitialised

end module halocomm
