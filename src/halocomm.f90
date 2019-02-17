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

  !> List of points and to this process sending processes
  integer(int64), allocatable :: GetList(:,:)
  !> List of points and from this process recieving processes
  integer(int64), allocatable :: PutList(:,:)
contains
  !> @brief Initialises module
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    implicit none

    call CheckDependencies
    
    ! Initialise list of lattice points which are to be recieved from which other process
    call InitGetList(GetList)

    ! DONE
    IsInitialised = .TRUE.
    
  end subroutine InitModule
  
  !> @brief Checks previous necessary initialisations
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine CheckDependencies
    use, intrinsic :: iso_fortran_env
    use lattice, only: IsModuleInitialised_Lattice => IsModuleInitialised
    use mpi
    implicit none

    integer :: proc, mpierr
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
  end subroutine CheckDependencies
  
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

  !>@brief Initiales list of processes and points which send to this process
  !! @returns list of processes and points which send to this process
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.02.2019
  !! @version 1.0
  impure subroutine InitGetList(GetList)
    use lattice, only: GetLocalLatticeIndices_includingHalo_Allocatable, GetProc
    use mpiinterface, only: NumProcs, ThisProc
    use arrayoperations, only: RemoveDuplicates, Sort
    implicit none
    integer(int64), allocatable, intent(out) :: GetList(:,:)

    integer(int64), allocatable :: latticeindices(:)
    integer, allocatable :: procs(:), procs_copy(:)
    integer(int8) :: neighbours, i
    integer(int64), allocatable :: points_per_proc(:)

    call GetLocalLatticeIndices_includingHalo_Allocatable(latticeindices)
    allocate(procs(size(latticeindices)))
    procs = GetProc(latticeindices)
    allocate(procs_copy(size(procs)))
    procs_copy=procs

    call RemoveDuplicates(procs_copy,points_per_proc)
    call Sort(procs_copy)

    ! Count number of different neighbours
    neighbours = size(procs_copy,int8) - 1_int8

    if(ThisProc()==0) then
       do i=1,neighbours
          print*,procs_copy(i),points_per_proc(i)
       end do
    end if
  end subroutine InitGetList
end module halocomm
