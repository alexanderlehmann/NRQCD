!----------------------------------------------------------------------
! PETSc-Interface
!----------------------------------------------------------------------
!
! MODULE: PETScInterface
!>@brief Initialisation, finalization, index routines, call to solver etc.
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 11.03.2019
!!@version 1.0
! REVISION HISTORY:
! 11 03 2019 - Initial version
!----------------------------------------------------------------------
module PETScInterface
#include <petsc/finclude/petscksp.h>
  use petscksp
  
  use, intrinsic :: iso_fortran_env
  use precision, only: fp

  implicit none
  
  PRIVATE

  public InitModule, FinalizeModule, &
       GetPETScIndex_G, GetLatticeIndex_PETSc, GetPETScIndex_M, GetMemoryIndex_PETSc


  !> Module name
  character(len=14), parameter, public ::  modulename='PETScInterface'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.

  !> Local PETSc indices
  PetscInt, allocatable :: LocalPETScIndices(:)
contains ! Module procedures

  !>@brief Initialises module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 08.03.2019
  !!@version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: mpistop

    use mpiinterface, only: thisproc, numprocs, syncall, intmpi
    use lattice

    implicit none

    ! *** TEST VARIABLES

    integer(int64) :: LatticeIndex
    ! ******************

    if(isInitialised) then
       call MPISTOP('Error in init of '//modulename//': already initialised.')
    else
       
       call CheckObligatoryInitialisations

       ! Initialise PETSc-solver
       call InitPETScSolver

       call InitPETScIndices(LocalPETScIndices)

       
       ! ***   Test   ***
       do LatticeIndex=1,GetLatticeSize()
          if(ThisProc()==GetProc_G(LatticeIndex)) then
             write(output_unit,*) LatticeIndex, GetPETScIndex_G(LatticeIndex)
             
          end if
          call flush(output_unit)
          call syncall
       end do
       
       ! DONE
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule

  !>@brief Finalizes module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 08.03.2019
  !!@version 1.0
  impure subroutine FinalizeModule
    use mpiinterface, only: mpistop
    implicit none

    if(isInitialised) then
       ! Clean up solver
       call FinalizePETScSolver

       IsInitialised = .FALSE.
    else
       call MPISTOP('Error in finalization of '//modulename//': is not initialised.')
    end if
  end subroutine FinalizeModule
  
  !>@brief Initialises PETSc-solver
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  impure subroutine InitPETScSolver
    use mpiinterface, only: mpistop
    implicit none

    PetscErrorCode PETScierr

    call PETScInitialize(PETSC_NULL_CHARACTER,PETScierr)

    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': Initialisation of PETSc failed.',&
            errorcode = PETScIerr)
    end if
  end subroutine InitPETScSolver

  !>@brief Finalizes PETSc-solver
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  impure subroutine FinalizePETScSolver
    use mpiinterface, only: mpistop
    implicit none

    PetscErrorCode PETScierr
    
    call PETScFinalize(PETScierr)

    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in finalization of '//&
            modulename//': Finalization of PETSc failed.',&
            errorcode = PETScIerr)
    end if
  end subroutine FinalizePETScSolver

  !>@brief Checks previous necessary initialisations
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 08.03.2019
  !!@version 1.0
  impure subroutine CheckObligatoryInitialisations
    use, intrinsic :: iso_fortran_env
    use Lattice,  only: IsLatticeInitialised  => IsModuleInitialised, LatticeName  => modulename
    use HaloComm, only: IsHaloCommInitialised => IsModuleInitialised, HaloCommName => modulename
    use mpiinterface, only: mpistop
    implicit none

    if(.not.IsLatticeInitialised()) then
       call mpistop('Error in init of '//modulename//': '//LatticeName//' is not initialised.')
    end if
    
    if(.not.IsHaloCommInitialised()) then
       call mpistop('Error in init of '//modulename//': '//HaloCommName//' is not initialised.')
    end if
  end subroutine CheckObligatoryInitialisations

  !>@brief Initialises local PETSc indices
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure subroutine InitPETScIndices(LocalPETScIndices)
    use lattice, only: GetMemorySize, GetLatticeIndex_M
    implicit none
    PetscInt, allocatable, intent(out) :: LocalPETScIndices(:)

    integer(int64) :: MemoryIndex
    
    allocate(LocalPETScIndices(GetMemorySize()))
    forall(MemoryIndex=1:GetMemorySize())
       LocalPETScIndices(MemoryIndex) = GetPETScIndex_G(GetLatticeIndex_M(MemoryIndex))
    end forall
  end subroutine InitPETScIndices
    
  !>@brief Returns PETSc index from lattice index
  !!@returns PETSc index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetPETScIndex_G(LatticeIndex)
    use lattice, only: GetLatticePosition, GetProc_G, &
         nDim, GetLocalLowerLatticeBoundary, GetLocalUpperLatticeBoundary, &
         GetIndex_fromPosition, GetLocalLatticeSize
    use mpiinterface, only: intmpi
    implicit none
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex

    ! Lattice related variables
    integer(int64), dimension(nDim) :: LocalLowerLatticeBoundaries, LocalUpperLatticeBoundaries
    integer(int64), dimension(nDim) :: LatticePosition
    integer(int64), dimension(nDim) :: LocalExtensions
    integer(intmpi) :: proc

    ! PETSc related variables
    integer(int64), dimension(nDim) :: LocalPosition

    LatticePosition = GetLatticePosition(LatticeIndex)
    proc = GetProc_G(LatticeIndex)

    LocalLowerLatticeBoundaries = GetLocalLowerLatticeBoundary([1_int8:nDim],proc)
    LocalUpperLatticeBoundaries = GetLocalUpperLatticeBoundary([1_int8:nDim],proc)
    LocalExtensions = LocalUpperLatticeBoundaries - LocalLowerLatticeBoundaries + 1
    
    LocalPosition = LatticePosition - LocalLowerLatticeBoundaries + 1

    GetPETScIndex_G = GetIndex_fromPosition(LocalPosition,LocalExtensions) &
         + proc * GetLocalLatticeSize()
  end function GetPETScIndex_G

  !>@brief Returns PETSc index from memory index
  !!@returns PETSc index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetPETScIndex_M(MemoryIndex)
    implicit none
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    GetPETScIndex_M = LocalPETScIndices(MemoryIndex)
  end function GetPETScIndex_M

  !>@brief Returns lattice index from PETSc index
  !!@returns lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetLatticeIndex_PETSc(PETScIndex)
    use lattice, only: GetLatticeIndex_M
    implicit none
    !> PETSc index
    integer(int64), intent(in) :: PETScIndex
    GetLatticeIndex_PETSc = GetLatticeIndex_M(GetMemoryIndex_PETSc(PETScIndex))
  end function GetLatticeIndex_PETSc

  !>@brief Returns memory index from PETSc index
  !!@returns memory index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetMemoryIndex_PETSc(PETScIndex)
    implicit none
    !> PETSc index
    integer(int64), intent(in) :: PETScIndex

    GetMemoryIndex_PETSc = FindLoc(&
         LocalPETScIndices,dim=1,& ! Where to look
         value = PETScIndex)       ! What to look for
  end function GetMemoryIndex_PETSc
  
end module PETScInterface
