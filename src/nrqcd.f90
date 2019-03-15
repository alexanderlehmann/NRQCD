!-------------------------------------------------------------------------
! NRQCD, Non-Relativistic Quantumchromodynamics simulation of heavy quarks
!-------------------------------------------------------------------------
!
! MODULE: nrqcd
!>@brief Heavy quarks, described by NRQCD in a real-time evolution scheme
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 06.03.2019
!!@version 1.0
! REVISION HISTORY:
! 06 03 2019 - Initial Version
!-------------------------------------------------------------------------
module nrqcd
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  use gaugeconfiguration_su3, only: &
       SU3GaugeConfiguration => GaugeConfiguration
  use su3, only: &
       SU3Generators => Generators, &
       nColours      => nSUN, &
       nGluons       => nGen
  use su2, only: &
       SU2Generators => Generators, &
       nSpins        => nSUN

#include <petsc/finclude/petscksp.h>
  use petscksp
  
  implicit none

  PRIVATE

  public NRQCDField, InitModule, FinalizeModule

  !> Module name
  character(len=5), parameter, public ::  modulename='nrqcd'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.

  !> Number of degrees of freedom per site
  integer(int8), parameter, public :: nDof = nSpins*nColours
  !> Number of Wilson coefficients
  integer(int8), parameter, public :: nWilsonCoefficients = 4

  !> Constructor for NRQCDField
  interface NRQCDField
     module procedure NRQCDField_Constructor
  end interface NRQCDField
  
  !>@brief Heavy quark, described by NRQCD, containing quark- and antiquark-propagators
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 13.03.2019
  !!@version 1.0
  type NRQCDField
     !> Heavy-quark-propagator
     complex(fp), allocatable, private :: QuarkProp(:,:,:)

     !> Heavy antiquark-propagator
     complex(fp), allocatable, private :: AntiQprop(:,:,:)

     !> PETSc system matrix
     Mat :: SystMat
     !> Solution vector
     Vec :: x
     !> RHS vector
     Vec :: rhs
     !> Solver context
     KSP :: ksp
   contains ! Member functions
     ! Allocation and deallocation
     procedure, private :: Allocate
     procedure, private :: Deallocate
     
     ! Destructor
     procedure, public  :: Destructor

     ! Assignment
     procedure, public :: Assign
     generic :: assignment(=) => Assign

     ! Boundary-communication
     procedure, public :: CommunicateBoundary
     
     ! Initialisation routines
     procedure, public :: InitSinglePoint

     ! Norm
     procedure, public :: GetNorm_Quark
     procedure, public :: GetNorm_AntiQ
     
  end type NRQCDField

  !> Local PETSc indices
  PetscInt, allocatable :: LocalSpatialPETScIndices(:)
  !> Non-zero entries per row
  PetscInt, parameter :: nz = 25*ndof
  !> Local number of diagonal non-zero entries per row
  PetscInt, allocatable :: d_nnz(:)
  !> Local number of off-diagonal non-zero entries per row
  PetscInt, allocatable :: o_nnz(:)
  !> Total number of rows and cols
  PetscInt :: SystMatSize
  !> Local number of rows
  PetscInt :: nLocalRows
  !> Local number of cols and rows
  PetscInt :: LocalSystMatSize
  !> Local minimum of row index (= diagonal section)
  PetscInt :: LocalMinRow
  !> Local maximum of row index (= diagonal section)
  PetscInt :: LocalMaxRow
  !> Local maximum of spatial row index (= diagonal section)
  PetscInt :: LocalSpatialMaxRow
  !> Local minimum of spatial row index (= diagonal section)
  PetscInt :: LocalSpatialMinRow
  !> MPI-rank
  PetscMPIInt :: Rank
  !> Total number of ranks
  PetscMPIInt :: NumRanks
  
contains ! Module procedures

  !>@brief Initialises module
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 13.03.2019
  !!@version 1.0
  impure subroutine InitModule
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: mpistop
    implicit none

    if(isInitialised) then
       call MPISTOP('Error in init of '//modulename//': already initialised.')
    else
       
       call CheckObligatoryInitialisations

       ! Initialise PETSc-solver
       call InitPETScSolver
       
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

  !>@brief Assignment
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 13.03.2019
  !!@version 1.0
  impure subroutine assign(to,from)
    implicit none
    !> left hand side, to
    class(NRQCDField), intent(out) :: to
    !> right hand side, from
    class(NRQCDField), intent(in)  :: from

    allocate(to%QuarkProp,source=from%QuarkProp)
    allocate(to%AntiQProp,source=from%AntiQProp)

    to%ksp     = from%ksp
    to%SystMat = from%SystMat
    to%x       = from%x
    to%rhs     = from%rhs
  end subroutine assign
  
  !>@brief Constructor
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 13.03.2019
  !!@version 1.0
  impure function NRQCDField_Constructor()
    use mpiinterface, only: intmpi,ThisProc, NumProcs,SyncAll,mpistop
    implicit none
    
    type(NRQCDField) :: NRQCDField_Constructor

    PetscMPIInt :: PETScIerr
    
    call NRQCDField_Constructor%Allocate

    ! PETSc variables for linear solver
    ! 1. System matrix
    call MatCreateAIJ(PETSC_COMM_WORLD,&
         LocalSystMatSize,LocalSystMatSize,&
         SystMatSize,SystMatSize,&
         0_int64,d_nnz,0_int64,o_nnz,&
         NRQCDField_Constructor%SystMat,&
         PETScIerr)
    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': MatCreateAIJ failed.',&
            errorcode = PETScIerr)
    end if

    ! 2. Solution vector
    call VecCreateMPI(PETSC_COMM_WORLD,&
         LocalSystMatSize,SystMatsize,&
         NRQCDField_Constructor%x,PETScIerr)
    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': VecCreateMPI for x failed.',&
            errorcode = PETScIerr)
    end if

    ! 3. RHS vector
    call VecCreateMPI(PETSC_COMM_WORLD,&
         LocalSystMatSize,SystMatsize,&
         NRQCDField_Constructor%rhs,PETScIerr)
    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': VecCreateMPI for rhs failed.',&
            errorcode = PETScIerr)
    end if

    ! 4. KSP-solver object
    call KSPCreate(PETSC_COMM_WORLD,NRQCDField_Constructor%ksp,PETScIerr)
    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': KSPCreate failed.',&
            errorcode = PETScIerr)
    end if
  end function NRQCDField_Constructor

  
  !>@brief Destructor
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 15.03.2019
  !!@version 1.0
  impure subroutine Destructor(object)
    implicit none
    class(NRQCDField), intent(inout) :: object
    
    PetscErrorCode :: PETScIerr
    
    call VecDestroy(object%x,PETScIerr)
    call VecDestroy(object%rhs,PETScIerr)
    call MatDestroy(object%SystMat,PETScIerr)
    call KSPDestroy(object%ksp,PETScIerr)
    
    call Object%Deallocate
  end subroutine Destructor
  
  !>@brief Allocation of NRQCD heavy quark
  !!@details Allocates quark and antiquark propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 06.03.2019
  !!@version 1.0
 pure subroutine Allocate(HeavyField)
    use lattice, only: nDim, GetMemorySize
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField

    allocate(HeavyField%QuarkProp(nDof,nDof,GetMemorySize()))
    allocate(HeavyField%AntiQProp(nDof,nDof,GetMemorySize()))
  end subroutine Allocate
  
  !>@brief Deallocation of NRQCD heavy quark
  !!@details Deallocates quark and antiquark propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 06.03.2019
  !!@version 1.0
 pure subroutine Deallocate(HeavyField)
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField
    if(allocated(HeavyField%QuarkProp))  deallocate(HeavyField%QuarkProp)
    if(allocated(HeavyField%AntiQProp))  deallocate(HeavyField%AntiQProp)    
  end subroutine Deallocate

  !>@brief Communication routine for boundary values in NRQCD heavy field
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 06.03.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary(HeavyField)
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(inout) :: HeavyField
    call CommunicateBoundary_Propagator(HeavyField%QuarkProp)
    call CommunicateBoundary_Propagator(HeavyField%AntiQProp)
  end subroutine CommunicateBoundary

  !>@brief Communication routine for boundary values of propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 06.03.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_Propagator(Propagator)
    use precision, only: fp
    use HaloComm, only: HaloComm_CommunicateBoundary => CommunicateBoundary
    implicit none
    complex(fp), intent(inout) :: Propagator(:,:,:)
    call HaloComm_CommunicateBoundary(Propagator)
  end subroutine CommunicateBoundary_Propagator

  !>@brief Initialises the propagator and anti-propagator
  !! at a single spatial point for a single colour- and spin-combination
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine InitSinglePoint(Heavyfield,spin,colour,LatticeIndex)
    use mpiinterface, only: ThisProc
    use lattice, only: GetProc_G, GetMemoryIndex
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField
    !> Spin-component of created quark-antiquark-pair
    integer(int8),     intent(in)  :: spin
    !> Colour-component of created quark-antiquark-pair
    integer(int8),     intent(in)  :: colour
    !> Lattice index
    integer(int64),    intent(in)  :: LatticeIndex
    
    !call HeavyField%Allocate
    HeavyField = NRQCDField()
    
    HeavyField%QuarkProp = 0
    HeavyField%AntiQProp = 0

    if(ThisProc()==GetProc_G(LatticeIndex)) then
       HeavyField%QuarkProp(&
            GetSpinColourIndex(Spin,Colour),&
            GetSpinColourIndex(Spin,Colour),&
            GetMemoryIndex(LatticeIndex)) = 1
       
       HeavyField%AntiQProp(&
            GetSpinColourIndex(Spin,Colour),&
            GetSpinColourIndex(Spin,Colour),&
            GetMemoryIndex(LatticeIndex)) = 1
    end if

    call HeavyField%CommunicateBoundary
  end subroutine InitSinglePoint

  !>@brief Combined spin-colour-index
  !!@returns Combined spin-colour-index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  elemental integer(int8) function GetSpinColourIndex(spin,colour)
    implicit none
    !> Spin-component \f$s\in{1,2}\equiv\{+1/2,-1/2\}\f$
    integer(int8), intent(in) :: spin
    !> Colour-component \f$c\in{1,2,3}\equiv\{\text{red},\text{green},\text{blue}\}\f$
    integer(int8), intent(in) :: colour
    GetSpinColourIndex = colour + ncolours*(spin-1_int8)
  end function GetSpinColourIndex

  !>@brief Raises a matrix from colour-space to colour-spin-space
  !!@returns Kronecker-product of colour-matrix (3x3) with unit matrix (2x2) to colour-spin-space
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  pure function C2CS(input)
    use precision, only: fp
    use matrixoperations, only: GetUnitMatrix, GetKronProd
    implicit none
    complex(fp), intent(in) :: input(ncolours,ncolours)
    complex(fp)             :: C2CS(ndof,ndof)

    complex(fp), parameter ::&
         Unitmatrix(nspins,nspins) &
         = reshape([&
         cmplx(1,0,fp),cmplx(0,0,fp),&
         cmplx(0,0,fp),cmplx(1,0,fp)],&
         [nspins,nspins])

    C2CS = GetKronProd(input,unitmatrix)
  end function C2CS

  !>@brief Raises a matrix from spin-space to colour-spin-space
  !!@returns Kronecker-product of spin-matrix (2x2) with unit matrix (3x3) to colour-spin-space
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  pure function S2CS(input)
    use precision, only: fp
    use matrixoperations, only: GetUnitMatrix, GetKronProd
    implicit none
    complex(fp), intent(in) :: input(nspins,nspins)
    complex(fp)             :: S2CS(ndof,ndof)

    complex(fp), parameter ::&
         Unitmatrix(ncolours,ncolours) &
         = reshape([&
         cmplx(1,0,fp),cmplx(0,0,fp),cmplx(0,0,fp),&
         cmplx(0,0,fp),cmplx(1,0,fp),cmplx(0,0,fp),&
         cmplx(0,0,fp),cmplx(0,0,fp),cmplx(1,0,fp)],&
         [ncolours,ncolours])

    S2CS = GetKronProd(unitmatrix,input)
  end function S2CS

  !>@brief Returns link represented in colour-spin-space
  !!@returns Link represented in colour-spin-space
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  pure function GetLinkCS_M(GaugeConf,i,MemoryIndex)
    use precision, only: fp
    implicit none
    !> Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),               intent(in) :: i
    !> Memory index
    integer(int64),              intent(in) :: MemoryIndex
    !> Link variable in colour x spin - space
    complex(fp) :: GetLinkCS_M(ndof,ndof)

    complex(fp) :: link(ncolours,ncolours)

    link = GaugeConf%GetLink_M(i,MemoryIndex)
    GetLinkCS_M = C2CS(link)
  end function GetLinkCS_M

  !>@brief Returns link represented in colour-spin-space
  !!@returns Link represented in colour-spin-space
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  pure function GetLinkCS_G(GaugeConf,i,LatticeIndex)
    use precision, only: fp
    implicit none
    !> Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),               intent(in) :: i
    !> Lattice index
    integer(int64),              intent(in) :: LatticeIndex
    !> Link variable in colour x spin - space
    complex(fp) :: GetLinkCS_G(ndof,ndof)

    complex(fp) :: link(ncolours,ncolours)

    link = GaugeConf%GetLink_G(i,LatticeIndex)
    GetLinkCS_G = C2CS(link)
  end function GetLinkCS_G

  !>@brief Norm of the quark-propagator
  !!@returns Norm of the quark-propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure real(fp) function GetNorm_Quark(HeavyField)
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(in) :: HeavyField
    GetNorm_Quark = GetNorm_Propagator(HeavyField%QuarkProp)
  end function GetNorm_Quark

  !>@brief Norm of the anti-quark-propagator
  !!@returns Norm of the anti-quark-propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure real(fp) function GetNorm_AntiQ(HeavyField)
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(in) :: HeavyField
    GetNorm_AntiQ = GetNorm_Propagator(HeavyField%AntiQProp)
  end function GetNorm_AntiQ
  
  !>@brief Norm of the anti-quark-propagator
  !!@returns Norm of the anti-quark-propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure real(fp) function GetNorm_Propagator(propagator)
    use mpi
    use matrixoperations, only: GetTrace
    use mpiinterface, only: intmpi, ThisProc, GetRealSendType
    use lattice, only : GetProc_M, GetMemorySize
    implicit none
    !> NRQCD-Quark- or Antiquark-Propagator
    complex(fp), intent(in) :: propagator(:,:,:)

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    complex(fp) :: prop_squared(ndof,ndof)

    integer(int64) :: MemoryIndex
    
    ! 1. Calculation of local contribution
    local_contribution = 0
    do concurrent(MemoryIndex=1:GetMemorySize(), ThisProc()==GetProc_M(MemoryIndex))
       prop_squared = real(matmul(&
            propagator(:,:,MemoryIndex),&
            conjg(transpose(propagator(:,:,MemoryIndex)))),fp)
       
       local_contribution = local_contribution &
            + GetTrace(prop_squared)
    end do

    ! 2. MPI-Sum over all partitions
    call MPI_ALLREDUCE(&
         local_contribution,&
         GetNorm_Propagator,&
         1_intmpi,&
         GetRealSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)

    GetNorm_Propagator = sqrt(GetNorm_Propagator)

  end function GetNorm_Propagator


  ! **PETSc**
  !>@brief Initialises PETSc-solver
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  impure subroutine InitPETScSolver
    use mpiinterface, only: intmpi,mpistop,thisproc,numprocs
    use lattice, only: GetLatticeSize, GetLocalLatticeSize, GetMemorySize, GetNeib_G, GetProc_G
    use arrayoperations, only: RemoveDuplicates
    implicit none

    PetscErrorCode :: PETScierr

    integer(int64) :: SpatialRow, SpatialCol, LatticeIndex, i
    
    integer(int64), allocatable :: neibs(:)
    
    integer(int64) :: dofMinRow,dofMaxRow
    
    call PETScInitialize(PETSC_NULL_CHARACTER,PETScierr)

    if(PETScIerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': Initialisation of PETSc failed.',&
            errorcode = PETScIerr)
    end if

    call InitLocalSpatialPETScIndices(LocalSpatialPETScIndices)
    
    SystMatSize = GetLatticeSize()*ndof
    SystMatSize = SystMatSize
    LocalSystMatSize = GetLocalLatticeSize()*ndof

    LocalSpatialMinRow = 1+ ThisProc()   *GetLocalLatticeSize()
    LocalSpatialMaxRow =   (ThisProc()+1)*GetLocalLatticeSize()
    LocalMinRow = 1 + ThisProc()   *GetLocalLatticeSize()*nDof
    LocalMaxRow =    (ThisProc()+1)*GetLocalLatticeSize()*nDof
    
    ! Determine where the rows and columns belong to and thus the diagonal and off-diagonal entries
    
    rank = ThisProc()
    numranks = NumProcs()

    allocate(d_nnz(LocalMinRow:LocalMaxRow))
    d_nnz = 0
    allocate(o_nnz,source=d_nnz)

    ! For each spatial row determine the neighbouring points (regarding Laplacian^2),
    ! the correspoding "PETSc"-column and the process to which the column belongs
    ! --> Counting number of diagonal and off-diagonal entries
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       LatticeIndex = GetLatticeIndex_PETSc(SpatialRow)
       
       !1. Computing all neighbouring lattice sites
       if(allocated(neibs)) deallocate(neibs)
       allocate(neibs(25))
       
       ! ∇⁰
       neibs(1) = GetNeib_G(0_int8,LatticeIndex)
       ! ∇¹ and ∇²
       neibs(2) = GetNeib_G(+1_int8,LatticeIndex)
       neibs(3) = GetNeib_G(+2_int8,LatticeIndex)
       neibs(4) = GetNeib_G(+3_int8,LatticeIndex)
       neibs(5) = GetNeib_G(-1_int8,LatticeIndex)
       neibs(6) = GetNeib_G(-2_int8,LatticeIndex)
       neibs(7) = GetNeib_G(-3_int8,LatticeIndex)
       ! ∇⁴
       neibs(8) = GetNeib_G(+1_int8,GetNeib_G(+1_int8,LatticeIndex))
       neibs(9) = GetNeib_G(+1_int8,GetNeib_G(+2_int8,LatticeIndex))
       neibs(10) = GetNeib_G(+1_int8,GetNeib_G(+3_int8,LatticeIndex))
       neibs(11) = GetNeib_G(+2_int8,GetNeib_G(+2_int8,LatticeIndex))
       neibs(12) = GetNeib_G(+2_int8,GetNeib_G(+3_int8,LatticeIndex))
       neibs(13) = GetNeib_G(+3_int8,GetNeib_G(+3_int8,LatticeIndex))
       neibs(14) = GetNeib_G(-1_int8,GetNeib_G(-1_int8,LatticeIndex))
       neibs(15) = GetNeib_G(-1_int8,GetNeib_G(-2_int8,LatticeIndex))
       neibs(16) = GetNeib_G(-1_int8,GetNeib_G(-3_int8,LatticeIndex))
       neibs(17) = GetNeib_G(-2_int8,GetNeib_G(-2_int8,LatticeIndex))
       neibs(18) = GetNeib_G(-2_int8,GetNeib_G(-3_int8,LatticeIndex))
       neibs(19) = GetNeib_G(-3_int8,GetNeib_G(-3_int8,LatticeIndex))
       neibs(20) = GetNeib_G(+1_int8,GetNeib_G(-2_int8,LatticeIndex))
       neibs(21) = GetNeib_G(+1_int8,GetNeib_G(-3_int8,LatticeIndex))
       neibs(22) = GetNeib_G(+2_int8,GetNeib_G(-3_int8,LatticeIndex))
       neibs(23) = GetNeib_G(-1_int8,GetNeib_G(+2_int8,LatticeIndex))
       neibs(24) = GetNeib_G(-1_int8,GetNeib_G(+3_int8,LatticeIndex))
       neibs(25) = GetNeib_G(-2_int8,GetNeib_G(+3_int8,LatticeIndex))

       !2. Reducing them (in case of very small lattice extensions indices will repeat)
       call RemoveDuplicates(neibs)
       
       !3. Computing number of non-zeros entries in diagonal and off-diagonal section
       dofMinRow = 1 + ndof*(SpatialRow-1)
       dofMaxRow = ndof*SpatialRow-1
       do concurrent(i=1:size(neibs))
          if(GetProc_G(neibs(i))==ThisProc()) then
             d_nnz(dofMinRow:dofMaxRow) = d_nnz(dofMinRow:dofMaxRow) + nDof
          else
             o_nnz(dofMinRow:dofMaxRow) = o_nnz(dofMinRow:dofMaxRow) + nDof
          end if
       end do
    end do
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
    
    deallocate(LocalSpatialPETScIndices)
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
  pure subroutine InitLocalSpatialPETScIndices(LocalSpatialPETScIndices)
    use lattice, only: GetMemorySize, GetLatticeIndex_M
    implicit none
    PetscInt, allocatable, intent(out) :: LocalSpatialPETScIndices(:)

    integer(int64) :: MemoryIndex
    
    allocate(LocalSpatialPETScIndices(GetMemorySize()))
    forall(MemoryIndex=1:GetMemorySize())
       LocalSpatialPETScIndices(MemoryIndex) = GetSpatialPETScIndex_G(GetLatticeIndex_M(MemoryIndex))
    end forall
  end subroutine InitLocalSpatialPETScIndices
    
  !>@brief Returns PETSc index from lattice index
  !!@returns PETSc index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetSpatialPETScIndex_G(LatticeIndex)
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

    GetSpatialPETScIndex_G = GetIndex_fromPosition(LocalPosition,LocalExtensions) &
         + proc * GetLocalLatticeSize() !&
         !- 1 ! PETSc-Index convention is C-like
  end function GetSpatialPETScIndex_G

  !>@brief Returns PETSc index from memory index
  !!@returns PETSc index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetSpatialPETScIndex_M(MemoryIndex)
    implicit none
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    GetSpatialPETScIndex_M = LocalSpatialPETScIndices(MemoryIndex)
  end function GetSpatialPETScIndex_M

  !>@brief Returns lattice index from PETSc index
  !!@returns lattice index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetLatticeIndex_PETSc(SpatialPETScIndex)
    use lattice, only: GetLatticeIndex_M
    implicit none
    !> PETSc index
    integer(int64), intent(in) :: SpatialPETScIndex
    GetLatticeIndex_PETSc = GetLatticeIndex_M(GetMemoryIndex_PETSc(SpatialPETScIndex))
  end function GetLatticeIndex_PETSc

  !>@brief Returns memory index from PETSc index
  !!@returns memory index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  pure integer(int64) function GetMemoryIndex_PETSc(SpatialPETScIndex)
    implicit none
    !> PETSc index
    integer(int64), intent(in) :: SpatialPETScIndex

    GetMemoryIndex_PETSc = FindLoc(&
         LocalSpatialPETScIndices,dim=1,& ! Where to look
         value = SpatialPETScIndex)       ! What to look for
  end function GetMemoryIndex_PETSc


  
end module nrqcd
