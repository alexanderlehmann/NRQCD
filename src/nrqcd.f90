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

  public NRQCDField, InitModule, FinalizeModule,S2CS,C2CS,GetLinkCS_M,GetLinkCS_G

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
     complex(fp), allocatable, public :: QuarkProp(:,:,:)

     !> Heavy antiquark-propagator
     complex(fp), allocatable, public :: AntiQprop(:,:,:)
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
     procedure, public :: InitSinglePointSingleDoF

     ! Norm
     procedure, public :: GetNorm_Quark
     procedure, public :: GetNorm_AntiQ

     ! Update methods
     generic, public :: Update => Update_CrankNicholson
     procedure, private :: Update_CrankNicholson

     ! Mesoncorrelators
     procedure, public :: GetMesoncorrelator_3S1_ZeroMomentum

     ! Access function
     procedure, public :: GetQuarkProp_G
     procedure, public :: GetAntiqProp_G
  end type NRQCDField


  ! ** PETSc variables **
  !> Preconditioner
  PC  :: PreCondMat
  !> PETSc system matrix
  Mat :: SystMat, SystMat_herm, SystMat_precond
  !> Solution vector
  Vec :: PETScX
  !> RHS vector
  Vec :: PETScRHS
  !> Prec-conditioned RHS vector
  Vec :: PETScRHS_precond
  !> Solver context
  KSP :: ksp
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
  impure subroutine InitModule(Mass,WilsonCoeffs,StepWidth,kspTol)
    use, intrinsic :: iso_fortran_env
    use mpiinterface, only: mpistop
    implicit none
    !> Bare heavy quark mass
    real(fp), intent(in) :: Mass
    !> Wilson coefficients
    complex(fp), intent(in) :: WilsonCoeffs(nWilsonCoefficients)
    !> Step width in units of \f$a_0\f$
    real(fp), intent(in) :: StepWidth
    !> Tolerance for iterative solver
    real(fp), intent(in) :: kspTol
    
    if(isInitialised) then
       call MPISTOP('Error in init of '//modulename//': already initialised.')
    else
       
       call CheckObligatoryInitialisations

       ! Initialise PETSc-solver
       call InitPETScSolver(Mass,WilsonCoeffs,StepWidth,kspTol)
       
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
  end subroutine assign
  
  !>@brief Constructor
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 13.03.2019
  !!@version 1.0
  impure function NRQCDField_Constructor()
    use mpiinterface, only: ThisProc, NumProcs,SyncAll,mpistop

    implicit none
    
    type(NRQCDField) :: NRQCDField_Constructor

    PetscMPIInt :: Petscerr

    integer(int8) :: i
    
    call NRQCDField_Constructor%Allocate
  end function NRQCDField_Constructor

  
  !>@brief Destructor
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 15.03.2019
  !!@version 1.0
  impure subroutine Destructor(object)
    implicit none
    class(NRQCDField), intent(inout) :: object
    
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
  !! at a single spatial point with a normalized unit matrix in colour x spin -space
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.04.2019
  !!@version 1.0
  impure subroutine InitSinglePoint(Heavyfield,LatticeIndex_quark,LatticeIndex_antiq)
    use mpiinterface, only: ThisProc
    use lattice, only: GetProc_G, GetMemoryIndex
    use matrixoperations, only: GetUnitMatrix
    use precision, only: fp
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField
    !> Lattice index of created quark
    integer(int64),    intent(in)  :: LatticeIndex_quark,LatticeIndex_antiq
    
    HeavyField = NRQCDField()
    
    HeavyField%QuarkProp = 0
    HeavyField%AntiQProp = 0

    if(ThisProc()==GetProc_G(LatticeIndex_quark)) then
       HeavyField%QuarkProp(&
            :,:,GetMemoryIndex(LatticeIndex_quark)) = GetUnitMatrix(nDoF)/sqrt(real(nDoF,fp))
    end if
    if(ThisProc()==GetProc_G(LatticeIndex_antiq)) then
       HeavyField%AntiQProp(&
            :,:,GetMemoryIndex(LatticeIndex_antiq)) = GetUnitMatrix(nDoF)/sqrt(real(nDoF,fp))
    end if

    call HeavyField%CommunicateBoundary
  end subroutine InitSinglePoint

  !>@brief Initialises the propagator and anti-propagator
  !! at a single spatial point for a single colour- and spin-combination
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine InitSinglePointSingleDoF(Heavyfield,spin_quark,colour_quark,LatticeIndex_quark,spin_antiq,colour_antiq,LatticeIndex_antiq)
    use mpiinterface, only: ThisProc
    use lattice, only: GetProc_G, GetMemoryIndex
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField
    !> Spin-component of created quark
    integer(int8),     intent(in)  :: spin_quark,spin_antiq
    !> Colour-component of created quark
    integer(int8),     intent(in)  :: colour_quark,colour_antiq
    !> Lattice index of created quark
    integer(int64),    intent(in)  :: LatticeIndex_quark,LatticeIndex_antiq
    
    !call HeavyField%Allocate
    HeavyField = NRQCDField()
    
    HeavyField%QuarkProp = 0
    HeavyField%AntiQProp = 0

    if(ThisProc()==GetProc_G(LatticeIndex_quark)) then
       HeavyField%QuarkProp(&
            GetSpinColourIndex(Spin_quark,Colour_quark),&
            GetSpinColourIndex(Spin_quark,Colour_quark),&
            GetMemoryIndex(LatticeIndex_quark)) = 1
    end if
    if(ThisProc()==GetProc_G(LatticeIndex_antiq)) then
       HeavyField%AntiQProp(&
            GetSpinColourIndex(Spin_antiq,Colour_antiq),&
            GetSpinColourIndex(Spin_antiq,Colour_antiq),&
            GetMemoryIndex(LatticeIndex_antiq)) = 1
    end if

    call HeavyField%CommunicateBoundary
  end subroutine InitSinglePointSingleDoF

  !>@brief Access to quark propagator
  !!@returns Quark propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 02.04.2019
  !!@version 1.0
  pure function GetQuarkProp_G(HeavyField,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(in) :: HeavyField
    !> Lattice index
    integer(int64),    intent(in)  :: LatticeIndex
    !> Quark propagator
    complex(fp), dimension(ndof,ndof) :: GetQuarkProp_G

    GetQuarkProp_G = HeavyField%QuarkProp(:,:,GetMemoryIndex(LatticeIndex))
  end function GetQuarkProp_G

  !>@brief Access to antiquark propagator
  !!@returns Antiquark propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 02.04.2019
  !!@version 1.0
  pure function GetAntiQProp_G(HeavyField,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(in) :: HeavyField
    !> Lattice index
    integer(int64),    intent(in)  :: LatticeIndex
    !> Quark propagator
    complex(fp), dimension(ndof,ndof) :: GetAntiQProp_G

    GetAntiQProp_G = HeavyField%AntiQProp(:,:,GetMemoryIndex(LatticeIndex))
  end function GetAntiQProp_G

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
  impure subroutine InitPETScSolver(Mass,WilsonCoeffs,StepWidth,kspTol)
    use mpiinterface, only: intmpi,mpistop,thisproc,numprocs
    use lattice, only: GetLatticeSize, GetLocalLatticeSize, GetMemorySize, GetNeib_G, GetProc_G
    use arrayoperations, only: RemoveDuplicates
    implicit none
    !> Bare heavy quark mass
    real(fp), intent(in) :: Mass
    !> Wilson coefficients
    complex(fp), intent(in) :: WilsonCoeffs(nWilsonCoefficients)
    !> Step width in units of \f$a_0\f$
    real(fp), intent(in) :: StepWidth
    !> Tolerance for iterative solver
    real(fp), intent(in) :: kspTol
    
    PetscErrorCode :: Petscerr
    PetscReal :: Tol

    PetscInt :: SpatialRow, SpatialCol
    integer(int64) :: LatticeIndex, i
    
    integer(int64), allocatable :: neibs(:)
    
    integer(int64) :: dofMinRow,dofMaxRow

    ! Dummy variables for first matrix building
    type(SU3GaugeConfiguration) :: DummyConf

    PetscScalar :: TestScalar
    
    !call PETScInitialize(PETSC_NULL_CHARACTER,Petscerr)

    !if(Petscerr /= 0) then
    !   call MPIStop(&
    !        errormessage = 'Error in initialization of '//&
    !        modulename//': Initialisation of PETSc failed.',&
    !        errorcode = Petscerr)
    !end if

    TestScalar = cmplx(1,1,fp)
    TestScalar = conjg(TestScalar)

    call InitLocalSpatialPETScIndices(LocalSpatialPETScIndices)
    
    SystMatSize = GetLatticeSize()*ndof
    SystMatSize = SystMatSize
    LocalSystMatSize = GetLocalLatticeSize()*ndof

    LocalSpatialMinRow = 1+ ThisProc()   *GetLocalLatticeSize()
    LocalSpatialMaxRow =   (ThisProc()+1)*GetLocalLatticeSize()
    LocalMinRow = 1 + (LocalSpatialMinRow-1)*nDof
    LocalMaxRow = (ThisProc()+1)*GetLocalLatticeSize()*nDof
    
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
       dofMaxRow = ndof*SpatialRow
       do i=1,size(neibs)
          if(GetProc_G(neibs(i))==ThisProc()) then
             d_nnz(dofMinRow:dofMaxRow) = d_nnz(dofMinRow:dofMaxRow) + nDof
          else
             o_nnz(dofMinRow:dofMaxRow) = o_nnz(dofMinRow:dofMaxRow) + nDof
          end if
       end do
    end do

    ! PETSc variables for linear solver
    ! 1. KSP-solver object
    call KSPCreate(PETSC_COMM_WORLD,ksp,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': KSPCreate failed.',&
            errorcode = Petscerr)
    end if
    !call KSPSetType(ksp,KSPGMRES,Petscerr)
    
    ! 2. System matrix
    call MatCreateAIJ(PETSC_COMM_WORLD,&
         LocalSystMatSize,LocalSystMatSize,&
         SystMatSize,SystMatSize,&
         PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,&
         SystMat,&
         Petscerr)
    call MatCreateAIJ(PETSC_COMM_WORLD,&
         LocalSystMatSize,LocalSystMatSize,&
         SystMatSize,SystMatSize,&
         PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,&
         SystMat_herm,&
         Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': MatCreateAIJ failed.',&
            errorcode = Petscerr)
    end if
    
    ! Tell PETSc that the matrices are structurally (non-zero-pattern-wise) symmetric
    call MatSetOption(SystMat,        MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE,PETScErr)
    call MatSetOption(SystMat_herm,   MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE,PETScErr)
    ! Tell PETSc to keep the non-zero-pattern, thus never "compressing" it
    call MatSetOption(SystMat,        MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,PETScErr)
    call MatSetOption(SystMat_herm,   MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,PETScErr)
    ! Tell PETSc to keep the symmetry pattern for ever
    call MatSetOption(SystMat,        MAT_SYMMETRY_ETERNAL,PETSC_TRUE,PETScErr)
    call MatSetOption(SystMat_herm,   MAT_SYMMETRY_ETERNAL,PETSC_TRUE,PETScErr)

    call DummyConf%ColdInit
    call BuildSystemMatrix(SystMat_herm,DummyConf,1._fp,WilsonCoeffs,StepWidth)
    call BuildSystemMatrix(SystMat,     DummyConf,1._fp,WilsonCoeffs,StepWidth)
    
    call MatMatMultSymbolic(SystMat_herm,SystMat,PETSC_DEFAULT_REAL,SystMat_precond,PETScErr)
    
    call MatSetOption(SystMat_precond,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE,PETScErr)
    call MatSetOption(SystMat_precond,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,PETScErr)
    call MatSetOption(SystMat_precond,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,PETScErr)
    
    ! 3. Pre-conditioner
    call PCCreate(PETSC_COMM_WORLD,PreCondMat,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': PCCreate failed.',&
            errorcode = Petscerr)
    end if
    call PCSetType(PreCondMat,PCJACOBI,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': PCSetType to PCJACOBI failed.',&
            errorcode = Petscerr)
    end if
    
    ! 4. Set Tolerance
    Tol = kspTol
    call KSPSetTolerances(ksp,Tol,&
         PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': KSPSetTOLERANCES failed.',&
            errorcode = Petscerr)
    end if

    ! 5. Solution vector
    call VecCreateMPI(PETSC_COMM_WORLD,&
         LocalSystMatSize,SystMatsize,&
         PETScX,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': VecCreateMPI for x failed.',&
            errorcode = Petscerr)
    end if

    ! 6. RHS vector
    call VecCreateMPI(PETSC_COMM_WORLD,&
         LocalSystMatSize,SystMatsize,&
         PETScRHS,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': VecCreateMPI for rhs failed.',&
            errorcode = Petscerr)
    end if
    
    ! 7. Preconditioned RHS vector
    call VecCreateMPI(PETSC_COMM_WORLD,&
         LocalSystMatSize,SystMatsize,&
         PETScRHS_precond,Petscerr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in initialization of '//&
            modulename//': VecCreateMPI for preconditioned rhs failed.',&
            errorcode = Petscerr)
    end if

    call DummyConf%Deallocate
  end subroutine InitPETScSolver

  !>@brief Finalizes PETSc-solver
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 11.03.2019
  !!@version 1.0
  impure subroutine FinalizePETScSolver
    use mpiinterface, only: mpistop
    implicit none

    PetscErrorCode Petscerr
    
    call VecDestroy(PETScX,Petscerr)
    call VecDestroy(PETScRHS,Petscerr)
    call VecDestroy(PETScRHS_PreCond,Petscerr)

    call MatDestroy(SystMat,Petscerr)
    call MatDestroy(SystMat_herm,Petscerr)
    call MatDestroy(SystMat_PreCond,Petscerr)
    call KSPDestroy(ksp,Petscerr)
    
    deallocate(LocalSpatialPETScIndices)
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
  pure PetscInt function GetSpatialPETScIndex_G(LatticeIndex)
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
    PetscInt, intent(in) :: SpatialPETScIndex
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
    PetscInt, intent(in) :: SpatialPETScIndex

    GetMemoryIndex_PETSc = FindLoc(&
         LocalSpatialPETScIndices,dim=1,& ! Where to look
         value = SpatialPETScIndex)       ! What to look for
  end function GetMemoryIndex_PETSc

  !>@brief Adds submatrix to system matrix at spatial PETSc row and column indices
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.03.2019
  !!@version 1.0
  impure subroutine AddMatrixToSystem(SpatialRow,SpatialCol,SubMat,SystMat)
    use mpiinterface, only: MPIStop
    implicit none
    !> Spatial row index
    PetscInt, intent(in) :: SpatialRow
    !> Spatial column index
    PetscInt, intent(in) :: SpatialCol
    !> Submatrix
    complex(fp), intent(in) :: submat(nDof,nDof)
    !> System matrix
    Mat, intent(inout) :: SystMat
    
    PetscErrorCode :: PETScErr

    PetscInt :: Cols(nDof), Rows(nDof)
    PetscInt, parameter :: nRow=nDof, nCol=nDof
    PetscScalar :: v(nRow*nCol)

    Rows = [1_int8:nDof] + nDof*(SpatialRow-1) - 1 ! -1 for C indexing
    Cols = [1_int8:nDof] + nDof*(SpatialCol-1) - 1 ! -1 for C indexing
    v = reshape(SubMat,[nDof**2])
    
    call MatSetValues(SystMat,nRow,Rows,nCol,Cols,v,ADD_VALUES,PETScErr)
    if(PETScErr /= 0) then
       call MPIStop(&
            errormessage = 'Error in AddMatrixToSystem of '//&
            modulename//': MatSetValues failed.',&
            errorcode = PETScErr)
    end if
  end subroutine AddMatrixToSystem

  !>@brief Resets system matrix to zero
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.03.2019
  !!@version 1.0
  impure subroutine ResetSystemMatrix(SystMat)
    use mpiinterface, only: MPIStop
    implicit none
    !> System matrix
    Mat, intent(inout) :: SystMat

    PetscErrorCode :: PETScErr
    call MatZeroEntries(SystMat,PETScErr) !retains the non-zero structure
    if(PETScErr /= 0) then
       call MPIStop(&
            errormessage = 'Error in ResetSystemMatrix of '//&
            modulename//': MatZeroEntries failed.',&
            errorcode = PETScErr)
    end if
  end subroutine ResetSystemMatrix
  
  !>@brief Performs update step using Crank-Nicholson scheme
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.03.2019
  !!@version 1.0
  impure subroutine Update_CrankNicholson(HeavyField,GaugeConf,Mass,WilsonCoeffs,StepWidth)
    use mpiinterface, only: mpistop, ThisProc

use lattice
    use mpiinterface, only: numprocs, intmpi, syncall
    implicit none
    !> NRQCD heavy field
    class(NRQCDField),           intent(inout) :: HeavyField
    !> Gauge configuration
    type(SU3GaugeConfiguration), intent(in)    :: GaugeConf
    !> Bare mass of heavy quark field
    real(fp),                    intent(in)    :: Mass
    !> Wilson coefficients
    complex(fp),                 intent(in)    :: WilsonCoeffs(nWilsonCoefficients)
    !> Step width in units of \f$a_0\f$
    real(fp), optional,          intent(in)    :: StepWidth

    real(fp) :: dt

    complex(fp) :: wilsonCoeffs_conjg(nWilsonCoefficients)
    real(fp)    :: mass_negative

    integer(int8) :: i
    PetscErrorCode :: PETScErr

    integer(intmpi) :: proc
    integer(int64) :: row

    ! ..--** START: Optional parameters **--..
    if(present(StepWidth)) then
       dt = StepWidth
    else
       ! Default value
       dt = 1._real64 ! in units of a0
    end if
    ! ..--**  END : Optional parameters **--..
    
    ! Quark propagator
    ! Calling the buildsystemmatrix-routine on all matrices ensures that they have the same
    ! non-zero-pattern
    call BuildSystemMatrix(SystMat_herm,GaugeConf,+Mass,WilsonCoeffs,dt)
    call BuildSystemMatrix(SystMat,GaugeConf,+Mass,WilsonCoeffs,dt)
    
    ! ... and hermitian conjugate
    call MatHermitianTranspose(SystMat,MAT_INPLACE_MATRIX,SystMat,PETScErr)
    
    ! ... perform custom preconditioning step
    !call MatMatMultNumeric(SystMat_herm,SystMat,SystMat_PreCond,PETScErr)
    !call KSPSetOperators(ksp,SystMat_precond,SystMat_precond,PETScErr)
    
    call KSPSetOperators(ksp,SystMat,SystMat,PETScErr)
    
    call KSPSetUp(ksp,Petscerr)

    do i=1,nDof
       ! Load from configuration
       call LoadPropagatorIntoPETSc(HeavyField%QuarkProp,i,PETScX)

       ! Perform half explicit step
       call MatMult(SystMat_herm,PETScX,PETScRHS,PETScErr)

       call LoadPropagatorFromPETSc(HeavyField%QuarkProp,i,PETScRHS)

       !**: START preconditioning
       ! Perform preconditioning step on RHS
       !call MatMult(SystMat_herm,PETScRHS,PETScRHS_precond,PETScErr)       
       ! Perform half implicit step
       !call KSPSolve(ksp,PETScRHS_precond,PETScX,PETScErr)
       !**: END preconditioning

       call KSPSolve(ksp,PETScRHS,PETScX,PETScErr)

       ! Load into configuration
       call LoadPropagatorFromPETSc(HeavyField%QuarkProp,i,PETScX)
    end do

    ! Anti quark propagator
    wilsonCoeffs_conjg = conjg(WilsonCoeffs)
    mass_negative = -Mass
    ! Calling the buildsystemmatrix-routine on all matrices ensures that they have the same
    ! non-zero-pattern
    call BuildSystemMatrix(SystMat_herm,GaugeConf,mass_negative,WilsonCoeffs_conjg,dt)
    call BuildSystemMatrix(SystMat,GaugeConf,mass_negative,WilsonCoeffs_conjg,dt)
    
    ! ... and hermitian conjugate
    call MatHermitianTranspose(SystMat,MAT_INPLACE_MATRIX,SystMat,PETScErr)
    
    ! ... perform custom preconditioning step
    !call MatMatMultNumeric(SystMat_herm,SystMat,SystMat_PreCond,PETScErr)
    !call KSPSetOperators(ksp,SystMat_precond,SystMat_precond,PETScErr)
    
    call KSPSetOperators(ksp,SystMat,SystMat,PETScErr)
    
    call KSPSetUp(ksp,Petscerr)
    
    do i=1,nDof
       ! Load from configuration
       call LoadPropagatorIntoPETSc(HeavyField%AntiQProp,i,PETScX)
       
       ! Perform half explicit step
       call MatMult(SystMat_herm,PETScX,PETScRHS,PETScErr)
       
       ! Perform preconditioning step on RHS
       !call MatMult(SystMat_herm,PETScRHS,PETScRHS_precond,PETScErr)
       
       ! Perform half implicit step
       !call KSPSolve(ksp,PETScRHS_precond,PETScX,PETScErr)

       call KSPSolve(ksp,PETScRHS,PETScX,PETScErr)
       
       ! Load into configuration
       call LoadPropagatorFromPETSc(HeavyField%AntiQProp,i,PETScX)
    end do

    call HeavyField%CommunicateBoundary
  end subroutine Update_CrankNicholson

  !>@brief Builds system matrix
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.03.2019
  !!@version 1.0
  impure subroutine BuildSystemMatrix(SystMat,GaugeConf,Mass,WilsonCoeffs,StepWidth)
    use mpiinterface, only: SyncAll
    use matrixoperations, only: GetUnitMatrix, GetKronProd
    use lattice, only: nDim,GetLatticeSpacing,GetNeib_G

    use mpiinterface, only: ThisProc, MPIStop
    implicit none
    !> System matrix
    Mat, intent(inout) :: SystMat
    !> Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeConf
    !> Bare heavy field mass
    real(fp), intent(in)    :: Mass
    !> Wilson coefficients
    complex(fp), intent(in) :: WilsonCoeffs(nWilsonCoefficients)
    !> Step width in units of \f$a_0\f$
    real(fp),    intent(in) :: StepWidth

    PetscInt :: SpatialRow, SpatialCol
    PetscErrorCode :: PETScErr

    integer(int8) :: i,j,k, LeviCivita
    integer(int64) :: neib_p, neib_m, latticeindex
    complex(fp), dimension(nDof,nDof) :: submatrix, link_CS
    complex(fp), dimension(nColours,nColours) :: link, link_p, link_m, efield, efield_p, efield_m, derivEfield
    
    call ResetSystemMatrix(SystMat)

    ! O(dt^0)
    ! Unit matrix
    submatrix = GetUnitMatrix(nDof)
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       call AddMatrixToSystem(SpatialRow,SpatialRow,submatrix,SystMat)
    end do

    ! O(dt^1)
    ! .. -c1/2m D^2
    !goto 1 ! Leave out diffusion term which corresponds to p²/(2m)
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       do i=1,nDim
          ! + delta_{x,y}/(ai**2)/m
          SpatialCol = SpatialRow
          
          submatrix = + cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *WilsonCoeffs(1)/mass/GetLatticeSpacing(i)**2 &
               *GetUnitMatrix(nDof)
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)
          
          ! - delta_{x+î,y}U(x,i)/(ai**2)/2m
          neib_p = GetNeib_G(+i,GetLatticeIndex_PETSc(SpatialRow))
          SpatialCol = GetSpatialPETScIndex_G(neib_p)
          link_CS = GetLinkCS_G(GaugeConf,i,GetLatticeIndex_PETSc(SpatialRow))
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *(-WilsonCoeffs(1))/2/mass/GetLatticeSpacing(i)**2 &
               *Link_CS
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)
          
          ! - delta_{x-î,y}U(x-î,i)†/(ai**2)/2m
          neib_m = GetNeib_G(-i,GetLatticeIndex_PETSc(SpatialRow))
          SpatialCol = GetSpatialPETScIndex_G(neib_m)
          link_CS = GetLinkCS_G(GaugeConf,i,neib_m)
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *(-WilsonCoeffs(1))/2/mass/GetLatticeSpacing(i)**2 &
               *conjg(transpose(Link_CS))
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)
          
       end do
    end do
    
    ! .. -c3/8m^2.[D_i,g.at.E_i]
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       do i=1,nDim
          ! + delta_{x,y}/2m ( F(x,jk)-F(x,kj) )/2
          SpatialCol = SpatialRow
          LatticeIndex = GetLatticeIndex_PETSc(SpatialCol)

          efield_p = GaugeConf%GetElectricField(i,GetNeib_G(+i,LatticeIndex))
          efield_m = GaugeConf%GetElectricField(i,GetNeib_G(-i,LatticeIndex))
          link     = GaugeConf%GetLink_G(i,LatticeIndex)
          link_m   = GaugeConf%GetLink_G(i,GetNeib_G(-i,LatticeIndex))

          DerivEfield = (&
               + matmul(matmul(                 link,    efield_p),conjg(transpose(link  ))) &
               - matmul(matmul( conjg(transpose(link_m)),efield_m),                link_m)   &
               )/(2*GetLatticeSpacing(i))
          
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *(-WilsonCoeffs(3))/8/mass**2 &
               *C2CS(DerivEfield)
          
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)
       end do
    end do

    ! .. -c2.g/2m B.sigma
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       do i=1,nDim
          ! + delta_{x,y}/2m ( F(x,jk)-F(x,kj) )/2
          SpatialCol = SpatialRow
          LatticeIndex = GetLatticeIndex_PETSc(SpatialCol)
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *(-WilsonCoeffs(2))/2/mass &
               *GetMagneticField_CS_G(GaugeConf,i,LatticeIndex)
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)
       end do
    end do
    
    ! .. -c4.i/8m^2.eps_ijk sigma_k x {D_i,g.at.E_j}
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       do i=1,nDim
          LatticeIndex = GetLatticeIndex_PETSc(SpatialRow)
          neib_p = GetNeib_G(+i,LatticeIndex)
          neib_m = GetNeib_G(-i,LatticeIndex)
          
          ! Indices of ε:
          ! positive permutation (i,j,k)
          j = mod(i,  ndim)+1
          k = mod(i+1,ndim)+1
          LeviCivita=+1
          
          ! - ì ε_{ijk} delta_{x+î,y}/(ai**2)/m
          SpatialCol = GetSpatialPETScIndex_G(neib_p)

          efield   = GaugeConf%GetElectricField(j,LatticeIndex)
          efield_p = GaugeConf%GetElectricField(j,neib_p)
          link     = GaugeConf%GetLink_G(i,LatticeIndex)

          DerivEfield = ( matmul(efield,link) + matmul(link,efield_p) )/2/GetLatticeSpacing(i)
          
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *LeviCivita&
               *(-WilsonCoeffs(4))/8/mass**2*cmplx(0._fp,1._fp,fp) &
               *GetKronProd(DerivEfield,SU2Generators(:,:,k))

          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)

          ! + ì ε_{ijk} delta_{x+î,y}/(ai**2)/m
          SpatialCol = GetSpatialPETScIndex_G(neib_m)

          efield   = GaugeConf%GetElectricField(j,LatticeIndex)
          efield_m = GaugeConf%GetElectricField(j,neib_m)
          link_m   = conjg(transpose(GaugeConf%GetLink_G(i,neib_m)))

          DerivEfield = ( matmul(efield,link_m) + matmul(link_m,efield_m) )/2/GetLatticeSpacing(i)
          
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *LeviCivita&
               *(+WilsonCoeffs(4))/8/mass**2*cmplx(0._fp,1._fp,fp) &
               *GetKronProd(DerivEfield,SU2Generators(:,:,k))
          
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)

          ! Indices of ε:
          ! negative permutation (i,j,k)
          j = mod(i+1,ndim)+1
          k = mod(i,  ndim)+1
          LeviCivita=-1
          
          ! - ì ε_{ijk} delta_{x+î,y}/(ai**2)/m
          SpatialCol = GetSpatialPETScIndex_G(neib_p)

          efield   = GaugeConf%GetElectricField(j,LatticeIndex)
          efield_p = GaugeConf%GetElectricField(j,neib_p)
          link     = GaugeConf%GetLink_G(i,LatticeIndex)

          DerivEfield = ( matmul(efield,link) + matmul(link,efield_p) )/2/GetLatticeSpacing(i)
          
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *LeviCivita&
               *(-WilsonCoeffs(4))/8/mass**2*cmplx(0._fp,1._fp,fp) &
               *GetKronProd(DerivEfield,SU2Generators(:,:,k))

          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)

          ! + ì ε_{ijk} delta_{x+î,y}/(ai**2)/m
          SpatialCol = GetSpatialPETScIndex_G(neib_m)

          efield   = GaugeConf%GetElectricField(j,LatticeIndex)
          efield_m = GaugeConf%GetElectricField(j,neib_m)
          link_m   = conjg(transpose(GaugeConf%GetLink_G(i,neib_m)))

          DerivEfield = ( matmul(efield,link_m) + matmul(link_m,efield_m) )/2/GetLatticeSpacing(i)
          
          submatrix = cmplx(0._fp,-StepWidth/2*GetLatticeSpacing(0_int8),fp)&
               *LeviCivita&
               *(+WilsonCoeffs(4))/8/mass**2*cmplx(0._fp,1._fp,fp) &
               *GetKronProd(DerivEfield,SU2Generators(:,:,k))
          
          call AddMatrixToSystem(SpatialRow,SpatialCol,submatrix,SystMat)
       end do
    end do
    
    call SyncAll
    call MatAssemblyBegin(SystMat,MAT_FINAL_ASSEMBLY,PETScErr)
    call MatAssemblyEnd(SystMat,MAT_FINAL_ASSEMBLY,PETScErr)
    call SyncAll
  end subroutine BuildSystemMatrix

  !>@brief Loads propagator into PETSc-vectors
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.03.2019
  !!@version 1.0
  impure subroutine LoadPropagatorIntoPETSc(Propagator,PropagatorCol,PETScVec)
    use mpiinterface, only: SyncAll, MPIStop
    use lattice, only: GetLocalLatticeSize
    implicit none
    !> Propagator
    complex(fp),   intent(in) :: Propagator(:,:,:)
    !> Column of propagator to be read
    integer(int8), intent(in) :: PropagatorCol
    !> PETSc vector
    Vec, intent(inout) :: PETScVec

    PetscErrorCode :: PETScErr

    PetscScalar :: v(nDof)
    PetscInt, parameter :: nRows=nDof
    PetscInt :: rows(nRows),SpatialRow
    
    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       rows = [1:nDof] + (SpatialRow-1)*nDof - 1

       v = Propagator(:,PropagatorCol,GetMemoryIndex_PETSc(SpatialRow))
       call VecSetValues(PETScVec,nRows,rows,v,INSERT_VALUES,PETScErr)
       if(Petscerr /= 0) then
          call MPIStop(&
               errormessage = 'Error in LoadPropagatorIntoPETSc of '//&
               modulename//': VecSetValues failed.',&
               errorcode = Petscerr)
       end if
    end do

    call SyncAll
    call VecAssemblyBegin(PETScVec,PETScErr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in LoadPropagatorIntoPETSc of '//&
            modulename//': VecAssemblyBegin failed.',&
            errorcode = Petscerr)
    end if
    
    call VecAssemblyEnd(PETScVec,PETScErr)
    if(Petscerr /= 0) then
       call MPIStop(&
            errormessage = 'Error in LoadPropagatorIntoPETSc of '//&
            modulename//': VecAssemblyEnd failed.',&
            errorcode = Petscerr)
    end if
  end subroutine LoadPropagatorIntoPETSc

  !>@brief Loads propagator from PETSc-vectors
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 18.03.2019
  !!@version 1.0
  impure subroutine LoadPropagatorFromPETSc(Propagator,PropagatorCol,PETScVec)
    use mpiinterface, only: MPIStop
    implicit none
    !> Propagator
    complex(fp),   intent(out) :: Propagator(:,:,:)
    !> Column of propagator to be read
    integer(int8), intent(in)  :: PropagatorCol
    !> PETSc vector
    Vec, intent(in) :: PETScVec

    PetscErrorCode :: PETScErr

    PetscScalar :: v(nDof)
    PetscInt, parameter :: nRows=nDof
    PetscInt :: rows(nRows), SpatialRow

    do SpatialRow=LocalSpatialMinRow,LocalSpatialMaxRow
       rows = [1:nDof] + (SpatialRow-1)*nDof - 1

       call VecGetValues(PETScVec,nRows,rows,v,PETScErr)
       if(Petscerr /= 0) then
          call MPIStop(&
               errormessage = 'Error in LoadPropagatorFromPETSc of '//&
               modulename//': VecGetValues failed.',&
               errorcode = Petscerr)
       end if
       
       Propagator(:,PropagatorCol,GetMemoryIndex_PETSc(SpatialRow)) = v
    end do
  end subroutine LoadPropagatorFromPETSc

  !>@brief Computes \f$^3\text{S}_1\f$-Quarkonium-correlator at \f$\vec{p}=\vec{0}\f$
  !!@returns \f$^3\text{S}_1\f$-Quarkonium-correlator at \f$\vec{p}=\vec{0}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 19.03.2019
  !!@version 1.0
  impure complex(fp) function GetMesoncorrelator_3S1_ZeroMomentum(HeavyField)
    use mpi
    use mpiinterface, only: intmpi, ThisProc, GetComplexSendType
    use lattice, only: nDim, GetProc_M, GetMemorySize
    use matrixoperations, only: GetTrace
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(in) :: HeavyField

    integer(intmpi) :: mpierr
    complex(fp) :: PauliMatrix_CS(nDof,nDof), Correlator(nDof,nDof)
    integer(int64) :: MemoryIndex
    integer(int8) :: i
    complex(fp) :: LocalValue
    
    LocalValue = 0
    do concurrent(MemoryIndex=1:GetMemorySize(),GetProc_M(MemoryIndex)==ThisProc())
       Correlator = 0
       do concurrent(i=1:nDim)
          PauliMatrix_CS = S2CS(SU2Generators(:,:,i))

          !Correlator = Correlator &
          !     + cmplx(0,1,fp)*matmul(matmul(matmul(&
          !     PauliMatrix_CS,&
          !     HeavyField%QuarkProp(:,:,MemoryIndex)),&
          !     PauliMatrix_CS),&
          !     conjg(HeavyField%AntiQProp(:,:,MemoryIndex)))

          
          Correlator = Correlator &
               + cmplx(0,1,fp)*matmul(matmul(matmul(&
               PauliMatrix_CS,&
               HeavyField%QuarkProp(:,:,MemoryIndex)),&
               conjg(PauliMatrix_CS)),&
               conjg(transpose(HeavyField%AntiQProp(:,:,MemoryIndex))))
       end do
       !LocalValue = LocalValue &
       !     + GetTrace(Correlator)*nspins**2/ncolours

       
       LocalValue = LocalValue + GetTrace(Correlator)/(2*ncolours)
    end do

    call MPI_ALLREDUCE(&
         LocalValue,&
         GetMesoncorrelator_3S1_ZeroMomentum,&
         1_intmpi,&
         GetComplexSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)
  end function GetMesoncorrelator_3S1_ZeroMomentum

  !>@brief Computes chromo-magnetic field
  !!@returns chromo-magnetic field
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.03.2019
  !!@version 1.0
  pure function GetMagneticField_CS_G(GaugeConf,i,LatticeIndex)
    use lattice, only: nDim
    implicit none
    !> Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeConf
    !> Spatial direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex

    complex(fp), dimension(nDoF,nDoF) :: GetMagneticField_CS_G

    complex(fp), dimension(nColours,nColours) :: MagneticField
    
    ! Spatial directions as positive permutation (i,j,k)
    integer(int8) :: j,k
    
    j = modulo(i,       nDim)+1_int8
    k = modulo(i+1_int8,nDim)+1_int8

    MagneticField = ( GaugeConf%GetFieldStrengthTensor_G(j,k,LatticeIndex)&
         + GaugeConf%GetFieldStrengthTensor_G(k,j,LatticeIndex) ) / 2
    
    GetMagneticField_CS_G = C2CS(MagneticField)
  end function GetMagneticField_CS_G
end module nrqcd
