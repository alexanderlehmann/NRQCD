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

  implicit none

  PRIVATE

  public NRQCDField

  !> Module name
  character(len=5), parameter, public ::  modulename='nrqcd'

  !> Number of degrees of freedom per site
  integer(int8), parameter, public :: nDof = nSpins*nColours
  !> Number of Wilson coefficients
  integer(int8), parameter, public :: nWilsonCoefficients = 4

  !>@brief Heavy quark, described by NRQCD, containing quark- and antiquark-propagators
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 06.03.2019
  !!@version 1.0
  type NRQCDField
     !> Heavy-quark-propagator
     complex(fp), allocatable, public :: QuarkProp(:,:,:)

     !> Heavy antiquark-propagator
     complex(fp), allocatable, public :: AntiQprop(:,:,:)

   contains ! Member functions
     ! Allocation and deallocation
     procedure, public :: Allocate
     procedure, public :: Deallocate

     ! Boundary-communication
     procedure, public :: CommunicateBoundary
     
     ! Initialisation routines
     procedure, public :: InitSinglePoint

     ! Norm
     procedure, public :: GetNorm_Quark
     procedure, public :: GetNorm_AntiQ
     
  end type NRQCDField

contains ! Module procedures
  
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
    
    call HeavyField%Allocate
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
end module nrqcd
