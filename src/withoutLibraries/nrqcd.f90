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
       nSpins        => nSUN, &
       nSpinGenerators => nGen
  
  implicit none

  PRIVATE

  public InitModule, FinalizeModule,S2CS,C2CS,GetLinkCS_M,GetLinkCS_G

  !> Module name
  character(len=5), parameter, public ::  modulename='nrqcd'
  
  !> Contains information, whether module is initialised
  logical :: IsInitialised = .false.

  !> Number of degrees of freedom per site
  integer(int8), parameter, public :: nDof = nSpins*nColours
  !> Number of Wilson coefficients
  integer(int8), parameter, public :: nWilsonCoefficients = 4

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
       
       ! DONE
       IsInitialised = .TRUE.
    end if
  end subroutine InitModule
  
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

       IsInitialised = .FALSE.
    else
       call MPISTOP('Error in finalization of '//modulename//': is not initialised.')
    end if
  end subroutine FinalizeModule

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
  
end module nrqcd
