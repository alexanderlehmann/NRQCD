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
  !!@date 06.03.2019
  !!@version 1.0
  pure subroutine InitSinglePoint(Heavyfield,spin,colour,latticeindex)
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField
    !> Spin-component of created quark-antiquark-pair
    integer(int8),     intent(in)  :: spin
    !> Colour-component of created quark-antiquark-pair
    integer(int8),     intent(in)  :: colour
    !> Lattice index
    integer(int64),    intent(in)  :: latticeindex
    
    call HeavyField%Allocate
    HeavyField%QuarkProp = 0
    HeavyField%AntiQProp = 0

    
  end subroutine InitSinglePoint
end module nrqcd
