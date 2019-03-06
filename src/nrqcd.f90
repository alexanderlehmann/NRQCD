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

  !> Number of degrees of freedom per site
  integer(int8), parameter, public :: nDof = nSpins*nColours
  !> Number of Wilson coefficients
  integer(int8), parameter, public :: nWilsonCoeffs = 4

  !>@brief Heavy quark, described by NRQCD, containing quark- and antiquark-propagators
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 26.01.2019
  !!@version 1.0
  type NRQCDField
     !> Heavy-quark-propagator
     complex(real64), allocatable, public :: QuarkProp(:,:,:)

     !> Heavy antiquark-propagator
     complex(real64), allocatable, public :: AntiQprop(:,:,:)

   contains ! Member functions
     ! Allocation and deallocation
     procedure, public :: Allocate
     procedure, public :: Deallocate

     
  end type NRQCDField

contains ! Module procedures
  
  !>@brief Allocation of NRQCD heavy quark
  !!@details Allocates quark and antiquark propagator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 05.03.2019
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
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 05.03.2019
  !!@version 1.0
 pure subroutine Deallocate(HeavyField)
    implicit none
    !> NRQCD heavy field
    class(NRQCDField), intent(out) :: HeavyField
    if(allocated(HeavyField%QuarkProp))  deallocate(HeavyField%QuarkProp)
    if(allocated(HeavyField%AntiQProp))  deallocate(HeavyField%AntiQProp)    
  end subroutine Deallocate
end module nrqcd
