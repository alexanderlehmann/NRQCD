!----------------------------------------------------------------------
! RTLQCD, Real-Time-Lattice-QCD Simulation of Gauge Fields
!----------------------------------------------------------------------
!
!>@brief SU(3)-gauge-link-configuration in temporal gauge
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 22.01.2019
!!@version 1.0
!----------------------------------------------------------------------
module gaugeconfiguration_su3
  use, intrinsic :: iso_fortran_env
  use precision, only: fp


  use mpiinterface

  
  use su3

  implicit none

  PRIVATE

  public gaugeconfiguration
  
  !> Module name
  character(len=22), parameter, public ::  modulename='gaugeconfiguration_su3'
  
  !>@brief Gauge configuration in temporal gauge
  !!@details Electric field is safe in lattice units as
  !!\f$\tilde E_{i,\vec{v}}=ga_ta_S(i)E_{i}(\vec{v})\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  type GaugeConfiguration
     !> Gauge links \f$U_\mu(x)\f$
     complex(fp), allocatable, private :: links(:,:,:,:)
     !> Electric field
     real(fp),    allocatable, private :: efield(:,:,:)
     
   contains ! Member functions
     ! Allocation and deallocation
     procedure, public :: Allocate
     procedure, public :: Deallocate

     ! Boundary value communication
     procedure, public :: CommunicateBoundary
     procedure, public :: CommunicateBoundary_Links
     procedure, public :: CommunicateBoundary_Efield
     
     ! Return and setting of member variables
     procedure, public :: GetLink
     procedure, private:: SetLink
     procedure, private:: GetEfield
     procedure, private:: SetEfield

     ! Initialisation routines
     procedure, public :: HotInit
     procedure, public :: ColdInit
     procedure, public :: TransversePolarisedOccupiedInit_Box
     procedure, private :: TransversePolarisedOccupiedInit

     ! Gauss law deviation
     procedure, public :: GetDeviationFromGausslaw

     ! Energy
     procedure, public :: GetEnergy

     ! Field strength tensor and plaquettes
     procedure, public :: GetFieldStrengthTensor
     procedure, public :: GetPlaquette
     procedure, private :: GetSpatialPlaquette
     procedure, private :: GetTemporalPlaquette

     ! Return of semi-physical fields
     procedure, public :: GetGaugefield_AlgebraCoordinate
     procedure, public :: GetGaugefield_AlgebraCoordinates
     generic, public :: GetElectricField => &
          GetElectricField_AlgebraCoordinate, GetElectricField_AlgebraMatrix
     procedure, private :: GetElectricField_AlgebraMatrix
     procedure, private :: GetElectricField_AlgebraCoordinate
     
     ! AA- and EE-correlator
     procedure, public :: GetTransverseAACorrelator
     procedure, public :: GetTransverseEECorrelator
     !procedure, public :: GetAACorrelator

     ! Update routines
     generic, public :: Update => Update_Leapfrog
     procedure, private :: Update_Leapfrog
     procedure, private :: Update_Efield_Leapfrog
     procedure, private :: Update_Links_Leapfrog
  end type GaugeConfiguration


  
contains ! Module procedures
  
  !>@brief Allocation of gaugeconfiguration
  !!@details Allocates link variables and electric field in gaugeconfiguration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 IMpure subroutine Allocate(GaugeConf)
    use lattice, only: nDim, GetMemorySize
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf

    allocate(GaugeConf%links(nSUN,nSUN,nDim,&
         GetMemorySize()))
    allocate(GaugeConf%efield(nGen,nDim,&
         GetMemorySize()))
  end subroutine Allocate
  
  !>@brief Deallocation of gaugeconfiguration
  !!@details Deallocates link variables and electric field in gaugeconfiguration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 IMpure subroutine Deallocate(GaugeConf)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    if(allocated(GaugeConf%links))  deallocate(GaugeConf%links)
    if(allocated(GaugeConf%efield)) deallocate(GaugeConf%efield)
  end subroutine Deallocate

  !>@brief Access routine to links
  !!@returns Link at given index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 IMpure function GetLink(GaugeConf,i,LatticeIndex)
    use matrixoperations, only: GetUnitmatrix
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Link variable
    complex(fp) :: GetLink(nSUN,nSUN)
    if(i/=0) then
       GetLink = GaugeConf%Links(:,:,i,GetMemoryIndex(LatticeIndex))
    else
       GetLink = GetUnitmatrix(nSUN)
    end if
  end function GetLink

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 IMpure subroutine SetLink(GaugeConf,i,LatticeIndex,Link)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Value for link variable
    complex(fp),    intent(in) :: Link(nSUN,nSUN)
    GaugeConf%Links(:,:,i,GetMemoryIndex(LatticeIndex)) = Link
  end subroutine SetLink

  !>@brief Access routine to the electric field
  !!@returns Efield at given index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 IMpure elemental function GetEfield(GaugeConf,a,i,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Efield
    real(fp) :: GetEfield
    if(i/=0) then
       GetEfield = GaugeConf%Efield(a,i,GetMemoryIndex(LatticeIndex))
    else
       GetEfield = 0._fp
    end if
  end function GetEfield

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 IMpure subroutine SetEfield(GaugeConf,a,i,LatticeIndex,Efield)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Value for link variable
    real(fp),       intent(in) :: Efield
    GaugeConf%Efield(a,i,GetMemoryIndex(LatticeIndex)) = Efield
  end subroutine SetEfield

  !>@brief Communication routine for boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary(GaugeConf)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    call GaugeConf%CommunicateBoundary_Links
    call GaugeConf%CommunicateBoundary_Efield
  end subroutine CommunicateBoundary
  
  !>@brief Communication routine for boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_Efield(GaugeConf)
    use HaloComm, only: HaloComm_CommunicateBoundary => CommunicateBoundary
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    call HaloComm_CommunicateBoundary(GaugeConf%efield)
  end subroutine CommunicateBoundary_Efield
  
  !>@brief Communication routine for boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_Links(GaugeConf)
    use HaloComm, only: HaloComm_CommunicateBoundary => CommunicateBoundary
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    call HaloComm_CommunicateBoundary(GaugeConf%links)
  end subroutine CommunicateBoundary_Links

  !>@brief Initialises the links with unit matrices and electric field with zeroes
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  IMpure subroutine ColdInit(GaugeConf)
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf

    integer(int64) :: MemoryIndex
    integer(int8)  :: i
    
    call GaugeConf%Allocate

    GaugeConf%Efield = 0
    forall(i=1:size(GaugeConf%Links,3), MemoryIndex=1:size(GaugeConf%Links,4))
       GaugeConf%Links(:,:,i,MemoryIndex) = GetUnitMatrix(nSUN)
    end forall
  end subroutine ColdInit

  impure subroutine HotInit(GaugeConf)
    use random, only: GetRandomUniformReal
    use lattice, only: ndim, GetLatticeIndex, GetProc, GetMemorySize, GetLatticeSpacing
    implicit none
    class(GaugeConfiguration), intent(out) :: GaugeConf

    integer(int64) :: MemoryIndex, LatticeIndex
    integer(int8)  :: i
    real(fp) :: r(ngen)

    
    call GaugeConf%Allocate
    
    do MemoryIndex=1,GetMemorySize()
       LatticeIndex = GetLatticeIndex(MemoryIndex)
       if(ThisProc()==GetProc(LatticeIndex)) then
          do i=1,ndim
             ! Link
             r = GetRandomUniformReal(int(ngen,int64))*GetLatticeSpacing(i)
             GaugeConf%Links(:,:,i,MemoryIndex) = GetGroupExp(r)

             ! E-field
             r = GetRandomUniformReal(int(ngen,int64))*GetLatticeSpacing(i)*GetLatticeSpacing(0_int8)
             GaugeConf%Efield(:,i,MemoryIndex) = r
          end do
       end if
    end do

    call GaugeConf%CommunicateBoundary
  end subroutine HotInit

  !>@brief Initialises the configuration as in the shape of a box highly occupied,
  !! transverse polarised fields
  !!@details Initialises the configuration with gaussian initial conditions which fills
  !! the gluon occupation up to a saturation scale \f$q_s\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  impure subroutine TransversePolarisedOccupiedInit_Box(GaugeConf,SaturationScale,Amplitude,Coupling,&
       aa_correlator,ee_correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use lattice, only: GetNorm2Momentum,GetMemorySize,GetLatticeIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> Saturation scale \f$q_s\f$
    real(fp),                  intent(in)  :: SaturationScale
    !> Occupation of the box for \f$p<q_s\f$
    real(fp),                  intent(in)  :: Amplitude
    !> Coupling
    real(fp),                  intent(in)  :: Coupling
    !> AA-Correlator
    real(fp), optional, allocatable, intent(out) :: aa_correlator(:)
    !> EE-Correlator
    real(fp), optional, allocatable, intent(out) :: ee_correlator(:)
    
    real(fp), allocatable :: Occupation(:)
    integer(int64) :: LatticeIndex, MemoryIndex
    
    allocate(Occupation(GetMemorySize()))
    !forall(is=1:size(LocalLatticeIndices))
    do MemoryIndex=1,GetMemorySize()
       Occupation(MemoryIndex)&
            = GetBoxOccupation(&
            GetNorm2Momentum(GetLatticeIndex(MemoryIndex)),& !|p|
            SaturationScale,&
            Amplitude,&
            Coupling)
    end do
    !end forall

    if(present(aa_correlator).and.present(ee_correlator)) then
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation,aa_correlator,ee_correlator)
    elseif(present(aa_correlator).and..not.present(ee_correlator)) then
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation,aa_correlator=aa_correlator)
    elseif(.not.present(aa_correlator).and.present(ee_correlator)) then
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation,ee_correlator=ee_correlator)
    else
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation)
    end if

    deallocate(Occupation)
    
  end subroutine TransversePolarisedOccupiedInit_Box

  !>@brief Initialises the configuration as highly occupied, transverse polarised fields
  !!@details Initialises the configuration according to a given gluon occupation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  impure subroutine TransversePolarisedOccupiedInit(GaugeConf,Occupation,aa_correlator,ee_correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use tolerances, only: GetZeroTol
    use mpiinterface, only: intmpi, ThisProc, mpistop, SyncAll
    use lattice, only: nDim, GetMemorySize, GetMemoryIndex, GetLatticeSize,&
         GetProc, GetVolume, GetPolarisationVectors, GetNorm2Momentum, GetMomentum,&
         GetLatticeSpacing, GetLatticeIndex, GetLocalLatticeSize
    use random, only: IsModuleInitialised_random=>IsModuleInitialised,&
         GetRandomNormalCmplx_specificProcess, modulename_random=>modulename
    use xpfft, only: p2x, x2p
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> Gauge particle occupation
    real(fp),                  intent(in)  :: Occupation(:)
    !> AA-Correlator
    real(fp), optional, allocatable, intent(out) :: aa_correlator(:)
    !> EE-Correlator
    real(fp), optional, allocatable, intent(out) :: ee_correlator(:)

    ! Gauge and electric field
    complex(fp), allocatable :: afield(:,:,:), efield(:,:,:)

    ! Prefactors
    complex(fp) :: prefactor_afield, prefactor_efield

    ! Number of transverse polarisations
    integer(int8), parameter :: nPol_transverse = ndim - 1_int8
    ! Transverse polarisation vectors
    complex(fp), dimension(nDim,nDim) :: PolarisationVectors
    ! Random numbers for fluctuations in mode decomposition
    complex(real64) :: r_afield(npol_transverse), r_efield(npol_transverse)
    ! Momentum
    complex(fp), dimension(nDim) :: Momentum(nDim)
    ! Norm of momentum
    real(fp) :: MomentumNorm

    ! Lattice
    integer(int64) :: LatticeIndex, MemoryIndex, is
    integer(int8) :: i,a
    ! MPI
    integer(intmpi) :: RecvProc
    integer(intmpi), parameter :: SendProc=0_intmpi

    ! Link
    complex(fp) :: Link(nSUN,nSUN)
    ! E-field
    real(fp) :: Efield_site(ngen)
    ! A-field
    real(fp) :: Afield_site(ngen)
    character(len=105) :: errormessage

    complex(fp) :: afield_p(ndim),efield_p(ndim)
    
    ! Check if random numbers are initialised
    if(.not.isModuleInitialised_random()) then
       errormessage = 'Error in TransversePolarisedOccupiedInit of '//modulename&
               //':'// modulename_random //' is not initialised.'
       call MPISTOP(errormessage)
    end if
    
    ! ..--** START: Allocations **--..
    allocate(afield(GetMemorySize(),nGen,nDim))
    allocate(efield(GetMemorySize(),nGen,nDim))
    ! ..--**  END : Allocations **--..

    ! ..--** START: Drawing random numbers and setting the field in p-space **--..
    !
    !        Note : All processes are doing this simultaneously together
    !               in order to ensure that the results do not depend on
    !               the number of used processes
    !
    GlobalLattice: do LatticeIndex=1,GetLatticeSize()
       RecvProc = GetProc(LatticeIndex)

       if(ThisProc()==RecvProc .or. ThisProc()==SendProc) then

          MomentumNorm = GetNorm2Momentum(LatticeIndex)
          if(MomentumNorm > GetZeroTol()) then
             if(ThisProc()==RecvProc) then
                Momentum = GetMomentum(LatticeIndex)
                call GetPolarisationVectors(Momentum,PolarisationVectors)

                MemoryIndex = GetMemoryIndex(LatticeIndex)

                prefactor_afield = &
                     sqrt(GetVolume()*Occupation(MemoryIndex)/MomentumNorm/2)

                prefactor_efield = &
                     cmplx(0,sqrt(GetVolume()*Occupation(MemoryIndex)*MomentumNorm/2),fp)
             end if ! is recieving process

             generator: do a=1_int8,ngen
                r_afield = GetRandomNormalCmplx_specificProcess(int(ngen,int64),SendProc,RecvProc)

                if(ThisProc()==RecvProc) then
                   r_efield = r_afield
                   
                   dims: do i=1,ndim !concurrent(i=1_int8:ndim)
                      afield(MemoryIndex,a,i) = &
                           prefactor_afield*(&
                                ! First transversal polarisation
                           + r_afield(1)*PolarisationVectors(i,1) &
                                ! Second transversal polarisation
                           + r_afield(2)*PolarisationVectors(i,2) )

                      efield(MemoryIndex,a,i) = &
                           prefactor_efield*(&
                                ! First transversal polarisation
                           + r_efield(1)*PolarisationVectors(i,1) &
                                ! Second transversal polarisation
                           + r_efield(2)*PolarisationVectors(i,2) )
                   end do dims
                end if ! is recieving process
             end do generator ! Generators
          
          else !p=0
             if(ThisProc()==RecvProc) then
                MemoryIndex = GetMemoryIndex(LatticeIndex)
                afield(MemoryIndex,:,:) = 0
                efield(MemoryIndex,:,:) = 0
             end if !recieving process
          end if !p>0
       end if !recieving or sending process

       call SyncAll
    end do GlobalLattice
    ! ..--**  END : Drawing random numbers and setting the field in p-space **--..

    ! Symmetrisation A(p) = A(-p)* in order to make A(x) real
    do i=1,ndim
       do a=1,ngen
          call SymmetriseInPspace(afield(:,a,i))
          call SymmetriseInPspace(efield(:,a,i))
       end do !a
    end do !i

    if(present(aa_correlator)) then
       allocate(aa_correlator(GetLocalLatticeSize()))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex(MemoryIndex)
          if(GetProc(LatticeIndex)==ThisProc()) then
             is = is+1
             aa_correlator(is) = 0
             do a=1,ngen
                afield_p = afield(MemoryIndex,a,:)
                aa_correlator(is) = aa_correlator(is) &
                     + GetTransverseField_usingProjector(afield_p,LatticeIndex)&
                     /GetVolume()/nGen
             end do
          end if
       end do
    end if
    if(present(ee_correlator)) then
       allocate(ee_correlator(GetLocalLatticeSize()))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex(MemoryIndex)
          if(GetProc(LatticeIndex)==ThisProc()) then
             is = is+1
             ee_correlator(is) = 0
             do a=1,ngen
                efield_p = efield(MemoryIndex,a,:)
                ee_correlator(is) = ee_correlator(is) &
                     + GetTransverseField_usingProjector(efield_p,LatticeIndex)&
                     /GetVolume()/nGen
             end do
          end if
       end do
    end if
    
    !..--** START: FFT p-->x **--..
    do i=1,ndim
       do a=1,ngen
          call p2x(afield(:,a,i))
          call p2x(efield(:,a,i))
       end do !a
    end do !i
    !..--**  END : FFT p-->x **--..

    
    !..--** START: Writing fields to configuration **--..
    call GaugeConf%Allocate
    !do concurrent(&
    !     is=1:size(LocalLatticeIndices),&
    !     i =1:ndim)
    do MemoryIndex=1,GetMemorySize()
       do i=1,ndim
          LatticeIndex = GetLatticeIndex(MemoryIndex)
          if(GetProc(LatticeIndex)==ThisProc()) then
             ! Link
             afield_site = real(afield(MemoryIndex,:,i),fp) &
                                ! Translation to lattice units
                  *GetLatticeSpacing(i)
             GaugeConf%Links(:,:,i,MemoryIndex) = GetGroupExp(afield_site)
             
             ! E-field
             efield_site = real(efield(MemoryIndex,:,i),fp) &
                                ! Translation to lattice units
                  *GetLatticeSpacing(0_int8)*GetLatticeSpacing(i)
             GaugeConf%Efield(:,i,MemoryIndex) = efield_site
          end if
       end do
    end do
    !..--**  END : Writing fields to configuration **--..
    call GaugeConf%CommunicateBoundary()
  contains
    impure subroutine SymmetriseInPspace(field)
      use precision, only: fp
      use lattice, only: nDim, GetProc, GetNegativeLatticeIndex, GetLatticeSize
      use mpiinterface, only: ThisProc, SyncAll, intmpi, GetComplexSendType
      use mpi
      implicit none
      complex(fp), intent(inout) :: field(:)


      ! MPI
      complex(fp) :: a
      integer(intmpi) :: status(mpi_status_size)
      integer(intmpi), parameter :: buffersize=1_intmpi
      integer(intmpi) :: dest, source
      integer(intmpi) :: tag
      integer(intmpi) :: mpierr

      ! Indices
      integer(intmpi) :: PositiveProc, NegativeProc
      integer(int64)  :: PositiveLatticeIndex, NegativeLatticeIndex
      integer(int64)  :: PositiveMemoryIndex, NegativeMemoryIndex
      integer(int64)  :: LatticeSize

      LatticeSize = GetLatticeSize()

      do PositiveLatticeIndex=1,LatticeSize
         ! Getting negative lattice index, corresponding to -p
         NegativeLatticeIndex = GetNegativeLatticeIndex(PositiveLatticeIndex)

         ! Getting process numbers for positive and negative momentum
         PositiveProc = GetProc(PositiveLatticeIndex)
         NegativeProc = GetProc(NegativeLatticeIndex)

         ! Setting field such that f(p) = f(-p)*

         if(PositiveProc == NegativeProc) then
            ! positive and negative process are the same --> No MPI-communication needed
            if(ThisProc()==PositiveProc) then
               PositiveMemoryIndex = GetMemoryIndex(PositiveLatticeIndex)
               NegativeMemoryIndex = GetMemoryIndex(NegativeLatticeIndex)
               
               a = field(PositiveMemoryIndex)
               if(PositiveMemoryIndex==NegativeMemoryIndex) then
                  field(PositiveMemoryIndex) = real(a,fp)
               else
                  field(PositiveMemoryIndex) = a
                  field(NegativeMemoryIndex) = conjg(a)
               end if
            end if
         else 
            ! positive and negative process are not the same --> MPI-communication
            tag  = PositiveLatticeIndex
            source = PositiveProc
            dest = NegativeProc
            
            if(ThisProc()==PositiveProc) then
               PositiveMemoryIndex = GetMemoryIndex(PositiveLatticeIndex)
               a = field(PositiveMemoryIndex)

               call MPI_Send(&
                    a,                   & ! What to send
                    buffersize,          & ! How many points to send
                    GetComplexSendType(),& ! What type to send
                    dest,                & ! Destination
                    tag,                 & ! Tag
                    MPI_COMM_WORLD,      & ! Communicator
                    mpierr)                ! Error code
            elseif(ThisProc()==NegativeProc) then
               call MPI_Recv(&
                    a,                   & ! What to send
                    buffersize,          & ! How many points to send
                    GetComplexSendType(),& ! What type to send
                    source,              & ! Destination
                    tag,                 & ! Tag
                    MPI_COMM_WORLD,      & ! Communicator
                    status,              & ! Status
                    mpierr)                ! Error code
               
               NegativeMemoryIndex = GetMemoryIndex(NegativeLatticeIndex)
               field(NegativeMemoryIndex) = conjg(a)
            else
               ! do nothing
            end if
         end if

         call SyncAll
      end do
    end subroutine SymmetriseInPspace
  end subroutine TransversePolarisedOccupiedInit

  !>@brief Box occupation
  !!@details Step function with amplitude until saturation scale. Afterwards coupling\f$^2/2\f$
  !!@returns Box occupation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
 IMpure elemental real(fp) function GetBoxOccupation(Momentum,SaturationScale,Amplitude,Coupling)
    use precision, only: fp
    implicit none
    !> Lattice momentum norm
    real(fp), intent(in) :: Momentum
    !> Saturation scale \f$q_s\f$
    real(fp), intent(in)  :: SaturationScale
    !> Occupation of the box for \f$p<q_s\f$
    real(fp), intent(in)  :: Amplitude
    !> Coupling
    real(fp), intent(in)  :: Coupling
    if( Momentum < SaturationScale ) then
       GetBoxOccupation = Amplitude
    else
       GetBoxOccupation = Coupling**2/2
    end if
  end function GetBoxOccupation
    

  !>@brief Total deviation from Gauss-law of the total MPI-distributed configuration
  !!@returns Total deviation from Gauss-law of the total MPI-distributed configuration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure real(fp) function GetDeviationFromGaussLaw(GaugeConf)
    use mpiinterface, only: intmpi, GetRealSendType
    use mpi
    use lattice, only: nDim, GetLatticeSpacing,GetNeib,GetLatticeIndex, GetNeib, GetMemorySize
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    real(fp), dimension(ngen) :: efield, efield_neib
    complex(fp), dimension(nsun,nsun) :: Mefield, Mefield_neib, derivative, Link_neib
    
    integer(int8)  :: i, a
    integer(int64) :: MemoryIndex, neib, latticeindex
    integer(int64), allocatable :: LocalLatticeIndices(:)

    ! 1. Calculation of local contribution
    local_contribution = 0
    !do concurrent (is=1:size(LocalLatticeIndices))
    do MemoryIndex=1,GetMemorySize()
       LatticeIndex = GetLatticeIndex(MemoryIndex)
       !do concurrent(i=1_int8:ndim)
       do i=1,ndim
          efield = GaugeConf%GetEfield([1_int8:ngen],i,LatticeIndex)
          
          neib = GetNeib(-i,LatticeIndex)
          efield_neib = GaugeConf%GetEfield([1_int8:ngen],i,Neib)

          Mefield = GetAlgebraMatrix(efield)
          Mefield_neib = GetAlgebraMatrix(efield_neib)

          Link_neib = GaugeConf%GetLink(i,Neib)

          derivative = (Mefield &
               - matmul(matmul(&
               conjg(transpose(Link_Neib)),&
               Mefield_neib),&
               Link_neib))&
               /GetLatticeSpacing(i)**2/GetLatticeSpacing(0)
          do concurrent(a=1_int8:ngen)
             local_contribution = local_contribution &
                  + 2*Abs(Aimag(GetTraceWithGenerator(a,derivative)))
          end do
       end do
    end do

    ! 2. MPI-Sum over all partitions
    call MPI_ALLREDUCE(&
         local_contribution,&
         GetDeviationFromGaussLaw,&
         1_intmpi,&
         GetRealSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)
  end function GetDeviationFromGaussLaw

  !>@brief Returns field strength tensor
  !!@details: Field strength tensor in physical units \f$g F_{ij}(\vec{v})^aT^a\f$
  !!in order \f$O(a_i^2,a_j^2,a_i\cdot a_j)\f$
  !!@returns Spatial field strength tensor \f$F_{i,j,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 23.02.2019
  !!@version 1.0
 IMpure function GetFieldStrengthTensor(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetLatticeSpacing
    use matrixoperations, only: LogU
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Field strength tensor
    complex(fp) :: GetFieldStrengthTensor(nSUN,nSUN)

    complex(fp) :: Plaquette(nSUN,nSUN)

    Plaquette = GaugeConf%GetPlaquette(i,j,LatticeIndex)
    
    GetFieldStrengthTensor = &
         LogU(Plaquette)/GetLatticeSpacing(i)/GetLatticeSpacing(j)/cmplx(0,1,fp)
  end function GetFieldStrengthTensor

  !> @brief Returns spatial plaquette \f$U_{ij,\vec{v}} =U_{i,\vec{v}}\cdot U_{j,\vec{v}+\hat{i}}
  !! \cdot U_{i,\vec{v}+\hat{j}}^\dagger\cdot U_{j,\vec{v}}^\dagger\f$
  !! @returns Spatial plaquette
  !! \f$U_{i,j,v}=U_{i,v}\cdot U_{j,v+\hat{i}}\cdot U_{i,v+\hat{j}}^\dagger\cdot U_{j,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 23.02.2019
  !! @version 1.0
 IMpure function GetSpatialPlaquette(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetNeib
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Spatial direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Spatial plaquette
    complex(fp) :: GetSpatialPlaquette(nSUN,nSUN)

    complex(fp), dimension(nSUN,nSUN) :: link1,link2,link3,link4

    link1 = GaugeConf%GetLink(i,LatticeIndex)                             ! U_i(x)
    link2 = GaugeConf%GetLink(j,GetNeib(i,LatticeIndex))                  ! U_j(x+î)
    link3 = conjg(transpose(GaugeConf%GetLink(i,GetNeib(j,LatticeIndex))))! U_i(x+ĵ)†
    link4 = conjg(transpose(GaugeConf%GetLink(j,LatticeIndex)))           ! U_j(x)†

    GetSpatialPlaquette = matmul(matmul(matmul(&
         link1,&  ! U_i(x)
         link2),& ! U_j(x+î)
         link3),& ! U_i(x+ĵ)†
         link4)   ! U_j(x)†
  end function GetSpatialPlaquette

  !> @brief Returns temporal plaquette
  !! @details The temporal plaquette is computed via the chromo-electric field
  !! due to the knowledge of the evolution scheme. Using that
  !! \f$U_{i,\vec{v}}(t+a_t) =\exp\left(i \tilde{E}_{C,i,\vec{v}}(t)^aT^a\right)
  !! \cdot U_{i,\vec{v}}(t)\f$
  !! one can exactly compute the temporal plaquette with
  !! \f$U_{0j,\vec{v}}(t)=U_{j,\vec{v}}(t+a_t)
  !! \cdot U_{j,\vec{v}}(t)^\dagger
  !! =\exp\left(i \tilde{E}_{C,i,\vec{v}}(t)^aT^a\right)\f$
  !! @returns Temporal plaquette
  !! \f$U_{0,j,v}=U_{0,v}\cdot U_{j,v+\hat{0}}
  !! \cdot U_{0,v+\hat{j}}^\dagger\cdot U_{j,v}^\dagger
  !! = U_{j,v+\hat{0}}\cdot U_{j,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 23.02.2019
  !! @version 1.0
 IMpure function GetTemporalPlaquette(GaugeConf,i,j,LatticeIndex)
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Temporal plaquette
    complex(fp) :: GetTemporalPlaquette(nSUN,nSUN)

    real(fp) :: efield(ngen)
    
    if(     i==0 .and. j/=0 ) then
       efield = +GaugeConf%GetEfield([1_int8:ngen],j,LatticeIndex)
       GetTemporalPlaquette = GetGroupExp(efield)
    elseif( i/=0 .and. j==0 ) then
       efield = -GaugeConf%GetEfield([1_int8:ngen],i,LatticeIndex)
       GetTemporalPlaquette = GetGroupExp(efield)
    else
       GetTemporalPlaquette = GetUnitMatrix(nSUN)
    end if
  end function GetTemporalPlaquette
  
  !> @brief Returns plaquette \f$U_{jk,\vec{v}}
  !! =U_{j,\vec{v}}\cdot U_{k,\vec{v}+\hat{j}}
  !! \cdot U_{j,\vec{v}+\hat{k}}^\dagger\cdot U_{k,\vec{v}}^\dagger\f$
  !! @returns Plaquette
  !! \f$U_{j,k,v}=U_{j,v}\cdot U_{k,v+\hat{j}}
  !! \cdot U_{j,v+\hat{k}}^\dagger\cdot U_{k,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 23.02.2019
  !! @version 1.0
 IMpure function GetPlaquette(GaugeConf,i,j,LatticeIndex)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Temporal plaquette
    complex(fp) :: GetPlaquette(nSUN,nSUN)
    if(i==0 .or. j==0) then
       GetPlaquette = GaugeConf%GetTemporalPlaquette(i,j,LatticeIndex)
    else
       GetPlaquette = GaugeConf%GetSpatialPlaquette(i,j,LatticeIndex)
    end if
  end function GetPlaquette

  !>@brief Energy functional
  !!@returns Energy functional
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure real(fp) function GetEnergy(GaugeConf)
    use precision, only: fp
    use matrixoperations, only: GetTrace
    use lattice, only: nDim, GetLatticeSpacing, GetMemorySize, GetLatticeIndex
    use mpiinterface, only: intmpi, GetRealSendType
    use mpi
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    real(fp), dimension(ngen) :: Efield
    complex(fp), dimension(nsun,nsun) :: Plaquette

    real(fp) :: PotentialTerm
    integer(int8)  :: i, j
    integer(int64) :: MemoryIndex, latticeindex

    ! 1. Calculation of local contribution
    local_contribution = 0
    !do concurrent (is=1:size(LocalLatticeIndices))
    do MemoryIndex=1,GetMemorySize()
       LatticeIndex = GetLatticeIndex(MemoryIndex)
       !do concurrent(i=1_int8:ndim)
       do i=1,ndim
          efield = GaugeConf%GetEfield([1_int8:ngen],i,LatticeIndex)
          
          local_contribution = local_contribution &
               ! Electric energy
               + sum(efield**2)/(GetLatticeSpacing(i)*GetLatticeSpacing(0_int8))**2/2
          ! Magnetic energy
          !do concurrent(j=i+1_int8:ndim)
          do j=i+1_int8,ndim
             Plaquette = GaugeConf%GetPlaquette(i,j,LatticeIndex)
             PotentialTerm = (1-real(GetTrace(Plaquette),fp)/nSUN)&
                  /GetLatticeSpacing(i)/GetLatticeSpacing(j)

             Local_contribution = Local_Contribution &
                  + 2*nSUN*PotentialTerm
          end do
       end do
    end do

    ! 2. MPI-Sum over all partitions
    call MPI_ALLREDUCE(&
         local_contribution,&
         GetEnergy,&
         1_intmpi,&
         GetRealSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)
  end function GetEnergy

  !>@brief Computes the transverse AA-correlator (MPI-distributed)
  !!@returns Transverse AA-correlator (MPI-distributed)
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
  impure subroutine GetTransverseAACorrelator(GaugeConf,Correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use lattice, only: nDim, GetLatticeIndex, GetMemoryIndex, GetProc, GetMemorySize,&
         GetLocalLatticeSize
    use mpiinterface, only: ThisProc
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in)  :: GaugeConf
    !> Transverse AA-correlator (MPI-distributed)
    real(fp), allocatable, intent(out) :: Correlator(:)

    integer(int8)  :: a, i
    integer(int64) :: MemoryIndex, LatticeIndex

    real(fp), allocatable :: field(:,:,:), correlator_a(:)
    
    ! Extracting the gauge field in position space
    allocate(field(GetMemorySize(),nDim,nGen))
    !forall(&
    !     LocalIndex=1_int64:Size(LocalLatticeIndices_includingHalo),&
    !     i=1_int8:nDim )
    do MemoryIndex=1,GetMemorySize()
       do i=1,ndim
          field(MemoryIndex,i,:)&
               = GaugeConf%GetGaugeField_AlgebraCoordinates(&
               i,GetLatticeIndex(MemoryIndex))
       end do
    end do
    !end forall
    
    ! Computing average over all generators
    allocate(correlator(GetLocalLatticeSize()))
    correlator = 0
    
    do a=1,nGen
       call GetTransverseCorrelator(field(:,:,a),correlator_a)
       correlator = correlator + correlator_a/nGen
    end do
    deallocate(correlator_a,field)
  end subroutine GetTransverseAACorrelator

  !>@brief Computes the transverse EE-correlator (MPI-distributed)
  !!@returns Transverse EE-correlator (MPI-distributed)
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
  impure subroutine GetTransverseEECorrelator(GaugeConf,Correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use lattice, only: nDim, GetLatticeIndex, GetMemoryIndex, GetProc, GetMemorySize,&
         GetLocalLatticeSize
    use mpiinterface, only: ThisProc
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in)  :: GaugeConf
    !> Transverse AA-correlator (MPI-distributed)
    real(fp), allocatable, intent(out) :: Correlator(:)

    integer(int8)  :: a, i
    integer(int64) :: MemoryIndex, LatticeIndex

    real(fp), allocatable :: field(:,:,:), correlator_a(:)
    
    ! Extracting the gauge field in position space
    allocate(field(GetMemorySize(),nDim,nGen))
    !forall(&
    !     LocalIndex=1_int64:Size(LocalLatticeIndices_includingHalo),&
    !     i=1_int8:nDim )
    do MemoryIndex=1,GetMemorySize()
       do i=1,ndim
          field(MemoryIndex,i,:)&
               = GaugeConf%GetElectricField_AlgebraCoordinate(&
               [1_int8:ngen],i,GetLatticeIndex(MemoryIndex))
       end do
    end do
    !end forall

    ! Computing average over all generators
    allocate(correlator(GetLocalLatticeSize()))
    correlator = 0
    
    do a=1,nGen
       call GetTransverseCorrelator(field(:,:,a),correlator_a)
       correlator = correlator + correlator_a/nGen
    end do
    deallocate(correlator_a,field)
  end subroutine GetTransverseEECorrelator

  !>@brief Computes the transverse correlator (MPI-distributed)
  !!@returns Transverse correlator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.01.2019
  !!@version 1.0
  impure subroutine GetTransverseCorrelator(x_field,correlator)
    use lattice, only: nDim, GetVolume, GetLatticeIndex, GetMemorySize, GetLocalLatticeSize, GetProc
    use xpfft,     only: x2p
    use precision, only: fp
    use mpiinterface, only: ThisProc
    implicit none
    !> Field in position space (e.g. A(x), E(x)), indexed as: field(latticeindex,direction)
    real(fp),              intent(in)  :: x_field(:,:)
    !> Transverse correlator (MPI-distributed)
    real(fp), allocatable, intent(out) :: correlator(:)


    integer(int8) :: i
    integer(int64) :: LatticeIndex, MemoryIndex, is

    complex(fp), allocatable :: p_field(:,:)
    complex(fp) :: field_site(nDim)

    ! Allocating fourier transform of input
    allocate(p_field(GetMemorySize(),nDim))
    p_field = cmplx(x_field,0,fp)

    ! Performing fourier transform for all spatial components
    do i=1,ndim
       call x2p(p_field(:,i))
    end do
    
    ! Allocating output
    allocate(correlator(GetLocalLatticeSize()))
    correlator = 0
    !do concurrent(LocalIndex=1:size(correlator))
    is=0
    do MemoryIndex=1,GetMemorySize()
       LatticeIndex = GetLatticeIndex(MemoryIndex)
       if(GetProc(LatticeIndex)==ThisProc()) then
          is=is+1
          field_site = p_field(MemoryIndex,:)
          Correlator(is) = &
               GetTransverseField_usingProjector&
               (field_site,LatticeIndex)/GetVolume()
       end if
    end do
    deallocate(p_field)
  end subroutine GetTransverseCorrelator
  IMpure real(fp) function GetTransverseField_usingProjector(field,LatticeIndex)
    use precision, only: fp
    use lattice, only: nDim, GetMomentum, GetTransverseProjector
    implicit none
    !> Field in p-space
    complex(fp),    intent(in) :: field(nDim)
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex

    integer(int8) :: i,j
    complex(fp) :: res, momentum(nDim)

    res = 0
    !do concurrent(i=1:ndim, j=1:ndim)
    do i=1,ndim
       do j=1,ndim
          momentum = GetMomentum(LatticeIndex)
          res = res + field(i)*GetTransverseProjector(i,j,momentum)*conjg(field(j))
       end do
    end do
    GetTransverseField_usingProjector = real(res,fp)
  end function GetTransverseField_usingProjector

  !>@brief Returns algebra-component of the gauge field
  !!@returns \f$A_{i,\vec{v}}^a\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 IMpure elemental function GetGaugefield_AlgebraCoordinate(GaugeConf,a,i,LatticeIndex)
    use precision, only: fp
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: LatticeIndex
    !> Algebra coordinates of electric field
    real(fp) :: GetGaugefield_AlgebraCoordinate

    real(fp) :: algebra_coordinates(ngen)

    algebra_coordinates = GaugeConf%GetGaugefield_AlgebraCoordinates(i,LatticeIndex)
    GetGaugefield_AlgebraCoordinate = algebra_coordinates(a)/GetLatticeSpacing(i)
  end function GetGaugefield_AlgebraCoordinate

  !>@brief Returns all algebra-components of the gauge field
  !!@returns  All \f$A_{i,\vec{v}}^a\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 IMpure function GetGaugefield_AlgebraCoordinates(GaugeConf,i,latticeindex)
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: latticeindex
    !> Algebra coordinates of electric field
    real(fp) :: GetGaugefield_AlgebraCoordinates(ngen)

    complex(fp) :: Link(nSUN,nSUN)

    Link = GaugeConf%GetLink(i,LatticeIndex)

    GetGaugefield_AlgebraCoordinates = GetGroupLog(Link)/GetLatticeSpacing(i)
  end function GetGaugefield_AlgebraCoordinates

  !>@brief Returns algebra-component of electric field
  !!@returns \f$E_{i,\vec{v}}^a\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 IMpure elemental function GetElectricField_AlgebraCoordinate(conf,a,i,latticeindex)
    use precision, only: fp
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(in) :: conf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: latticeindex
    !> Algebra coordinates of electric field
    real(fp) :: GetElectricField_AlgebraCoordinate
    
    GetElectricField_AlgebraCoordinate &
         = conf%GetEfield(a,i,latticeindex)/GetLatticeSpacing(i)/GetLatticeSpacing(0_int8)
  end function GetElectricField_AlgebraCoordinate

  !>@brief Returns electric field
  !!@returns \f$E_{i,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 IMpure function GetElectricField_AlgebraMatrix(conf,i,latticeindex)
    use precision, only: fp
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(in) :: conf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: latticeindex
    !> Electric field
    complex(fp) :: GetElectricField_AlgebraMatrix(nSUN,nSUN)
    
    real(fp) :: efield(ngen)
    efield = conf%GetEfield([1_int8:ngen],i,LatticeIndex)
    
    GetElectricField_AlgebraMatrix &
         = GetAlgebraMatrix(efield)/GetLatticeSpacing(i)/GetLatticeSpacing(0_int8)
  end function GetElectricField_AlgebraMatrix

  !>@brief Update routine using the Leapfrog algorithm
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 27.02.2019
  !!@version 1.0
  impure subroutine Update_Leapfrog(GaugeConf,StepWidth_)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Step width in units of \f$a_0\f$
    real(fp),                  intent(in), optional :: Stepwidth_

    real(fp) :: Stepwidth
    
    ! ..--** START: Optional parameters **--..
    if(present(Stepwidth_)) then
       Stepwidth = Stepwidth_
    else
       ! Default value
       Stepwidth = 1._fp ! in units of a0
    end if
    ! ..--**  END : Optional parameters **--..
    
    if(Stepwidth.gt.0) then
       call GaugeConf%Update_Links_Leapfrog(Stepwidth)
       !call GaugeConf%CommunicateBoundary_Links
       call GaugeConf%Update_Efield_Leapfrog(Stepwidth)
       call GaugeConf%CommunicateBoundary_Efield
    else
       call GaugeConf%Update_Efield_Leapfrog(Stepwidth)
       call GaugeConf%Update_Links_Leapfrog(Stepwidth)
       !call GaugeConf%CommunicateBoundary
       call GaugeConf%CommunicateBoundary_Efield
    end if
  end subroutine Update_Leapfrog

  !>@brief Update routine for the links as part of the Leapfrog algorithm
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 27.02.2019
  !!@version 1.0
  impure subroutine Update_Links_Leapfrog(GaugeConf,StepWidth)
    use lattice, only: ndim, GetMemorySize
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Step width in units of \f$a_0\f$
    real(fp),                  intent(in)    :: StepWidth
    
    integer(int8) :: i
    integer(int64):: MemoryIndex

    complex(fp) :: TimeEvolutionOperator(nSUN,nSUN)
    real(fp)    :: efield_times_dt(nGen)

    do MemoryIndex=1,GetMemorySize()
       do i=1,ndim
          efield_times_dt       = GaugeConf%Efield(:,i,MemoryIndex)*StepWidth
          TimeEvolutionOperator = GetGroupExp(efield_times_dt)
          
          GaugeConf%Links(:,:,i,MemoryIndex) =&
               matmul(TimeEvolutionOperator,GaugeConf%Links(:,:,i,MemoryIndex))
       end do
    end do
  end subroutine Update_Links_Leapfrog

  !>@brief Update routine for the electric field as part of the Leapfrog algorithm
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 27.02.2019
  !!@version 1.0
  impure subroutine Update_Efield_Leapfrog(GaugeConf,StepWidth)
    use lattice, only: ndim, GetMemorySize,GetLatticeIndex
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Step width in units of \f$a_0\f$
    real(fp),                  intent(in)    :: StepWidth
    
    integer(int8) :: i
    integer(int64):: MemoryIndex,LatticeIndex
    
    do MemoryIndex=1,GetMemorySize()
       LatticeIndex=GetLatticeIndex(MemoryIndex)
       do i=1,ndim
          call Update_Efield_Leapfrog_atSite_Direction(GaugeConf,StepWidth,i,latticeindex)
       end do
    end do
    
  contains
    pure subroutine Update_Efield_Leapfrog_atSite_Direction(GaugeConf,StepWidth,i,latticeindex)
      use lattice, only: ndim, GetLatticeSpacing, GetMemoryIndex
      implicit none
      !> Gauge link configuration
      type(GaugeConfiguration), intent(inout) :: GaugeConf
      !> Step width in units of \f$a_0\f$
      real(fp),                 intent(in)    :: StepWidth
      !> Direction
      integer(int8),            intent(in)    :: i
      !> Lattice index
      integer(int64),           intent(in)    :: latticeindex

      integer(int8) :: k,a
      integer(int64) :: MemoryIndex
      complex(fp) :: staplesum(nsun,nsun), link_times_staplesum(nsun,nsun)
      
      staplesum = 0
      do k=1,ndim
         if(k /= i) then
            staplesum = staplesum + &
                 StepWidth*(GetLatticeSpacing(0)/GetLatticeSpacing(k))**2 &
                 *(&
                 GetUStaple(GaugeConf,i,k,latticeindex) + &
                 GetDStaple(GaugeConf,i,k,latticeindex) &
                 )
         end if ! i /= k
      end do !k

      MemoryIndex = GetMemoryIndex(LatticeIndex)
      link_times_staplesum = matmul(GaugeConf%Links(:,:,i,MemoryIndex),staplesum)

      forall(a=1:ngen) GaugeConf%efield(a,i,MemoryIndex)&
           = GaugeConf%efield(a,i,MemoryIndex)&
           - 2*Aimag(GetTraceWithGenerator(a,link_times_staplesum))
    end subroutine Update_Efield_Leapfrog_AtSite_Direction

    !>@brief Returns staple
    !!@details The U-staple is defined as
    !! \f$ S^{\text{U}}_{ik,x}=U_{k,x+\hat{i}}\cdot U_{i,x+\hat{k}}^\dagger\cdot U_{k,x}^\dagger\f$
    !!@autho Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 27.02.2019
    !!@version 1.0
    pure function GetUStaple(GaugeConf,i,k,latticeindex)
      use lattice, only: GetNeib, GetMemoryIndex
      implicit none
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      !> Direction
      integer(int8),            intent(in) :: i,k
      !> Lattice index
      integer(int64),           intent(in) :: latticeindex
      !> Upwards staple
      complex(fp)                          :: GetUStaple(nsun,nsun)

      integer(int64) :: memoryindex,neib_i, neib_k

      memoryindex = GetMemoryIndex(LatticeIndex)
      neib_i      = GetMemoryIndex(GetNeib(+i,latticeindex))
      neib_k      = GetMemoryIndex(GetNeib(+k,latticeindex))

      GetUStaple = matmul(matmul(&
           GaugeConf%links(                :,:,k,neib_i),&
           conjg(transpose(GaugeConf%links(:,:,i,neib_k)))),&
           conjg(transpose(GaugeConf%links(:,:,k,MemoryIndex))))
    end function GetUStaple
    
    !>@brief Returns staple
    !!@details The U-staple is defined as \f$ S^{\text{D}}_{ik,x}=
    !! U_{k,x+\hat{i}-\hat{k}}^\dagger\cdot U_{i,x-\hat{k}}^\dagger\cdot U_{k,x-\hat{k}}\f$
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 27.02.2019
    !!@version 1.0
    pure function GetDStaple(GaugeConf,i,k,latticeindex)
      use lattice, only: GetNeib, GetMemoryIndex
      implicit none
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      !> Direction
      integer(int8),            intent(in) :: i,k
      !> Lattice index
      integer(int64),           intent(in) :: latticeindex
      !> Downwards staple
      complex(fp)                          :: GetDStaple(nsun,nsun)

      integer(int64) :: neib_k, neib_ik

      neib_k      = GetMemoryIndex(GetNeib(-k,latticeindex))
      neib_ik     = GetMemoryIndex(GetNeib(+i,neib_k))

      GetDStaple = matmul(matmul(&
           conjg(transpose(GaugeConf%links(:,:,k,neib_ik))),&
           conjg(transpose(GaugeConf%links(:,:,i,neib_k)))),&
           GaugeConf%links(:,:,k,neib_k))
    end function GetDStaple
  end subroutine Update_Efield_Leapfrog
end module gaugeconfiguration_su3
