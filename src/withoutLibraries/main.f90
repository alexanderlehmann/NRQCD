!----------------------------------------------------------------------
! PROGRAMS for Lattice-NRQCD
!----------------------------------------------------------------------
!
! MODULE: programs
!>@brief Program/mains
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 07.03.2019
!!@version 1.0
! REVISION HISTORY:
! 07 03 2019 - Initial Version
!----------------------------------------------------------------------
module programs
  PUBLIC

contains
  !>@brief Program for measuring the time-evolution of the energy and gauss law deviation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 23.04.2020
  !!@version 1.1
  impure subroutine MeasureEnergyAndGaussLawDeviation
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp)   :: CoMTime
    real(fp)   :: TimeRange

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf

    ! Monitoring variables
    integer(int64) :: it
    real(fp) :: time
    
    ! Observable output (Gauss, Energy)
    integer(int8) :: fileID_eg
    real(fp) :: gauss, energy

    call InitSimulation

    call GaugeConf%HotInit()
    
    if(ThisProc()==0) &
         fileID_eg = OpenFile(filename="energy_gauss.txt",st='REPLACE',fm='FORMATTED',act='WRITE')
    do it=0,TimeSteps
       time = it*GetLatticeSpacing(0_int8)
       gauss = GaugeConf%GetDeviationFromGaussLaw()
       energy= GaugeConf%GetEnergy()
       if(ThisProc()==0) write(fileID_eg,'(3(SP,E13.6,1X))') time,energy,gauss
       if(ThisProc()==0) write(output_unit,*) time
       call GaugeConf%Update
    end do
    if(ThisProc()==0) call CloseFile(fileID_eg)

    call EndSimulation
  contains

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! Center of mass time T
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') CoMTime
      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface, only: ThisProc, FinalizeModule_MPIinterface => FinalizeModule
      implicit none

      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureEnergyAndGaussLawDeviation
  
  impure subroutine MeasureWilsonLines_Equilibrium_StaticMeson
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use wilsonline
    use random
    use nrqcd
    use statistics
    use su3, only: nsun


    use matrixoperations

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: Beta
    real(fp) :: TimeRange

    logical :: MeasureEnergy
    character(len=80) :: EnergyFilename
    integer(int64) :: nefieldinit,nequilibriumtimesteps
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_initial, GaugeConf_t
    
    ! Counting
    integer :: i

    integer(int64) :: it

    integer(intmpi) :: mpierr
    
    ! Wilson loop parameters
    integer(int8)  :: messdir
    integer(int64) :: r
    
    integer(int64) :: meanMaxLatticeIndex=20
    integer(int8) :: k
    ! Indices
    integer(int64) :: LatticeIndex_distance, shiftstep
    integer(int8) :: colour

    real(fp) :: time
    integer(int64) :: TimePoints

    integer(intmpi) :: proc
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork

    ! Charge density of static meson
    real(fp), allocatable :: ChargeDensity(:,:)
    real(fp) :: tolerance_Eprojection
    
    real(fp) :: localcharge,totalcharge

    ! Observables
    integer(int64) :: index_quark
    complex(fp), allocatable :: WilsonLoops(:)
    
    character(len=80) :: &
         FileName_WilsonLoops
    
    call InitSimulation
    
    TimePoints = nint(TimeRange/LatticeSpacings(0))

    ! Initialising the gauge configuration for the given charge density
    call GaugeConf_initial%EquilibriumInit(&
         Beta,nefieldinit,nequilibriumtimesteps,tolerance_Eprojection,&
         ChargeDensity,MeasureEnergy,EnergyFilename)

    if(thisproc()==0) write(output_unit,*) 'Done: Equilibration'
    
    ! Performing update steps while computing energy tensor
    
    allocate(WilsonLoops(0:TimePoints))

    nwork = TimePoints+1

    iwork = 0
    GaugeConf_t = GaugeConf_initial
    ! Setting time of gauge conf to t1
    do it=0,TimeSteps,+1
       WilsonLoops(it) = GetWilsonLoop(GaugeConf_initial,GaugeConf_t,r,messdir,index_quark)
       
       call GaugeConf_t%Update
       iwork = iwork + 1
       if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
    end do

    call PrintObservables
    
    call EndSimulation

  contains

    !>@brief Output of observables
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 29.05.2019
    !!@version 1.0
    impure subroutine PrintObservables
      implicit none
      
      real(fp) :: t

      character(len=3) :: adv

      integer(int8) :: rmax, FileID_WilsonLoops

      rmax = r
      
      ! Opening files and printing header
      if(ThisProc()==0) then
         ! Opening files
         FileID_WilsonLoops = OpenFile(filename=FileName_WilsonLoops,&
              st='REPLACE',fm='FORMATTED',act='WRITE')

         write(FileID_WilsonLoops,'(A1,1X)',advance='no') 't'
         
         if(r<rmax) then
            adv='no'
         else
            adv='yes'
         end if

         write(FileID_WilsonLoops,'(A5,I0.2,A1,1X)',advance='no') 'Re(r=',r,')'
         write(FileID_WilsonLoops,'(A5,I0.2,A1,1X)',advance=trim(adv)) 'Im(r=',r,')'

         do it=lbound(WilsonLoops,1),ubound(WilsonLoops,1)
            time = it*LatticeSpacings(0)
            
            write(FileID_WilsonLoops,'(SP,E16.9,1X)',advance='no') time

            if(r<rmax) then
               adv='no'
            else
               adv='yes'
            end if

            write(FileID_WilsonLoops,'(2(SP,E16.9,1X))',advance=trim(adv)) &
                 real(WilsonLoops(it)),aimag(WilsonLoops(it))
         end do

         call CloseFile(FileID_WilsonLoops)
      end if
      
    end subroutine PrintObservables

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i, j
      
      real(fp) :: c_re, c_im
      character(len=100) :: FileName_ChargeDensity
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') messdir
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I2)') r
      
      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed
      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') Beta

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') tolerance_Eprojection
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nefieldinit
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nequilibriumtimesteps
      
      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') j
      if(j==0) then
         MeasureEnergy=.false.
         arg_count = arg_count +1;
      else
         MeasureEnergy=.true.
         arg_count = arg_count +1; call get_command_argument(arg_count,EnergyFilename);
      end if

      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_ChargeDensity);

      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);
      
      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      ! Read charge density
      call ReadChargeDensity(ChargeDensity,FileName_ChargeDensity)

      call SyncAll
    end subroutine InitSimulation
    
    impure subroutine ReadChargeDensity(ChargeDensity,FileName_ChargeDensity)
      use mpiinterface
      use io
      implicit none
      real(fp), allocatable, intent(out) :: ChargeDensity(:,:)
      character(len=*), intent(in) :: FileName_ChargeDensity

      integer(int64) :: position(nDim), latticeindex
      integer(int8) :: colour
      real(fp), dimension(nSUN) :: charges

      integer(int8) :: FileID_ChargeDensity
      integer(int64) :: nLines, iLine

      index_quark = -1

      
      allocate(ChargeDensity(nsun,GetMemorySize()))
      ChargeDensity = 0

      FileID_ChargeDensity = OpenFile(FileName_ChargeDensity,st='OLD',fm='FORMATTED',act='READ')
      ! Determine number of entries in chargedensity-file
      nLines = 0
      do
         READ(FileID_ChargeDensity,*,END=10)
         nLines = nLines + 1
      end do
10    REWIND(FileID_ChargeDensity)
      
      ! Read in charge density
      do iLine=1,nLines
         READ(FileID_ChargeDensity,fmt=*) position, charges

         LatticeIndex = GetLatticeIndex(position)
         if(index_quark==-1.and.all(Charges.ge.0)) then
            index_quark = LatticeIndex
         end if

         if(GetProc_G(LatticeIndex)==ThisProc()) then
            ChargeDensity(:,GetMemoryIndex(LatticeIndex)) = &
                 ChargeDensity(:,GetMemoryIndex(LatticeIndex)) + Charges
         end if
      end do
    end subroutine ReadChargeDensity

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      implicit none

      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureWilsonLines_Equilibrium_StaticMeson
end module programs



program simulation
  use, intrinsic :: iso_fortran_env
  use mpiinterface, only: MPIstop
  use programs
  implicit none

  integer(int8) :: mode

  character(len=3) :: arg

  call get_command_argument(1,arg);
  read(arg,'(I3)') mode

  select case(mode)
  case(2)
     call MeasureEnergyAndGaussLawDeviation
  case (52)
     call MeasureWilsonLines_Equilibrium_StaticMeson
  case default
     call MPIStop('Invalid simulation mode selected')
  end select
end program simulation
