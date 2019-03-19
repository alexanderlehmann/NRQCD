!-------------------------------------------------------------------------
! PROGRAMS for Lattice-NRQCD
!-------------------------------------------------------------------------
!
! MODULE: programs
!>@brief Program/mains
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 07.03.2019
!!@version 1.0
! REVISION HISTORY:
! 07 03 2019 - Initial Version
!-------------------------------------------------------------------------
module programs
  PUBLIC

contains
  !>@brief Program for measuring the heavy quarkonium correlator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine MeasureHeavyQuarkoniumCorrelators
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    use nrqcd

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: CoMTime
    real(fp) :: TimeRange
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_atT1
    type(NRQCDField)         :: HeavyField, HeavyField_atT1

    ! Counting
    integer :: i

    ! CMS coordinates
    real(fp) :: t,s,t1,t2
    integer(int64) :: is, it

    ! Output
    real(fp) :: norm_quark, norm_antiq
    complex(fp) :: mesoncorrelator
    integer(int8) :: FileID_Norm, FileID_Correlator
    
    
    call InitSimulation
    
    if(ThisProc()==0) then
       fileID_Correlator = OpenFile(filename="Correlator_3s1.txt",&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_Norm = OpenFile(filename="Norm.txt",&
            st='REPLACE',fm='FORMATTED',act='WRITE')
    end if
    
    ! Determination of CMS coordinates
    t  = CoMTime   ! fixed
    s  = TimeRange ! varies
    t1 = t-s/2     ! varies
    t2 = t+s/2     ! varies
    
    ! Initial value for gauge field at t=0
    call GaugeConf_atT1%ColdInit

    ! Evolving gauge configuration to t1 (possibly negative)
    TimeSteps = abs(NINT(t1/LatticeSpacings(0)))
    do it=1,TimeSteps
       if(ThisProc()==0) write(output_unit,*) int(it,int16),'of',int(TimeSteps,int16)
       call GaugeConf_atT1%Update(sign(+1._real64,t1))
    end do

    ! Initial value for heavy field at t1
    call HeavyField_atT1%InitSinglePoint(spin=1_int8,colour=1_int8,latticeindex=1_int64)
    svalues: do is=-nint(TimeRange/LatticeSpacings(0)),nint(+TimeRange/LatticeSpacings(0))-1
       s = is*LatticeSpacings(0)

       if(ThisProc()==0) write(output_unit,*) 's=',real(s,real32),'step',int(is,int16),'max=',&
            int(nint(+TimeRange/LatticeSpacings(0)),int16)

       ! Setting links for time-evolution at t1
       GaugeConf = GaugeConf_atT1

       ! Initialising quark-pair at t1 ...
       HeavyField = HeavyField_atT1

       ! ... and evolving to t2
       TimeSteps = abs(is)
       dtstep_in_s: do it=1,TimeSteps
          if(thisproc()==0) write(output_unit,*) int(it,int16),'of',int(TimeSteps,int16)
          call HeavyField%Update(GaugeConf,HeavyQuarkMass,WilsonCoefficients,real(sign(+1,is),fp))
          call GaugeConf%Update
       end do dtstep_in_s

       mesoncorrelator = HeavyField%GetMesonCorrelator_3s1_ZeroMomentum()
       norm_quark = HeavyField%GetNorm_Quark()
       norm_antiq = HeavyField%GetNorm_AntiQ()

       if(ThisProc()==0) then
          write(FileID_Correlator,'(3(SP,E16.9,1X))') &
               s,real(mesoncorrelator,real64),aimag(mesoncorrelator)
          write(FileID_Norm,'(3(SP,E19.12,1X))')&
             s,norm_quark,norm_antiq
       end if

       call GaugeConf_atT1%Update
    end do svalues

    ! Closing files
    if(ThisProc()==0) then
       call CloseFile(FileID_Correlator)
       call CloseFile(FileID_Norm)
    end if
    
    call HeavyField%Destructor
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
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      use NRQCD,              only: InitModule_NRQCD              => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      real(fp) :: c_re, c_im
      real(fp) :: kspTol
      
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

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      ! Tolerance for iterative PETSc solver
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') kspTol
      
      ! Heavy quark mass
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') HeavyQuarkmass
      
      ! Wilson coefficents
      do i=1,nWilsonCoefficients
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_re
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_im

         WilsonCoefficients(i) = cmplx(c_re,c_im,fp)
      end do

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances
      call InitModule_NRQCD(HeavyQuarkMass,WilsonCoefficients,1._fp,kspTol)

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      use NRQCD,          only: FinalizeModule_NRQCD          => FinalizeModule
      implicit none

      call FinalizeModule_NRQCD
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      STOP "Simulation completed"
    end subroutine EndSimulation
  end subroutine MeasureHeavyQuarkoniumCorrelators
  
  !>@brief Program for measuring the time-evolution of the energy and gauss law deviation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
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

    real(fp)   :: GluonSaturationScale !qs
    real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp)   :: GluonCoupling
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

    call GaugeConf%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
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

    call endsimulation
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
      use xpfft,              only: InitModule_xpFFT              => InitModule
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

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
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
      use mpiinterface, only: FinalizeModule_MPIinterface => FinalizeModule
      use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
      implicit none

      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      STOP "Simulation completed"
    end subroutine EndSimulation
  end subroutine MeasureEnergyAndGaussLawDeviation

  !>@brief Program for measuring the gluon distribution function
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine MeasureGluondistribution
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

    real(fp)   :: GluonSaturationScale !qs
    real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp)   :: GluonCoupling
    integer(int64) :: EnsembleSize
    real(fp)   :: CoMTime
    real(fp)   :: TimeRange

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf

    ! Monitoring variables
    integer(int64) :: it, idistmeas
    real(fp) :: time

    integer(int64) :: MemoryIndex, LatticeIndex, i, is

    ! MPI
    integer(intmpi) :: mpierr, mpistatus(mpi_status_size), tag, src, dest, buffersize

    real(fp), allocatable :: aa_correlator_opt(:),ee_correlator_opt(:)
    ! Ensemble variables
    real(fp), allocatable :: AA_correlator_ensemble(:,:,:), EE_correlator_ensemble(:,:,:)
    real(fp), allocatable :: AA_correlator(:,:), AA_correlator_stderror(:,:)
    real(fp), allocatable :: EE_correlator(:,:), EE_correlator_stderror(:,:)
    real(fp), allocatable :: GluonDistribution(:,:), GluonDistribution_stderror(:,:)
    real(fp), allocatable :: &
         GluonDistribution_AtTimePoint(:), GluonDistribution_AtTimePoint_stderror(:)
    integer :: iensemble
    real(fp), allocatable :: momenta(:)

    integer(intmpi) :: proc

    ! Output
    integer(int8) :: FileID
    character(len=128) :: filename,time_tag
    integer(intmpi), allocatable :: mpisendrequest(:), mpisendstatus(:,:)
    
    integer(int64) :: GaugeFixingMaxIterations
    real(fp) :: GaugefixingCoefficient, GaugeFixingTolerance
    integer(int64) :: Number_of_Measurements_of_Gluondistribution
    integer(int64) :: TimePoints_between_Measurement_of_Gluondistribution
    real(fp) :: TimeBetweenGluonMeasurements

    type(GaugeConfiguration) :: gaugeconf_gaugefixed

    call InitSimulation
    
    allocate(AA_correlator_ensemble(GetLocalLatticeSize(),&
         0:Number_of_Measurements_of_Gluondistribution,EnsembleSize))
    allocate(EE_correlator_ensemble(GetLocalLatticeSize(),&
         0:Number_of_Measurements_of_Gluondistribution,EnsembleSize))

    ensemble: do iensemble=1,EnsembleSize
       if(ThisProc()==0) write(output_unit,*)&
            int(iensemble,int16),'of',&
            int(EnsembleSize,int16),'configurations';&
            call flush(output_unit)

       !call GaugeConf%SemicoldInit
       call GaugeConf%TransversePolarisedOccupiedInit_Box(&
            GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling&
       !,aa_correlator_opt,ee_correlator_opt&
            )
       TimeEvolution: do it=0,TimeSteps

          if(modulo(it,TimePoints_between_Measurement_of_Gluondistribution)==0) then
             if(ThisProc()==0) write(output_unit,*)&
                  it,'of',&
                  TimeSteps,'time steps';&
                  call flush(output_unit)

             idistmeas = it/TimePoints_between_Measurement_of_Gluondistribution

             gaugeconf_gaugefixed = GaugeConf
             call gaugeconf_gaugefixed%CoulombGaugefixing(&
                  Tolerance_    =GaugeFixingTolerance,&
                  alpha_        =GaugefixingCoefficient,&
                  MaxIterations_=GaugefixingMaxIterations)
             call gaugeconf_gaugefixed%GetTransverseAACorrelator(aa_correlator_opt)
             call gaugeconf_gaugefixed%GetTransverseEECorrelator(ee_correlator_opt)

             AA_correlator_ensemble(:,idistmeas,iensemble) = aa_correlator_opt
             EE_correlator_ensemble(:,idistmeas,iensemble) = ee_correlator_opt

          end if

          call GaugeConf%Update
       end do TimeEvolution
    end do ensemble

    ! ..--** Start: Statistics and output **--..
    it = 0
    idistmeas = 0
    call StatisticsCorrelator(AA_correlator_ensemble,AA_correlator,AA_correlator_stderror,0)
    call StatisticsCorrelator(EE_correlator_ensemble,EE_correlator,EE_correlator_stderror,0)

    allocate(GluonDistribution(&
         lbound(AA_correlator,1):ubound(AA_correlator,1),&
         lbound(AA_correlator,2):ubound(AA_correlator,2)))
    allocate(GluonDistribution_stderror(&
         lbound(AA_correlator,1):ubound(AA_correlator,1),&
         lbound(AA_correlator,2):ubound(AA_correlator,2)))

    GluonDistribution = &
                                !AA_correlator
         sqrt(AA_correlator*EE_correlator)
    GluonDistribution_stderror = &
                                !AA_correlator_stderror
         sqrt(&
         + (AA_correlator_stderror*EE_correlator)**2 &
         + (EE_correlator_stderror*AA_correlator)**2&
         )

    ! Printing to file
    dest=0

    allocate(GluonDistribution_AtTimePoint(GetLocalLatticeSize()))
    allocate(GluonDistribution_AtTimePoint_stderror(GetLocalLatticeSize()))
    allocate(Momenta(GetLocalLatticeSize()))
    is=0
    allocate(mpisendrequest(3))
    do idistmeas=lbound(GluonDistribution,2),ubound(GluonDistribution,2)
       if(ThisProc()==dest) then
          write(filename,"(A18,I0.3,A1,I0.3,A1,I0.3)") 'gluondistribution_', &
               LatticeExtensions(1),'x',LatticeExtensions(2),'x',&
               LatticeExtensions(3)

          time = (idistmeas - lbound(GluonDistribution,2)) &
               * TimeBetweenGluonMeasurements
          write(time_tag,"(F12.3)") time
          time_tag = '_t' // trim(ADJUSTL(time_tag))
          filename = trim(filename) // trim(time_tag) // '.txt'

          fileID = OpenFile(filename=filename,st='REPLACE',fm='FORMATTED',act='WRITE')
       end if
       do src=0,NumProcs()-1
          if(ThisProc()==dest .or. ThisProc()==src) then

             ! Momenta
             buffersize = size(Momenta)
             tag = (0+NumProcs())*src
             if(ThisProc()==src) then
                is = 0
                do MemoryIndex=1,GetMemorySize()
                   LatticeIndex = GetLatticeIndex_M(MemoryIndex)
                   if(ThisProc()==GetProc_G(LatticeIndex)) then
                      is = is + 1
                      Momenta(is) = GetNorm2Momentum_G(LatticeIndex)
                   end if
                end do
                call mpi_isend(&
                     Momenta,                      & ! What to send
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     dest,                         & ! Recieving process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpisendrequest(1),            & ! Request handle
                     mpierr)                         ! Error code
             end if
             if(ThisProc()==dest) then
                call mpi_recv(&
                     Momenta,                      & ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if

             ! Gluon distribution
             buffersize = size(GluonDistribution_atTimePoint)
             tag = (1+NumProcs())*src
             if(ThisProc()==src) then
                GluonDistribution_atTimePoint = GluonDistribution(:,idistmeas)
                call mpi_isend(&
                     GluonDistribution_atTimePoint,& ! What to send
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     dest,                         & ! Recieving process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpisendrequest(2),            & ! Request handle
                     mpierr)                         ! Error code
             end if
             if(ThisProc()==dest) then
                call mpi_recv(&
                     GluonDistribution_atTimePoint,& ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if

             ! Standard error of gluon distribution
             buffersize = size(GluonDistribution_atTimePoint)
             tag = (2+NumProcs())*src
             if(ThisProc()==src) then
                GluonDistribution_atTimePoint_stderror = GluonDistribution_stderror(:,idistmeas)
                call mpi_isend(&
                     GluonDistribution_atTimePoint_stderror,& ! What to send
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     dest,                         & ! Recieving process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpisendrequest(3),            & ! Request handle
                     mpierr)                         ! Error code
             end if
             if(ThisProc()==dest) then
                call mpi_recv(&
                     GluonDistribution_atTimePoint_stderror,& ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if


             if(ThisProc()==dest) then
                do is=1,size(Momenta)
                   write(fileID,'(3(SP,E13.6,1X))') &
                        Momenta(is),&
                        GluonDistribution_AtTimePoint(is),&
                        GluonDistribution_AtTimePoint_stderror(is)
                end do
             end if
          end if

          call SyncAll
       end do
       if(ThisProc()==0) call CloseFile(FileID)
    end do
    call endsimulation
  contains
    pure subroutine StatisticsCorrelator(&
         correlator_ensemble,correlator,correlator_stderror,lbound_time)
      use precision, only: fp
      use statistics, only: GetMean, GetStdError

      implicit none
      real(fp),intent(in) :: correlator_ensemble(:,lbound_time:,:)
      real(fp),allocatable,intent(out):: correlator(:,:)
      real(fp),allocatable,intent(out):: correlator_stderror(:,:)
      integer, intent(in) :: lbound_time

      integer :: latticeindex, it

      allocate(correlator(&
           lbound(correlator_ensemble,1):ubound(correlator_ensemble,1),&
           lbound(correlator_ensemble,2):ubound(correlator_ensemble,2)))

      allocate(correlator_stderror(&
           lbound(correlator_ensemble,1):ubound(correlator_ensemble,1),&
           lbound(correlator_ensemble,2):ubound(correlator_ensemble,2)))

      forall(&
           it          =lbound(correlator_ensemble,2):ubound(correlator_ensemble,2),&
           latticeindex=lbound(correlator_ensemble,1):ubound(correlator_ensemble,1))&
           correlator(latticeindex,it) = GetMean(correlator_ensemble(latticeindex,it,:))

      if(size(correlator_ensemble,3)>1) then
         forall(&
              it          =lbound(correlator_ensemble,2):ubound(correlator_ensemble,2),&
              latticeindex=lbound(correlator_ensemble,1):ubound(correlator_ensemble,1))&
              correlator_stderror(latticeindex,it) = GetStdError(correlator_ensemble(latticeindex,it,:))
      else
         correlator_stderror = 0
      end if
    end subroutine StatisticsCorrelator

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
      use xpfft,              only: InitModule_xpFFT              => InitModule
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

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      ! Ensemble size (statistical average of classical statistical simulation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') EnsembleSize

      ! Time between gluon measurements
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeBetweenGluonMeasurements
      Number_of_Measurements_of_Gluondistribution = aint(TimeRange/TimeBetweenGluonMeasurements,int64)
      TimePoints_between_Measurement_of_Gluondistribution&
           = int(TimeBetweenGluonMeasurements/LatticeSpacings(0))

      ! Maximum number of iterations for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I7)') GaugefixingMaxIterations
      ! Tolerance for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') GaugefixingTolerance
      ! Tolerance for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F5.7)') GaugefixingCoefficient

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
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
      use mpiinterface, only: FinalizeModule_MPIinterface => FinalizeModule
      use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
      implicit none
      
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      STOP "Simulation completed"
    end subroutine EndSimulation
  end subroutine MeasureGluondistribution
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
  case(1)
     call MeasureGluondistribution
  case(2)
     call MeasureEnergyAndGaussLawDeviation
  case(3)
     call MeasureHeavyQuarkoniumCorrelators
  case default
     call MPIStop('Invalid simulation mode selected')
  end select
end program simulation
