program simulation
  use, intrinsic :: iso_fortran_env
  use precision
  use mpiinterface, only: ThisProc, NumProcs, MPIstop, intmpi, syncall
  use lattice, only: ndim
  use lattice
  use halocomm
  use xpfft
  use gaugeconfiguration_su3
  use mpi
  use io
  use mathconstants
  use tolerances
  implicit none

  ! Simulation parameters
  integer(int64) :: LatticeExtensions(ndim)!(ndim)
  real(fp)       :: LatticeSpacings(0:ndim)!(0:ndim)
  integer(int64) :: TimeSteps
  integer(int64) :: RandomNumberSeed

  real(fp)   :: GluonSaturationScale !qs
  real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
  real(fp)   :: GluonCoupling
  integer(int64) :: EnsembleSize
  real(fp)   :: Quarkmass
  real(fp)   :: CoMTime
  real(fp)   :: TimeRange
  !real(fp)  :: Wilsoncoeffs(nWilsonCoeffs)

  real(fp) :: GaugefixingCoefficient
  integer(int64) :: Number_of_Measurements_of_Gluondistribution
  integer(int64) :: TimePoints_between_Measurement_of_Gluondistribution
  real(fp) :: TimeBetweenGluonMeasurements
  
  ! Physical fields
  type(GaugeConfiguration) :: GaugeConf
  
  complex(fp), allocatable :: testcomm(:)

  ! Monitoring variables
  integer(int64) :: it, idistmeas

  integer(int64) :: MemoryIndex, LatticeIndex, i, is
  integer(int64), allocatable :: LocalLatticeIndices(:)
  
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
  
  ! Control observable output (Gauss, Energy)
  integer(int8) :: fileID_eg
  real(fp) :: gauss, energy
  real(fp) :: time

  complex(fp), allocatable :: data(:)
  real(fp) :: Momentum(ndim)
  integer(int64) :: MomentumIndices(ndim)
  
  call InitSimulation
  
  call GaugeConf%TransversePolarisedOccupiedInit_Box(&
       GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)

  !.................................
  !.... Measuring gluon distribution
  !..
  allocate(AA_correlator_ensemble(GetLocalLatticeSize(),&
       0:Number_of_Measurements_of_Gluondistribution,EnsembleSize))
  allocate(EE_correlator_ensemble(GetLocalLatticeSize(),&
       0:Number_of_Measurements_of_Gluondistribution,EnsembleSize))

  ensemble: do iensemble=1,EnsembleSize
     if(ThisProc()==0) print*,iensemble

     call GaugeConf%TransversePolarisedOccupiedInit_Box(&
          GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling &
          !,aa_correlator_opt,ee_correlator_opt)
          )
     idistmeas = 0
     call GaugeConf%GetTransverseAACorrelator(aa_correlator_opt)
     call GaugeConf%GetTransverseEECorrelator(ee_correlator_opt)
     
     AA_correlator_ensemble(:,idistmeas,iensemble) = aa_correlator_opt
     EE_correlator_ensemble(:,idistmeas,iensemble) = ee_correlator_opt
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

  allocate(GluonDistribution_AtTimePoint(GetLocalLatticeSize()))
  allocate(GluonDistribution_AtTimePoint_stderror(GetLocalLatticeSize()))
  allocate(Momenta(GetLocalLatticeSize()))
  is=0
  allocate(mpisendrequest(3))
  do src=0,NumProcs()-1
     if(ThisProc()==dest .or. ThisProc()==src) then

        ! Momenta
        buffersize = size(Momenta)
        tag = (0+NumProcs())*src
        if(ThisProc()==src) then
           do MemoryIndex=1,GetMemorySize()
              LatticeIndex = GetLatticeIndex(MemoryIndex)
              if(ThisProc()==GetProc(LatticeIndex)) then
                 is = is + 1
                 Momenta(is) = GetNorm2Momentum(LatticeIndex)
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

  call endsimulation



  ! Test of fourier transform
  allocate(data(getmemorysize()))

  ! Set data to signal of a sine
  MomentumIndices=([1,1,2]-1)
  
  Momentum = 2*pi*MomentumIndices/GetLatticeExtension([1_int8:nDim])
  do MemoryIndex=1,GetMemorySize()
     LatticeIndex = GetLatticeIndex(MemoryIndex)
     if(ThisProc()==GetProc(LatticeIndex)) then
        data(MemoryIndex) &
             != Exp(Cmplx(0,sum(Momentum*(GetLatticePosition(LatticeIndex)-1)),fp))
             = sin(sum(Momentum*(GetLatticePosition(LatticeIndex)-1)))
     end if
  end do

  call x2p(data)

  do LatticeIndex=1,GetLatticeSize()

     if(ThisProc()==GetProc(LatticeIndex)) then
        MemoryIndex = GetMemoryIndex(LatticeIndex)

        if(abs(data(MemoryIndex)).gt.GetZeroTol()) &
             print*,int(GetLatticePosition(LatticeIndex),int16),':',data(MemoryIndex)
     end if

     call Flush(Output_Unit)
     call SyncAll
  end do


  call endsimulation



  !.......................
  !.... Energy measurement
  !..
  !call GaugeConf%HotInit
  !if(ThisProc()==0) &
  !     fileID_eg = OpenFile(filename="energy_gauss.txt",st='REPLACE',fm='FORMATTED',act='WRITE')
  !it=0
  !time = it*GetLatticeSpacing(0_int8)
  !gauss = GaugeConf%GetDeviationFromGaussLaw()
  !energy= GaugeConf%GetEnergy()
  !if(ThisProc()==0) write(fileID_eg,'(3(SP,E13.6,1X))') time,energy,gauss

  !do it=1,TimeSteps
  !   call GaugeConf%Update
  !   time = it*GetLatticeSpacing(0_int8)
  !   gauss = GaugeConf%GetDeviationFromGaussLaw()
  !   energy= GaugeConf%GetEnergy()
  !   if(ThisProc()==0) write(fileID_eg,'(3(SP,E13.6,1X))') time,energy,gauss
  !   if(ThisProc()==0) write(output_unit,*) time
  !end do
  !if(ThisProc()==0) call CloseFile(fileID_eg)


  
  
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
    implicit none

    integer(int64) :: arg_count
    character(len=80) :: arg
    integer(int8) :: i
    
    !..--** Reading simulation parameters **--..
    arg_count = 0
    
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
end program simulation
