!----------------------------------------------------------------------
! Statistics tools
!----------------------------------------------------------------------
!
! MODULE: statistics
!> @brief
!! Collection of statistics tools
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 01.10.2018
!! @version
!! 1.0
! REVISION HISTORY:
! 01 10 2018 - Initial Version
!----------------------------------------------------------------------
module statistics
  use precision, only: fp
  USE,INTRINSIC :: ISO_FORTRAN_ENV ! defines kinds
  
  implicit none

  PRIVATE

  public :: GetMean,GetStdError,GetStdDeviation
  !> @brief Mean of input
  !! @returns
  !! Mean of input
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 01.10.2018
  !! @version
  !! 1.0
  interface GetMean
     module procedure GetMean_real
     module procedure GetMean_cmplx
  end interface GetMean
contains
  !> @brief Mean of double input
  !! @returns
  !! Mean of double input
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 01.10.2018
  !! @version
  !! 1.0
  pure function GetMean_real(data) result(mean)
    implicit none
    real(fp), intent(in) :: data(:)
    real(fp)             :: mean

    mean = sum(data)/size(data)
  end function GetMean_Real

  !> @brief Mean of double complex input
  !! @returns
  !! Mean of double complex input
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 01.10.2018
  !! @version
  !! 1.0
  pure function GetMean_cmplx(data) result(mean)
    implicit none
    complex(fp), intent(in) :: data(:)
    complex(fp)             :: mean

    mean = sum(data)/size(data)

  end function GetMean_cmplx

  !> @brief Standard error of double input
  !! @returns
  !! Standard error of double input
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 01.10.2018
  !! @version
  !! 1.0
  pure function GetStdError(data) result(stdError)
    implicit none
    real(fp), intent(in) :: data(:)
    real(fp)             :: stdError

    real(fp) :: mean

    stdError = GetStdDeviation(data)/sqrt(real(size(data),fp))

  end function GetStdError

  !> @brief Standard deviation of double input
  !! @returns
  !! Standard deviation of double input
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 01.10.2018
  !! @version
  !! 1.0
  pure function GetStdDeviation(data) result(stdDeviation)
    implicit none
    real(fp), intent(in) :: data(:)
    real(fp)             :: stdDeviation
    
    real(fp) :: mean

    mean     = GetMean(data)

    stdDeviation = sqrt(sum(((data-mean)**2)/(size(data)-1)))

  end function GetStdDeviation
end module statistics
