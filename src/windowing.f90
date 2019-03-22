module windowing
  USE,INTRINSIC :: ISO_FORTRAN_ENV ! defines kinds
  use precision, only: fp

  PRIVATE

  public ApplyHannWindow, ApplyStepWindow

contains

  pure subroutine ApplyHannWindow(data)
    use mathconstants, only: pi
    complex(fp), intent(inout) :: data(:)

    integer :: size_data, i

    size_data = size(data)

    ! data*sin(t/(T/2))**2
    ! with
    ! ... T = signal length in time domain
    do i=1,size_data
       data(i) = data(i) * sin(pi*i/size_data)**2
    end do
  end subroutine ApplyHannWindow

  pure subroutine ApplyStepWindow(data,cut_index)
    use mathconstants, only: pi
    complex(fp), intent(inout) :: data(:)
    integer, intent(in) :: cut_index

    integer :: size_data, i

    size_data = size(data)

    data(cut_index:size(data)) = 0
  end subroutine ApplyStepWindow
end module windowing
