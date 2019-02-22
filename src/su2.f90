!------------------------------------------------------------------------------
! SU(2)-gauge group
!------------------------------------------------------------------------------
!
! MODULE: su2
!> @brief
!! SU(2) gauge group
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 22.02.2019
!! @version
!! 1.2
! REVISION HISTORY:
! 14 07 2018 - Initial Version
! 11 01 2019 - Added functions
! 22 02 2019 - Using precision-module
!------------------------------------------------------------------------------
module su2
  use, intrinsic :: iso_fortran_env
  use precision, only: fp

  implicit none

  PRIVATE
  
  public :: &
       GetGenerator,&
       GetGroupExp, GetGroupLog,&
       GetWrappedAlgebraCoordinates,&
       GetAlgebraMatrix,&
       GetAlgebraCoordinate,&
       GetTraceWithGenerator
  
  !> Matrix size of the SU(N) gauge group elements
  integer(int8),parameter,public :: nsun = 2
  !> Number of generators of the gauge group
  integer(int8),parameter,public :: ngen = nsun**2-1

  !> @brief SU(2)-generators aka Pauli-matrices
  !! @details
  !! Pauli-matrices times \f$1/2\f$ following the convention
  !! \f$\Tr(T^a\cdot T^b)=2\delta_{a,b}\f$:
  !! \f{align*}{
  !! T_0=\begin{pmatrix}
  !! +1 & 0 \\
  !! 0  & +1 
  !! \end{pmatrix},&&
  !! \\
  !! T_1=\frac{1}{2}\sigma_1=\frac{1}{2}\cdot\begin{pmatrix}
  !! 0  & +1  \\
  !! +1 & 0 
  !! \end{pmatrix},&&
  !! T_2=\frac{1}{2}\sigma_2=\frac{1}{2}\cdot\begin{pmatrix}
  !! 0         & -\text{i} \\
  !! +\text{i} & 0 
  !! \end{pmatrix},&&
  !! T_3=\frac{1}{2}\sigma_3=\frac{1}{2}\cdot\begin{pmatrix}
  !! +1 & 0  \\
  !! 0 & -1
  !! \end{pmatrix}
  !! \f}
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date
  !! 14.07.2018
  !! @version
  !! 1.0
  complex(fp), dimension(nsun,nsun,1:ngen), parameter, public :: &
       Generators = reshape(&
       [&                      ! First generator
       cmplx(00,00,fp),        cmplx(+0.5_fp,00,fp), & 
       cmplx(+0.5_fp,00,fp),   cmplx(00,00,fp)&
       ,&                      ! Second generator
       & cmplx(00,00,fp),      cmplx(00,-0.5_fp,fp),& 
       cmplx(00,+0.5_fp,fp),   cmplx(00,00,fp)&
       ,&                      ! Third generator
       & cmplx(+0.5_fp,00,fp), cmplx(00,00,fp),& 
       cmplx(00,00,fp),        cmplx(-0.5_fp,00,fp)&
       ],shape=[nsun,nsun,ngen])

    !> @brief Exponential from Lie-algebra to Lie-group
  !! @returns
  !! Matrix exponential from su(N)-Lie-algebra to SU(N)-Lie-group
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 14.07.2018
  !! @version 1.0
  interface GetGroupExp
     module procedure GetGroupExp_fromAlgebraCoordinates
     module procedure GetGroupExp_fromAlgebraMatrix
  end interface GetGroupExp
contains
  !> @brief
  !! Returns the i'th generator
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 14.07.2018
  !! @version 1.0
  !! @returns
  !! i'th generator of the SU(N)-Lie-group
  pure function GetGenerator(i) result(res)
    implicit none
    !> Generator index
    integer(int8), intent(in) :: i
    !> Generator
    complex(fp), dimension(nsun,nsun) :: res

    if(i.lt.1 .or. i>ngen) then
       res = -1
    else
       res = Generators(:,:,i)
    end if
  end function GetGenerator
  
  ! ..--** Auxiliary Mathematical Routines **--..

  !> @brief
  !! Wraps \f$\mathfrak{su}(N)\f$-Lie-algebra parameters back to first branch
  !! @details
  !! The wrapping is performed is such a way that the result still corresponds to the same
  !! Lie-group element. Because this is specific to the considered gauge group, it has to be
  !! done in the following way:
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.08.2018
  !! @version 1.0
  pure function GetWrappedAlgebraCoordinates(input)
    use mathconstants, only: pi
    use matrixoperations, only: EigenH, GetHerm
    implicit none
    !> \f$\mathfrak{su}(N)\f$-Lie-algebra parameters
    real(fp), dimension(ngen), intent(in) :: input
    real(fp), dimension(ngen)             :: GetWrappedAlgebraCoordinates

    ! su(N)-Lie-algebra element
    complex(fp), dimension(nsun,nsun) :: algebra
    ! Similarity transform to diagonalise Lie-algebra element
    complex(fp), dimension(nsun,nsun) :: S,hermS,algebraS
    ! Eigenvalues of Lie-algebra element
    real(fp), dimension(nsun) :: ev
    ! Lie-group generator index
    integer(int8) :: a

    ! Get input times generators
    algebra = GetAlgebraMatrix(input)

    ! Diagonalize that
    call EigenH(algebra,ev,S)

    ! Wraps to [-pi,+pi[
    ev = modulo(ev+pi,2*pi)-pi

    ! Similarity transform back
    algebra = 0
    forall(a=1:nsun) algebra(a,a) = ev(a)

    hermS    = GetHerm(S)
    algebraS = matmul(algebra,S)
    algebra  = matmul(hermS,algebraS)

    ! Get projected group parameters via scalar product with generator
    forall(a=1:ngen) GetWrappedAlgebraCoordinates(a) = GetAlgebraCoordinate(a,algebra)
  end function GetWrappedAlgebraCoordinates

  !> @brief
  !! Returns matrix-exponential with real
  !! \f$\mathfrak{su}(N)\f$-coordinates
  !! @details
  !! Computes \f$\exp(\text{i}\alpha_a T^a)\f$
  !! for real input \f$\vec{\alpha}\f$
  !! @returns
  !! Exponential \f$\exp\left(\text{i}\alpha_a T^a\right)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.08.2018
  !! @version 1.0
  pure function GetGroupExp_fromAlgebraCoordinates(input) result(res)
    implicit none
    !> \f$\mathfrak{su}(N)\f$-coordinates \f$\alpha_a\f$
    !! in (hermitian) Lie-algebra
    real(fp),    dimension(ngen),     intent(in) :: input
    !> SU(N)-Lie-group element corresponding to input
    complex(fp), dimension(nsun,nsun)            :: res

    complex(fp), dimension(nsun,nsun) :: input_times_generators

    input_times_generators = GetAlgebraMatrix(input)
    res = GetGroupExp_fromAlgebraMatrix(input_times_generators)
  end function GetGroupExp_fromAlgebraCoordinates

  !> @brief
  !! Returns matrix-exponential for given
  !! \f$\mathfrak{su}(N)\f$-Lie-algebra element
  !! @returns
  !! Exponential \f$\exp\left(\text{i} A\right)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.08.2018
  !! @version 1.0
  pure function GetGroupExp_fromAlgebraMatrix(input) result(res)
    use matrixoperations, only: ExpAH
    implicit none
    !> (Hermitian) \f$\mathfrak{su}(N)\f$-Lie-algebra element \f$A\f$
    complex(fp),dimension(nsun,nsun),intent(in) :: input
    !> SU(N)-Lie-group element corresponding to input
    complex(fp),dimension(nsun,nsun)            :: res

    ! skew-hermitian matrix constructed from input
    complex(fp),dimension(nsun,nsun)            :: A

    A = input*cmplx(0._fp,1._fp,fp)

    res = ExpAH(A)
  end function GetGroupExp_fromAlgebraMatrix

  !> @brief
  !! Returns
  !! \f$\mathfrak{su}(N)\f$-Lie-algebra element \f$A=\alpha^aT^a\f$
  !! @returns
  !! \f$A=\alpha^aT^a\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.08.2018
  !! @version 1.0
  pure function GetAlgebraMatrix(AlgebraCoordinates) result(res)
    implicit none
    !> \f$\mathfrak{su}(N)\f$-Lie-algebra coordinates \f$\alpha^a\f$
    real(fp),    dimension(ngen),     intent(in) :: AlgebraCoordinates
    !> \f$\mathfrak{su}(N)\f$-Lie-algebra element \f$A\f$
    complex(fp), dimension(nsun,nsun)            :: res

    ! Generator index
    integer(int8) :: a

    res = 0._fp
    do a=1,ngen
       res = res + AlgebraCoordinates(a)*Generators(:,:,a)
    end do !a
  end function GetAlgebraMatrix

  !> @brief
  !! Returns matrix-logarithm for given
  !! \f$\mathfrak{su}(N)\f$-Lie-algebra element to SU(N)-Lie-group
  !! @details
  !! The output is the hermitian \f$\mathfrak{su}(N)\f$-Lie-algebra
  !! element \f$A\f$
  !! instead of the direct logarithm \f$\log\exp(iA)=iA\f$
  !! @returns
  !! Logarithm \f$\log_{\text{SU(N)}}(G)=\log\exp(iA)=iA\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.08.2018
  !! @version 1.0
  pure function GetGroupLog(input) result(res)
    use matrixoperations, only: ArgLogU
    implicit none
    !> SU(N)-Lie-group element \f$G\f$
    complex(fp), dimension(nsun,nsun), intent(in) :: input
    !> \f$\mathfrak{su}(N)\f$-Lie-algebra element
    real(fp),    dimension(ngen)                  :: res

    ! Logarithm of input
    complex(fp), dimension(nsun,nsun) :: LogInput
    ! Lie-group generator index
    integer(int8) :: a

    LogInput = ArgLogU(input)

    forall(a=1:ngen) res(a) = GetAlgebraCoordinate(a,LogInput)
  end function GetGroupLog

  !> @brief
  !! Returns \f$2\Re\Tr(T^a\cdot M)\f$
  !! @details
  !! This function is the scalar product with a generator
  !! and thus gives the group parameter
  !! @returns
  !! \f$2\Re\Tr(T^a\cdot M)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.08.2018
  !! @version 1.0
  pure real(fp) function GetAlgebraCoordinate(a,matrix)
    use matrixoperations, only: GetTrace
    implicit none
    !> Generator index
    integer(int8),         intent(in) :: a
    !> Input matrix \f$M\f$
    complex(fp), intent(in) :: matrix(nsun,nsun)
    GetAlgebraCoordinate = 2._fp*real(GetTraceWithGenerator(a,matrix),fp)
  end function GetAlgebraCoordinate
  
  !> @brief
  !! Returns \f$\Tr(T^a\cdot M)\f$
  !! @returns
  !! \f$\Tr(T^a\cdot M)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 14.01.2019
  !! @version 1.0
  pure complex(fp) function GetTraceWithGenerator(a,matrix)
    use matrixoperations, only: GetTrace
    implicit none
    !> Generator index
    integer(int8),         intent(in) :: a
    !> Input matrix \f$M\f$
    complex(fp), intent(in) :: matrix(nsun,nsun)
    
    complex(fp) :: generator_times_matrix(nsun,nsun)

    generator_times_matrix = matmul(Generators(:,:,a),matrix)
    GetTraceWithGenerator  = GetTrace(generator_times_matrix)
  end function GetTraceWithGenerator

  !> @brief
  !! Kronecker-Delta
  !! @returns
  !! Kronecker-Delta \f$\delta_{i,j}\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 26.11.2018
  !! @version 1.0
  pure elemental integer(int8) function GetKroneckerDelta(i,j)
    implicit none
    !> Index \f$i\f$
    integer(int8), intent(in) :: i
    !> Index \f$j\f$
    integer(int8), intent(in) :: j

    if(i==j) then
       GetKroneckerDelta = 1
    else
       GetKroneckerDelta = 0
    end if
  end function GetKroneckerDelta
end module su2
