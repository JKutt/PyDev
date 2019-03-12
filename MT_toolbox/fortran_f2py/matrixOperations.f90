subroutine error_on_output(m, d, g, w, abs_error, sum_error, nc)

      implicit none
    ! input parameter declarations
      complex, dimension(nc, 2), intent(in) :: g
      complex, dimension(nc, 1), intent(in) :: d
      complex, dimension(1, 2), intent(in) :: m
      complex, dimension(nc, 1), intent(in) :: w
      integer, intent(in) :: nc
    ! output parameter declarations
      real, dimension(nc, 1), intent(out) :: abs_error
      real, intent(out) :: sum_error
    ! inside parameter declaration
      complex, dimension(nc, 1) :: error
      complex, dimension(nc, 1) :: abs_error_weighted
      integer :: i
      sum_error = 0.
      do i=1,nc
        error(i, 1) = d(i, 1) - g(i, 1) * m(1, 1) - g(i, 2) * m(1, 2)
        abs_error(i, 1) = abs(error(i, 1) * conjg(error(i, 1)))
        abs_error_weighted(i, 1) = abs(conjg(error(i, 1)) * w(i, 1) * error(i, 1))
        sum_error = sum_error + abs_error_weighted(i, 1) ** 2
      enddo

    end subroutine error_on_output

subroutine prepare_matrices(e, b, br, w, g, gr, d, lhs, rhs, nc)
  ! ****************************
  !   Subroutine prepare_tensor_inversion
  ! ****************************
  ! Algebra to quickly construct matrices for lots of things
  ! Input parameters:
  ! e: 1-D vector containing output variables from the transfer function
  ! b: 2-D vector containing both input variables from the transfer function
  ! br: 2-D vector containing both references variables
  ! w: 1-D vector containing regression weights
  ! Output parameters:
  ! d: data output matrix
  ! g: matrice in regression system (basically b)
  ! gr: another matrice in regression system (basically br)
  ! lhs: left hand side of the regression system
  ! rhs: right hand side of the regression system
  ! *****************************

    implicit none
  ! input parameter declarations
    complex, dimension(nc, 2), intent(in) :: b, br
    complex, dimension(nc, 1), intent(in) :: e
    complex, dimension(nc, 1), intent(in) :: w
    integer, intent(in) :: nc
  ! output parameter declarations
    complex, dimension(nc, 1), intent(out) :: d
    complex, dimension(nc, 2), intent(out) :: g, gr
    complex, dimension(2, 2), intent(out) :: lhs
    complex, dimension(2, 1), intent(out) :: rhs
  ! inside parameter declaration
    complex, dimension(2, nc) :: brt
    complex, dimension(nc, nc) :: weights
    integer :: i

    ! create arrays
    do i=1,nc
      weights(i, :) = (/((cmplx(0, 0)), i=1,nc)/)
      weights(i, i) = cmplx(1, 0) * w(i, 1)
      g(i, 1) = b(i, 1)
      gr(i, 1) = br(i, 1)
      g(i, 2) = b(i, 2)
      gr(i, 2) = br(i, 2)
      d(i, 1) = e(i, 1)
    enddo

    brt = conjg(transpose(br))
    !rhs = matmul(brt, matmul(w, e))
    !lhs = matmul(brt, matmul(w, b))
    rhs = cmplx(0, 0)
    lhs = cmplx(0, 0)

    do i=1,nc
      rhs(1, 1) = rhs(1, 1) + brt(1, i) * w(i, 1) * e(i, 1)
      rhs(2, 1) = rhs(2, 1) + brt(2, i) * w(i, 1) * e(i, 1)
      lhs(1, 1) = lhs(1, 1) + brt(1, i) * w(i, 1) * b(i, 1)
      lhs(2, 1) = lhs(2, 1) + brt(2, i) * w(i, 1) * b(i, 1)
      lhs(1, 2) = lhs(1, 2) + brt(1, i) * w(i, 1) * b(i, 2)
      lhs(2, 2) = lhs(2, 2) + brt(2, i) * w(i, 1) * b(i, 2)
    enddo

  end subroutine prepare_matrices

subroutine prepare_tensor_inversion(e, b, br, w, lhs, rhs, nc)
  ! ****************************
  !   Subroutine prepare_tensor_inversion
  ! ****************************
  ! Algebra to quickly construct matrices before using linear inversion
  ! Input parameters:
  ! e: 1-D vector containing output variables from the transfer function
  ! b: 2-D vector containing both input variables from the transfer function
  ! br: 2-D vector containing both references variables
  ! Output parameters:
  ! lhs: left hand side of the regression system
  ! rhs: right hand side of the regression system
  ! *****************************

    implicit none
  ! input parameter declarations
    complex, dimension(nc, 2), intent(in) :: b, br
    complex, dimension(nc, 1), intent(in) :: e
    complex, dimension(nc, 1), intent(in) :: w
    integer, intent(in) :: nc
  ! output parameter declarations
    complex, dimension(2, 2), intent(out) :: lhs
    complex, dimension(2, 1), intent(out) :: rhs
  ! inside parameter declaration
    complex, dimension(2, nc) :: brt
    complex, dimension(nc, nc) :: weights
    integer :: i

    ! initializing weigths
    do i=1,nc
      weights(i, :) = (/((cmplx(0, 0)), i=1,nc)/)
      weights(i, i) = cmplx(1, 0) * w(i, 1)
    enddo
    brt = conjg(transpose(br))
    !rhs = matmul(brt, matmul(w, e))
    !lhs = matmul(brt, matmul(w, b))
    rhs = cmplx(0, 0)
    lhs = cmplx(0, 0)

    do i=1,nc
      rhs(1, 1) = rhs(1, 1) + brt(1, i) * w(i, 1) * e(i, 1)
      rhs(2, 1) = rhs(2, 1) + brt(2, i) * w(i, 1) * e(i, 1)
      lhs(1, 1) = lhs(1, 1) + brt(1, i) * w(i, 1) * b(i, 1)
      lhs(2, 1) = lhs(2, 1) + brt(2, i) * w(i, 1) * b(i, 1)
      lhs(1, 2) = lhs(1, 2) + brt(1, i) * w(i, 1) * b(i, 2)
      lhs(2, 2) = lhs(2, 2) + brt(2, i) * w(i, 1) * b(i, 2)
    enddo

  end subroutine prepare_tensor_inversion

subroutine dot_product_h(b, u, dot, nc)
  ! ****************************
  !   Subroutine dot_product_h
  ! ****************************
  ! Just a simple dot product with weights
  ! Input parameters:
  ! b: 2D vector containing input predictor
  ! u: N*N matrix whose diagonal contains weights
  ! Output parameters:
  ! dot: dot product with weights
  ! *****************************

  implicit none
! input parameter declarations
  complex, dimension(nc, 2), intent(in) :: b
  integer, intent(in) :: nc
  complex, dimension(nc, nc), intent(in) :: u
  ! output parameter declarations
  complex, dimension(2, 2), intent(out) :: dot
  ! inside parameter declaration
  integer :: i

  do i=1,nc
    dot(1, 1) = dot(1, 1) + conjg(b(i, 1)) * b(i, 1) * u(i, i)
    dot(1, 2) = dot(1, 2) + conjg(b(i, 1)) * b(i, 2) * u(i, i)
    dot(2, 1) = dot(2, 1) + conjg(b(i, 2)) * b(i, 1) * u(i, i)
    dot(2, 2) = dot(2, 2) + conjg(b(i, 2)) * b(i, 2) * u(i, i)
  enddo

end subroutine dot_product_h


subroutine hat_matrix(b, dot, hat, nc)
  ! ****************************
  !   Subroutine hatmatrix
  ! ****************************
  ! Buils the hat matrix of the predictor
  ! Input parameters:
  ! b: 2D vector containing input predictor
  ! dot: dot product with weights
  ! Output parameters:
  ! hat: N*N hat matrix
  ! *****************************

  implicit none
! input parameter declarations
  complex, dimension(nc, 2), intent(in) :: b
  complex, dimension(2, 2), intent(in) :: dot
  integer, intent(in) :: nc
  ! output parameter declarations
  complex, dimension(nc, nc), intent(out) :: hat
  ! inside parameter declaration
  complex, dimension(2, nc) :: array
  integer :: i,j

  do i=1,nc
    array(1, i) = conjg(b(i, 1)) * dot(1, 1) + conjg(b(i, 2)) * dot(1, 2)
    array(2, i) = conjg(b(i, 1)) * dot(2, 1) + conjg(b(i, 2)) * dot(2, 2)
  enddo

  do i=1,nc
    do j=1,nc
      hat(i, j) = b(i, 1) * array(1, j) + b(i, 2) * array(2, j)
    enddo
  enddo

end subroutine hat_matrix
