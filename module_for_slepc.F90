#include <slepc/finclude/slepceps.h>
module module_for_slepc
  use slepceps
  use library_rel_ed  
  implicit none

  type(t_params), pointer :: ptr_dane
  real(8), pointer       :: ptr_przekatna(:)
  integer                :: n_rozmiar

   contains

  subroutine WrapperMatMult(A, x, y, ierr)
    Mat :: A
    Vec :: x, y
    PetscErrorCode :: ierr
    PetscScalar, pointer :: tablica_x(:), tablica_y(:)

    
    call VecGetArrayRead(x, tablica_x, ierr)
    call VecGetArray(y, tablica_y, ierr)

    call matrix_vector_product(ptr_dane, tablica_x, tablica_y)

    call VecRestoreArrayRead(x, tablica_x, ierr)
    call VecRestoreArray(y, tablica_y, ierr)

    ierr = 0
  end subroutine WrapperMatMult

  subroutine WrapperMatGetDiagonal(A, diag, ierr)
    Mat :: A
    Vec :: diag
    PetscErrorCode :: ierr
    PetscScalar, pointer :: tablica_diag(:)
    integer :: i

    call VecGetArray(diag, tablica_diag, ierr)

    do i = 1, n_rozmiar
      tablica_diag(i) = ptr_przekatna(i)
    end do

    call VecRestoreArray(diag, tablica_diag, ierr)

    ierr = 0
  end subroutine WrapperMatGetDiagonal

end module module_for_slepc