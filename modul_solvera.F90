#include <slepc/finclude/slepceps.h>
module modul_solvera
  use slepceps
  use library_rel_ed  ! Podłączasz swój własny moduł z rutynami i typami
  implicit none

  ! 1. Definiujemy wskaźniki. Dzięki nim wrappery "zobaczą" dane z programu main
  ! Zastąp 'TwojTyp' faktyczną nazwą typu z library_rel_ed
  type(t_params), pointer :: ptr_dane
  real(8), pointer       :: ptr_przekatna(:)
  integer                :: n_rozmiar

   contains

  ! 2. Wrapper dla mnożenia macierz-wektor
  subroutine WrapperMatMult(A, x, y, ierr)
    Mat :: A
    Vec :: x, y
    PetscErrorCode :: ierr
    PetscScalar, pointer :: tablica_x(:), tablica_y(:)

    ! Pobranie tablic z wektorów PETSc
    call VecGetArrayRead(x, tablica_x, ierr)
    call VecGetArray(y, tablica_y, ierr)

    ! Wywołanie TWOJEJ rutyny z modułu library_rel_ed.
    ! Przekazujemy wskaźnik na Twoje dane oraz surowe tablice
    call matrix_vector_product(ptr_dane, tablica_x, tablica_y)

    ! Oddanie tablic do PETSc
    call VecRestoreArrayRead(x, tablica_x, ierr)
    call VecRestoreArray(y, tablica_y, ierr)

    ierr = 0
  end subroutine WrapperMatMult

  ! 3. Wrapper dla pobierania wektora przekątnej
  subroutine WrapperMatGetDiagonal(A, diag, ierr)
    Mat :: A
    Vec :: diag
    PetscErrorCode :: ierr
    PetscScalar, pointer :: tablica_diag(:)
    integer :: i

    call VecGetArray(diag, tablica_diag, ierr)

    ! Kopiujemy wartości z Twojego wektora (widzianego przez wskaźnik) do PETSc
    do i = 1, n_rozmiar
      tablica_diag(i) = ptr_przekatna(i)
    end do

    call VecRestoreArray(diag, tablica_diag, ierr)

    ierr = 0
  end subroutine WrapperMatGetDiagonal

end module modul_solvera