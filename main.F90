#include <slepc/finclude/slepceps.h>
program main
use slepceps
use library_rel_ed
use module_for_slepc
implicit none
type(t_params), target, save :: dane, dane_el_gf_a
integer:: n_alpha,n_beta,norb
integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
logical:: relativistic
integer :: i, j,k,l,p,q 
real(8) :: r
complex(8), allocatable :: hopping_alpha(:,:),hopping_beta(:,:),interaction_alpha(:,:,:,:),interaction_beta(:,:,:,:)
complex(8), allocatable :: interaction_temp(:,:)
real, allocatable :: orbital_energies(:)
complex(8), allocatable, target, save :: diagonal(:)
complex(8), allocatable :: ground_state(:),new_state(:)
integer :: krylov_size
integer :: N
integer              :: nconv
integer :: liczba_wartosci = 1
PetscScalar          :: wartosc_rzeczywista, wartosc_urojona
Vec                  :: wektor_petsc, wektor_urojony_petsc
PetscScalar, pointer :: tablica_wynikowa(:)
real :: start_time, end_time
real :: total_time
PetscErrorCode :: ierr
Mat            :: A
EPS            :: eps
call cpu_time(start_time)
call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
  !********* INPUT **********
relativistic=.true.
dane%relativistic = relativistic


norb=6

n_alpha=1
n_beta=1
  
n_RAS_spaces_occ=0
n_RAS_spaces_virt=0

allocate(RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt))
allocate(excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt))
! number of orbitals cannot be equal 0

!RAS_space_occ(1)=2
!RAS_space_occ(2)=0
active_space=norb
!RAS_space_virt(1)=2
!RAS_space_virt(2)=0
excit_array(:)=1
call check_orbital_space_declarations(norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)


allocate(orbital_energies(norb))
allocate(interaction_temp(norb**2,norb**2))
allocate(hopping_alpha(norb,norb))
allocate(hopping_beta(norb,norb))
allocate(dane%interaction_mix(norb,norb,norb,norb))
allocate(interaction_alpha(norb,norb,norb,norb))
allocate(interaction_beta(norb,norb,norb,norb))
allocate(dane%hso_ab(norb,norb))
allocate(dane%hso_ba(norb,norb))
open(unit = 1, file = "integrals/mo_energies.txt")
open(unit = 2, file = "integrals/h_upup.txt")
open(unit = 4, file = "integrals/h_dndn.txt")
open(unit = 3, file = "integrals/eri_mo.txt")
open(unit = 5, file = "integrals/h_updn.txt")
open(unit = 7, file = "integrals/h_dnup.txt")
do i=1,norb 
    read(1,*) orbital_energies(i)
    read(2,*) hopping_alpha(i,:)
    read(4,*) hopping_beta(i,:)
    read(5,*) dane%hso_ab(i,:)
    read(7,*) dane%hso_ba(i,:)
end do

do i=1,norb**2
    read(3,*) interaction_temp(i,:)
end do

p = 1
do i=1,norb
    do j=1,norb
        q = 1
        do k=1,norb
            do l=1,norb 
                interaction_alpha(i,j,k,l) = interaction_temp(p,q)
                q = q + 1
            end do
        end do
        p = p + 1
    end do
end do

interaction_beta = interaction_alpha
dane%interaction_mix = interaction_alpha

call check_hamiltonian(relativistic,norb,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,dane%interaction_mix,dane%hso_ab,dane%hso_ba)

call generate_dane(dane,norb,n_alpha,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta)

diagonal = dane%diag

N = dane%size_tot(1,2) + dane%size_tot(2,2) + dane%size_tot(3,2)
allocate(ground_state(N))
ptr_dane => dane
ptr_przekatna => diagonal


call MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, &
                      N, N, PETSC_NULL_INTEGER, A, ierr)

call MatShellSetOperation(A, MATOP_MULT, WrapperMatMult, ierr)
call MatShellSetOperation(A, MATOP_GET_DIAGONAL, WrapperMatGetDiagonal, ierr)

call EPSCreate(PETSC_COMM_WORLD, eps, ierr)
call EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr)

call EPSSetProblemType(eps, EPS_HEP, ierr)

call EPSSetType(eps,EPSJD, ierr)
  
  
call EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr)
call EPSSetDimensions(eps, liczba_wartosci, PETSC_DECIDE, PETSC_DECIDE, ierr)
call EPSSetFromOptions(eps, ierr)
  
  
call EPSSolve(eps, ierr)
if (ierr .ne. 0) then
    print *, "BŁĄD: EPSSolve nie powiodło się! Kod:", ierr
    stop
end if
call EPSGetConverged(eps, nconv, ierr)
print *, "Solver found ", nconv, " convergent eigenvalues."
call MatCreateVecs(A, wektor_petsc, wektor_urojony_petsc, ierr)

do i = 0,nconv - 1
    
    call EPSGetEigenpair(eps, i, wartosc_rzeczywista, wartosc_urojona, &
                         wektor_petsc, wektor_urojony_petsc, ierr)

    print *, "Eigenvalue no.", i + 1, "=", real(wartosc_rzeczywista)+0.7151043390810812



  end do
    
call EPSGetEigenpair(eps, 0, wartosc_rzeczywista, wartosc_urojona, &
                         wektor_petsc, wektor_urojony_petsc, ierr)




call VecGetArrayRead(wektor_petsc, tablica_wynikowa, ierr)
ground_state = tablica_wynikowa
call VecRestoreArrayRead(wektor_petsc, tablica_wynikowa, ierr)

call VecDestroy(wektor_petsc, ierr)
call VecDestroy(wektor_urojony_petsc, ierr)



call EPSDestroy(eps, ierr)
call MatDestroy(A, ierr)

call SlepcFinalize(ierr)

nullify(ptr_dane)
nullify(ptr_przekatna)
r=0
do i=1,N 
write(*,*)i, ground_state(i)
r=r+abs(ground_state(i))**2
end do
write(*,*) r

!HERE WE START GENERATING ELECTRONIC GREEN FUNCTION
!BLOCK n_alpha+1
krylov_size = 4
dane_el_gf_a%relativistic = relativistic
dane_el_gf_a%interaction_mix = interaction_alpha
call check_orbital_space_declarations(norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha+1,n_beta)
call generate_dane(dane_el_gf_a,norb,n_alpha+1,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta)
allocate(new_state(dane_el_gf_a%size_tot(1,2)+dane_el_gf_a%size_tot(2,2)+dane_el_gf_a%size_tot(3,2)))
call create_electron(1,ground_state,new_state,dane,dane_el_gf_a,1)





write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time, dane%size_tot(1,2),dane%size_tot(2,2),dane%size_tot(3,2),dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)
call flush(6)



end program main