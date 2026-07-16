#include <slepc/finclude/slepceps.h>
program main
use slepceps
use library_rel_ed
use module_for_slepc
implicit none

!********HERE ARE INPUT VARIABLES********

integer:: n_alpha,n_beta,n_orb
integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
logical:: relativistic
complex(8), allocatable :: hopping_alpha(:,:),hopping_beta(:,:),interaction_alpha(:,:,:,:),interaction_beta(:,:,:,:),interaction_mix(:,:,:,:), hso_ab(:,:), hso_ba(:,:)
real(8), allocatable :: orbital_energies(:)
integer :: krylov_size
integer :: n_energies = 1
real(8) :: nuclear_energy
!*********HERE ARE OTHER VARIABLES*********

type(parameters), target, save :: gs_params, gf_params
complex(8), allocatable :: interaction_temp(:,:)
integer :: i,j,k,l,p,q
complex(8), allocatable, target, save :: diagonal(:)
complex(8), allocatable :: ground_state(:),new_state(:)
integer :: N
integer              :: nconv
PetscScalar          :: real_part, imaginary_part
Vec                  :: vector_petsc, im_vector_petsc
PetscScalar, pointer :: tablica_wynikowa(:)
real :: start_time, end_time
real :: total_time
PetscErrorCode :: ierr
Mat            :: A
EPS            :: eps
call cpu_time(start_time)
call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)

!*******HERE USER DEFINES VARIABLES********

relativistic=.true.
n_orb = 6
n_alpha = 1
n_beta = 1

n_RAS_spaces_occ=0
n_RAS_spaces_virt=0

!----------------------------DO NOT TOUCH-------------------------------------
allocate(RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt))
allocate(excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt))
!-----------------------  YOU CAN TOUCH NOW-------------------------------------------

!RAS_space_occ(1)=2
!RAS_space_occ(2)=0
active_space=n_orb
!RAS_space_virt(1)=2
!RAS_space_virt(2)=0
excit_array(:)=1

!---------------------------DO NOT TOUCH--------------------------------------------
call check_orbital_space_declarations(n_orb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)

allocate(orbital_energies(n_orb))
allocate(hopping_alpha(n_orb,n_orb))
allocate(hopping_beta(n_orb,n_orb))
allocate(interaction_mix(n_orb,n_orb,n_orb,n_orb))
allocate(interaction_alpha(n_orb,n_orb,n_orb,n_orb))
allocate(interaction_beta(n_orb,n_orb,n_orb,n_orb))
allocate(hso_ab(n_orb,n_orb))
allocate(hso_ba(n_orb,n_orb))
!-------------------------YOU CAN TOUCH NOW-------------------------------------------

!******************HERE USER NEEDS TO PROVIDE MATRIX ELEMENTS FOR DIFFERENT PARTS OF HAMILTONIAN**********************

open(unit = 1, file = "integrals/mo_energies.txt")
open(unit = 2, file = "integrals/h_upup.txt")
open(unit = 4, file = "integrals/h_dndn.txt")
open(unit = 3, file = "integrals/eri_mo.txt")
open(unit = 5, file = "integrals/h_updn.txt")
open(unit = 7, file = "integrals/h_dnup.txt")
open(unit = 8, file = "integrals/nuc_energy.txt")
do i=1,n_orb 
    read(1,*) orbital_energies(i)
    read(2,*) hopping_alpha(i,:)
    read(4,*) hopping_beta(i,:)
    read(5,*) hso_ab(i,:)
    read(7,*) hso_ba(i,:)
end do
    read(8,*) nuclear_energy
allocate(interaction_temp(n_orb**2,n_orb**2))
do i=1,n_orb**2
    read(3,*) interaction_temp(i,:)
end do 

p = 1
do i=1,n_orb
    do j=1,n_orb
        q = 1
        do k=1,n_orb
            do l=1,n_orb 
                interaction_alpha(i,j,k,l) = interaction_temp(p,q)
                q = q + 1
            end do
        end do
        p = p + 1
    end do
end do

interaction_beta = interaction_alpha
interaction_mix = interaction_alpha


call check_hamiltonian(relativistic,n_orb,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix,hso_ab,hso_ba)
call generate_params(gs_params,relativistic,n_orb,n_alpha,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta, interaction_mix, hso_ab, hso_ba)
diagonal = gs_params%diag

N = gs_params%size_tot(1,2) + gs_params%size_tot(2,2) + gs_params%size_tot(3,2)
allocate(ground_state(N))
ptr_dane => gs_params
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
call EPSSetDimensions(eps, n_energies, PETSC_DECIDE, PETSC_DECIDE, ierr)
call EPSSetFromOptions(eps, ierr)
  
  
call EPSSolve(eps, ierr)
if (ierr .ne. 0) then
    print *, "BŁĄD: EPSSolve nie powiodło się! Kod:", ierr
    stop
end if
call EPSGetConverged(eps, nconv, ierr)
print *, "Solver found ", nconv, " convergent eigenvalues."
call MatCreateVecs(A, vector_petsc, im_vector_petsc, ierr)

do i = 0,nconv - 1
    
    call EPSGetEigenpair(eps, i, real_part, imaginary_part, &
                         vector_petsc, im_vector_petsc, ierr)

    print *, "Eigenvalue no.", i + 1, "=", real(real_part) + nuclear_energy


  end do
    
call EPSGetEigenpair(eps, 0, real_part, imaginary_part, &
                         vector_petsc, im_vector_petsc, ierr)




call VecGetArrayRead(vector_petsc, tablica_wynikowa, ierr)
ground_state = tablica_wynikowa
call VecRestoreArrayRead(vector_petsc, tablica_wynikowa, ierr)

call VecDestroy(vector_petsc, ierr)
call VecDestroy(im_vector_petsc, ierr)



call EPSDestroy(eps, ierr)
call MatDestroy(A, ierr)

call SlepcFinalize(ierr)

nullify(ptr_dane)
nullify(ptr_przekatna)


!HERE WE START GENERATING ELECTRONIC GREEN FUNCTION
!BLOCK n_alpha+1
krylov_size = 4


call generate_params(gf_params,relativistic,n_orb,n_alpha+1,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix, hso_ab, hso_ba)
allocate(new_state(gf_params%size_tot(1,2)+gf_params%size_tot(2,2)+gf_params%size_tot(3,2)))
call create_electron(1,ground_state,new_state,gs_params,gf_params,1)





write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time, gs_params%size_tot(1,2),gs_params%size_tot(2,2),gs_params%size_tot(3,2),gs_params%size_tot(1,2)+gs_params%size_tot(2,2)+gs_params%size_tot(3,2)
call flush(6)



end program main