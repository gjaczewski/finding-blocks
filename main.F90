#include <slepc/finclude/slepceps.h>
program main
use slepceps
use library_rel_ed
implicit none

!********HERE ARE INPUT VARIABLES********

integer:: n_alpha,n_beta,n_orb
integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
logical:: relativistic
complex(8), allocatable :: hopping_alpha(:,:),hopping_beta(:,:),interaction_alpha(:,:,:,:),interaction_beta(:,:,:,:),interaction_mix(:,:,:,:), hso_ab(:,:), hso_ba(:,:)
real(8), allocatable :: orbital_energies(:)
integer :: krylov_size
integer :: n_pairs = 5

real(8) :: nuclear_energy
!*********HERE ARE OTHER VARIABLES*********
PetscErrorCode :: ierr
type(parameters), target, save :: gs_params, params2
complex(8), allocatable :: gf1(:,:), gf2(:,:)
complex(8), allocatable :: interaction_temp(:,:)
integer :: i,j,k,l,p,q
!complex(8), allocatable, target, save :: diagonal(:)
complex(8), allocatable :: ground_state(:), new_states(:,:)
real(8), allocatable :: eigenenergies(:)
complex(8), allocatable :: eigenstates(:,:)
integer :: N
real(8) :: gs_energy
complex(8) :: omega
real :: start_time, end_time
real :: total_time
complex(8) :: r1,r2
call cpu_time(start_time)


!*******HERE USER DEFINES VARIABLES********

relativistic=.false.
n_orb = 2
n_alpha = 1
n_beta = 1

n_RAS_spaces_occ=0
n_RAS_spaces_virt=0

!----------------------------DO NOT TOUCH-------------------------------------
call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
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
interaction_alpha(:,:,:,:) = 0
interaction_beta = interaction_alpha
interaction_mix = interaction_alpha
interaction_mix(1,1,1,1) = 4
interaction_mix(2,2,2,2) = 4
nuclear_energy = 0 
hopping_alpha(1,1)=5
hopping_alpha(2,2)=5
hopping_alpha(1,2)=-1
hopping_alpha(2,1)=-1
hopping_beta = hopping_alpha

call generate_params(gs_params,relativistic,n_orb,n_alpha,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta, interaction_mix, hso_ab, hso_ba)
gs_params%nuclear_energy = nuclear_energy
N = gs_params%size_tot(1,2) + gs_params%size_tot(2,2) + gs_params%size_tot(3,2)
allocate(eigenenergies(n_pairs))
allocate(eigenstates(n_pairs,N))
allocate(ground_state(N))
call eigensystem(gs_params,n_pairs,eigenenergies,eigenstates,ierr)
ground_state = eigenstates(1,:)
gs_energy = eigenenergies(1)

deallocate(eigenenergies)

deallocate(eigenstates)

call generate_params(params2,relativistic,n_orb,n_alpha+1,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta, interaction_mix, hso_ab, hso_ba)


params2%nuclear_energy = nuclear_energy
N = params2%size_tot(1,2) + params2%size_tot(2,2) + params2%size_tot(3,2)
allocate(eigenenergies(2))
allocate(eigenstates(2,N))
call eigensystem(params2,2,eigenenergies,eigenstates,ierr)
allocate(gf1(n_orb,n_orb))
allocate(gf2(n_orb,n_orb))
!HERE WE START GENERATING ELECTRONIC GREEN FUNCTION
!BLOCK n_alpha+1
krylov_size = 1000
!call generate_params(gf_params,relativistic,n_orb,n_alpha+1,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix, hso_ab, hso_ba)
gf1(:,:) = 0
gf2(:,:) = 0
allocate(new_states(2,2))

omega = 1
!call create_electron(1,[(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(1.0d0,0.0d0)],new_state1,gs_params,params2,1)
do i=1,n_orb
    call create_electron(i,ground_state,new_states(i,:),gs_params,params2,1)
end do


do i=1,n_orb
    do j=1,n_orb
        do k=1,2
            r1=0
            r2=0
            do l=1,2
                r1 = r1 + conjg(new_states(i,l))*eigenstates(k,l)
                r2 = r2 + conjg(eigenstates(k,l))*new_states(j,l)
            end do
            gf2(i,j) = gf2(i,j) + r1*r2/(omega-eigenenergies(k)+gs_energy)
        end do
    end do
end do

write(*,*) "GF:"
write(*,*) gf2(1,1)
write(*,*) gf2(2,2)
write(*,*) gf2(1,2)
write(*,*) gf2(2,1)
write(*,*) "*********"

do i=1,n_orb
    do j=1,n_orb
        call calc_e_gf(i,j,1,1,params2,params2,gs_params,ground_state,gs_energy,omega,krylov_size,gf1(i,j))
    end do
end do

write(*,*) gf1(1,1)
write(*,*) gf1(2,2)
write(*,*) gf1(1,2)
write(*,*) gf1(2,1)
write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time
call flush(6)


call SlepcFinalize(ierr)
end program main