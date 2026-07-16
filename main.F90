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
integer :: n_pairs = 2
real(8) :: nuclear_energy
!*********HERE ARE OTHER VARIABLES*********

type(parameters), target, save :: gs_params, gf_params
complex(8), allocatable :: interaction_temp(:,:)
integer :: i,j,k,l,p,q
!complex(8), allocatable, target, save :: diagonal(:)
complex(8), allocatable :: ground_state(:),new_state(:)
real(8), allocatable :: eigenenergies(:)
complex(8), allocatable :: eigenstates(:,:)
integer :: N


real :: start_time, end_time
real :: total_time

call cpu_time(start_time)


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


call generate_params(gs_params,relativistic,n_orb,n_alpha,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta, interaction_mix, hso_ab, hso_ba)
!diagonal = gs_params%diag
gs_params%nuclear_energy = nuclear_energy
N = gs_params%size_tot(1,2) + gs_params%size_tot(2,2) + gs_params%size_tot(3,2)
allocate(eigenenergies(n_pairs))
allocate(eigenstates(n_pairs,N))
allocate(ground_state(N))
call eigensystem(gs_params,n_pairs,eigenenergies,eigenstates)


!HERE WE START GENERATING ELECTRONIC GREEN FUNCTION
!BLOCK n_alpha+1
krylov_size = 4


call generate_params(gf_params,relativistic,n_orb,n_alpha+1,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix, hso_ab, hso_ba)
allocate(new_state(gf_params%size_tot(1,2)+gf_params%size_tot(2,2)+gf_params%size_tot(3,2)))
!call create_electron(1,ground_state,new_state,gs_params,gf_params,1)





write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time, gs_params%size_tot(1,2),gs_params%size_tot(2,2),gs_params%size_tot(3,2),gs_params%size_tot(1,2)+gs_params%size_tot(2,2)+gs_params%size_tot(3,2)
call flush(6)



end program main