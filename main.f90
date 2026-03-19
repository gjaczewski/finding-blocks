program relativistic_ed

use library_rel_ed
implicit none
type(t_params), target :: dane
integer:: n_alpha,n_beta,norb
integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space, n_spaces
integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
integer, dimension(:,:), allocatable:: RAS_el_array_alpha,RAS_el_array_beta
logical:: relativistic
integer :: NRD_spin_alpha, NRD_spin_beta
integer :: i, j, temp, n_distributions_alpha,n_distributions_beta, n_distributions_alpha_p1,n_distributions_beta_p1, n_distributions_alpha_m1,n_distributions_beta_m1
integer :: max,n_combinations
integer, allocatable :: all_combinations(:,:)
integer:: NRDa(3,2),NRDb(3,2),sizea(3,2),sizeb(3,2)
complex, allocatable :: hopping_alpha(:,:),hopping_beta(:,:),interaction_alpha(:,:,:,:),interaction_beta(:,:,:,:),vector(:),vector_new(:)
real, allocatable :: diagonal(:)
real :: start_time, end_time
real :: total_time

call cpu_time(start_time)

  !********* INPUT **********
relativistic=.false.
  
norb=2

n_alpha=1
n_beta=1
  
n_RAS_spaces_occ=0
n_RAS_spaces_virt=0

allocate(RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt))
allocate(excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt))
! number of orbitals cannot be equal 0

!RAS_space_occ(1)=2
!RAS_space_occ(2)=0
active_space=2
!RAS_space_virt(1)=2
!RAS_space_virt(2)=0
excit_array(:)=1
  
allocate(hopping_alpha(norb,norb))
allocate(hopping_beta(norb,norb))
allocate(dane%interaction_mix(norb,norb,norb,norb))
allocate(interaction_alpha(norb,norb,norb,norb))
allocate(interaction_beta(norb,norb,norb,norb))
hopping_alpha(:,:) = 0
hopping_alpha(1,2) = -1
hopping_alpha(2,1) = -1
hopping_beta = hopping_alpha
interaction_alpha(:,:,:,:) = 0
interaction_beta(:,:,:,:) = 0
dane%interaction_mix(:,:,:,:) = 0
dane%interaction_mix(1,1,1,1) = 5
dane%interaction_mix(2,2,2,2) = 5
allocate(dane%hso_ab(norb,norb))
allocate(dane%hso_ba(norb,norb))
if (relativistic .eqv. .true.) then
    dane%hso_ab(:,:) = 0
    dane%hso_ab(1,2) = 11
    dane%hso_ab(2,1) = 11
    dane%hso_ab(1,1) = 23
    dane%hso_ab(2,2) = 23
    dane%hso_ba = dane%hso_ab
end if 
!***********************
!here we check if input is ok 

call check_orbital_space_declarations(norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)
call check_hamiltonian(relativistic,norb,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,dane%interaction_mix,dane%hso_ab,dane%hso_ba)

! calculate all combinations of assigning electrons to spaces

max = MAXVAL([RAS_space_occ,active_space,RAS_space_virt])
n_spaces = n_RAS_spaces_occ+1+n_RAS_spaces_virt
n_combinations = (max+1)**n_spaces
allocate(all_combinations(n_combinations,n_spaces))

do i=1,n_combinations
    temp = i 
    do j=1,n_spaces
        all_combinations(i,j) = mod(temp,max+1)
        temp = temp/(max+1)
    end do
end do

! here we calculate number of possible distributions NRD_spin_alpha, NRD_spin_beta and number of distributions for fixed spin

call count_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta)

allocate(RAS_el_array_alpha(NRD_spin_alpha,n_spaces))
allocate(RAS_el_array_beta(NRD_spin_beta,n_spaces))

write(*,*) "NRD_spin_alpha: ",NRD_spin_alpha
write(*,*) "NRD_spin_beta: ",NRD_spin_beta

! here we fill RAS_el_array_spin tables and calculate NRDa, NRDb
!NRD(:,1) denotes beginning of spin block, NRD(:,2) denotes end

call find_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,NRDa,NRDb)

write(*,*) "NRDa:"
write(*,*) NRDa(1,:)
write(*,*) NRDa(2,:)
write(*,*) NRDa(3,:)

write(*,*) "NRDb: "
write(*,*) NRDb(1,:)
write(*,*) NRDb(2,:)
write(*,*) NRDb(3,:)

! here we calculate sizea,sizeb,dane%size_tot

call calculate_space_size(relativistic,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,dane%size_tot)

write(*,*) "Space size alpha: ", sizea(1,2), sizea(2,2), sizea(3,2)
write(*,*) "Space size beta: ", sizeb(1,2), sizeb(2,2), sizeb(3,2)

write(*,*) "Total size:"
write(*,*) dane%size_tot(1,:)
write(*,*) dane%size_tot(2,:)
write(*,*) dane%size_tot(3,:)
  
  !dane%size_tot(1:3,1:2) - an array of pointers to the blocks of the vector (or Hamiltonian)
  
  !order of blocks: (n_alpha, n_beta); (n_alpha+1, n_beta-1); (n_alpha-1, n_beta+1)
  
  !dane%size_tot(1,:) - for sector (n_alpha, n_beta)
  !dane%size_tot(2,:) - for sector (n_alpha+1, n_beta-1)
  !dane%size_tot(3,:) - for sector (n_alpha-1, n_beta+1) 

  !dane%size_tot(:,1) - the beginning of such sector
  !dane%size_tot(:,2) - the size (number of elements) of/in such sector

!Here we check if any crucial blocks are empty

if (dane%size_tot(1,2) .eq. 0) then
    write(*,*) "MAIN BLOCK IS EMPTY"
    write(*,*) "QUITTING THE PROGRAM"
    !SHOULD WE INLCUDE VACUUM?
    stop
else if ((relativistic .eqv. .true.) .and. (dane%size_tot(2,2)+dane%size_tot(3,2) .eq. 0)) then
    write(*,*) "RELATIVISTIC BLOCKS ARE EMPTY"
    write(*,*) "QUITTING THE PROGRAM"
    stop
end if

! here we generate strings in first block (alpha and beta separately)

if (n_alpha .gt. 0)  then
  allocate(dane%str_a(sizea(1,2),n_alpha))
  call fill_spin_strings(n_alpha,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(1,1),NRDa(1,2),sizea(1,2),dane%str_a)
else if (n_alpha .eq. 0)  then
  allocate(dane%str_a(sizea(1,2),1))
  dane%str_a(sizea(1,2),:) = 0 !vacuum  
end if

if (n_beta .gt. 0) then
  allocate(dane%str_b(sizeb(1,2),n_beta))
  call fill_spin_strings(n_beta,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(1,1),NRDb(1,2),sizeb(1,2),dane%str_b)
else if  (n_beta .eq. 0)  then
  allocate(dane%str_b(sizeb(1,2),1))
  dane%str_b(sizeb(1,2),:) = 0 !vacuum  
end if

!here we generate creation-annihilation matrices for first block

allocate(dane%alpha_annihilation_creation_matrix(norb,norb,sizea(1,2),2))  
call fill_annihilation_creation_matrix(n_alpha,sizea(1,2),dane%str_a,norb,dane%alpha_annihilation_creation_matrix)

allocate(dane%beta_annihilation_creation_matrix(norb,norb,sizeb(1,2),2))
call fill_annihilation_creation_matrix(n_beta,sizeb(1,2),dane%str_b,norb,dane%beta_annihilation_creation_matrix)

!here we generate single spin part of hamiltonian for first block

allocate(dane%alpha_hamiltonian(sizea(1,2),sizea(1,2)))
allocate(dane%beta_hamiltonian(sizeb(1,2),sizeb(1,2)))

call fill_nr_single_spin_hamiltonian(dane%str_a,sizea(1,2),n_alpha,norb,hopping_alpha,interaction_alpha,dane%alpha_annihilation_creation_matrix,dane%alpha_hamiltonian)
call fill_nr_single_spin_hamiltonian(dane%str_b,sizeb(1,2),n_beta,norb,hopping_beta,interaction_beta,dane%beta_annihilation_creation_matrix,dane%beta_hamiltonian)

! here we prepare for matrix-vector multiplication

allocate(vector(dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)))
allocate(vector_new(dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)))
allocate(diagonal(dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)))

call generate_diagonal_elements(dane%alpha_hamiltonian,dane%str_a,sizea(1,2),n_alpha,dane%beta_hamiltonian,dane%str_b,sizeb(1,2),n_beta,dane%interaction_mix,norb,diagonal(1:dane%size_tot(1,2)))

!here we repeat for other blocks



!here we generate strings
    
    if (n_alpha .gt. 1)  then
        allocate(dane%str_a_m1(sizea(3,2),n_alpha-1))
        if (sizea(3,2) .ne. 0) then
            call fill_spin_strings(n_alpha-1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(3,1),NRDa(3,2),sizea(3,2),dane%str_a_m1)
        end if
    else if (n_alpha .eq. 1) then
        allocate(dane%str_a_m1(sizea(3,2),1))
        dane%str_a_m1(:,:) = 0 !vacuum
    else
        allocate(dane%str_a_m1(sizea(3,2),0))
    end if
      
    if (n_beta .gt. 1)  then
        allocate(dane%str_b_m1(sizeb(3,2),n_beta-1))
        if (sizeb(3,2) .ne. 0) then
            call fill_spin_strings(n_beta-1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(3,1),NRDb(3,2),sizeb(3,2),dane%str_b_m1)
        end if
    else if (n_beta .eq. 1) then
        allocate(dane%str_b_m1(sizeb(3,2),1))
        dane%str_b_m1(:,:) = 0 !vacuum
    else 
        allocate(dane%str_b_m1(sizeb(3,2),0))
    end if
    
    allocate(dane%str_a_p1(sizea(2,2),n_alpha+1))
    if (sizea(2,2) .ne. 0) then
     call fill_spin_strings(n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),dane%str_a_p1)
    end if
      
    allocate(dane%str_b_p1(sizeb(2,2),n_beta+1)) 
    if (sizeb(2,2) .ne. 0) then
        call fill_spin_strings(n_beta+1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(2,1),NRDb(2,2),sizeb(2,2),dane%str_b_p1)
    end if

!here we generate annihilation, creation-annihilation and single spin hamiltonians
      
    allocate(dane%alpha_annihilation_matrix(norb,sizea(1,2),2))
    allocate(dane%alpha_annihilation_creation_matrix_m1(norb,norb,sizea(3,2),2))
    allocate(dane%alpha_hamiltonian_m1(sizea(3,2),sizea(3,2)))

    if (sizea(3,2) .ne. 0) then      
        call fill_annihilation_results(sizea(3,2),sizea(1,2),n_alpha,dane%str_a_m1,dane%str_a,norb,dane%alpha_annihilation_matrix)
        call fill_annihilation_creation_matrix(n_alpha-1,sizea(3,2),dane%str_a_m1,norb,dane%alpha_annihilation_creation_matrix_m1)
        call fill_nr_single_spin_hamiltonian(dane%str_a_m1,sizea(3,2),n_alpha-1,norb,hopping_alpha,interaction_alpha,dane%alpha_annihilation_creation_matrix_m1,dane%alpha_hamiltonian_m1)
    end if
    
    allocate(dane%beta_annihilation_matrix(norb,sizeb(1,2),2))
    allocate(dane%beta_annihilation_creation_matrix_m1(norb,norb,sizeb(3,2),2))
    allocate(dane%beta_hamiltonian_m1(sizeb(3,2),sizeb(3,2)))
    
    if (sizeb(3,2) .ne. 0) then      
        call fill_annihilation_results(sizeb(3,2),sizeb(1,2),n_beta,dane%str_b_m1,dane%str_b,norb,dane%beta_annihilation_matrix)
        call fill_annihilation_creation_matrix(n_beta-1,sizeb(3,2),dane%str_b_m1,norb,dane%beta_annihilation_creation_matrix_m1)
        call fill_nr_single_spin_hamiltonian(dane%str_b_m1,sizeb(3,2),n_beta-1,norb,hopping_beta,interaction_beta,dane%beta_annihilation_creation_matrix_m1,dane%beta_hamiltonian_m1)
    end if
 

    allocate(dane%alpha_annihilation_matrix_p1(norb,sizea(2,2),2))
    allocate(dane%alpha_annihilation_creation_matrix_p1(norb,norb,sizea(2,2),2))
    allocate(dane%alpha_hamiltonian_p1(sizea(2,2),sizea(2,2)))
    
    if (sizea(2,2) .ne. 0) then
        call fill_spin_strings(n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),dane%str_a_p1)
        call fill_annihilation_results(sizea(1,2),sizea(2,2),n_alpha+1,dane%str_a,dane%str_a_p1,norb,dane%alpha_annihilation_matrix_p1)
        call fill_annihilation_creation_matrix(n_alpha+1,sizea(2,2),dane%str_a_p1,norb,dane%alpha_annihilation_creation_matrix_p1)
        call fill_nr_single_spin_hamiltonian(dane%str_a_p1,sizea(2,2),n_alpha+1,norb,hopping_alpha,interaction_alpha,dane%alpha_annihilation_creation_matrix_p1,dane%alpha_hamiltonian_p1)
    end if

    allocate(dane%beta_annihilation_matrix_p1(norb,sizeb(2,2),2))  
    allocate(dane%beta_annihilation_creation_matrix_p1(norb,norb,sizeb(2,2),2))
    allocate(dane%beta_hamiltonian_p1(sizeb(2,2),sizeb(2,2)))   

    if (sizeb(2,2) .ne. 0) then
        call fill_spin_strings(n_beta+1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(2,1),NRDb(2,2),sizeb(2,2),dane%str_b_p1)
        call fill_annihilation_results(sizeb(1,2),sizeb(2,2),n_beta+1,dane%str_b,dane%str_b_p1,norb,dane%beta_annihilation_matrix_p1)
        call fill_annihilation_creation_matrix(n_beta+1,sizeb(2,2),dane%str_b_p1,norb,dane%beta_annihilation_creation_matrix_p1)
        call fill_nr_single_spin_hamiltonian(dane%str_b_p1,sizeb(2,2),n_beta+1,norb,hopping_beta,interaction_beta,dane%beta_annihilation_creation_matrix_p1,dane%beta_hamiltonian_p1)
    end if

! here we generate hamiltonian diagonal elements
    if (dane%size_tot(2,2) .ne. 0) then
        call generate_diagonal_elements(dane%alpha_hamiltonian_p1,dane%str_a_p1,sizea(2,2),n_alpha+1,dane%beta_hamiltonian_m1,dane%str_b_m1,sizeb(3,2),n_beta-1,dane%interaction_mix,norb,diagonal(dane%size_tot(2,1):dane%size_tot(2,2)+dane%size_tot(1,2)))
    end if

    if (dane%size_tot(3,2) .ne. 0) then 
        call generate_diagonal_elements(dane%alpha_hamiltonian_m1,dane%str_a_m1,sizea(3,2),n_alpha-1,dane%beta_hamiltonian_p1,dane%str_b_p1,sizeb(2,2),n_beta+1,dane%interaction_mix,norb,diagonal(dane%size_tot(3,1):dane%size_tot(3,2)+dane%size_tot(2,2)+dane%size_tot(1,2)))
    end if
dane%n_strings_alpha = sizea(1,2)
dane%n_strings_beta = sizeb(1,2)
dane%norb = norb
dane%n_alpha = n_alpha
dane%n_beta = n_beta
dane%n_strings_alpha_p1 = sizea(2,2)
dane%n_strings_beta_p1 = sizeb(2,2)
dane%n_strings_alpha_m1 = sizea(3,2)
dane%n_strings_beta_m1 = sizeb(3,2)
dane%relativistic = relativistic

vector = [0,0,0,1]
vector_new(:) = 0
call matrix_vector_product(dane,vector,vector_new)
write(*,*) vector_new




write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time, dane%size_tot(1,2),dane%size_tot(2,2),dane%size_tot(3,2),dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)
call flush(6)


end program relativistic_ed