program relativistic_ed
implicit none
integer:: n_alpha,n_beta,norb
integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space, n_spaces
integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
integer, dimension(:,:), allocatable:: RAS_el_array_alpha,RAS_el_array_beta
logical:: relativistic
integer :: NRD_spin_alpha, NRD_spin_beta
integer :: i, j, temp, n_distributions_alpha,n_distributions_beta, n_distributions_alpha_p1,n_distributions_beta_p1, n_distributions_alpha_m1,n_distributions_beta_m1
integer :: max,n_combinations
integer, allocatable :: all_combinations(:,:)
integer:: NRDa(3,2),NRDb(3,2),sizea(3,2),sizeb(3,2),size_tot(3,2)
integer, dimension(:,:), allocatable :: str_a, str_b, str_a_p1, str_a_m1, str_b_p1, str_b_m1
integer, allocatable :: alpha_annihilation_matrix(:,:,:), alpha_annihilation_matrix_p1(:,:,:), beta_annihilation_matrix(:,:,:), beta_annihilation_matrix_p1(:,:,:)
integer, allocatable :: alpha_annihilation_creation_matrix(:,:,:,:), beta_annihilation_creation_matrix(:,:,:,:),alpha_annihilation_creation_matrix_p1(:,:,:,:), beta_annihilation_creation_matrix_p1(:,:,:,:),alpha_annihilation_creation_matrix_m1(:,:,:,:), beta_annihilation_creation_matrix_m1(:,:,:,:)
complex, allocatable :: alpha_hamiltonian(:,:), beta_hamiltonian(:,:),alpha_hamiltonian_p1(:,:), beta_hamiltonian_p1(:,:),alpha_hamiltonian_m1(:,:), beta_hamiltonian_m1(:,:)
complex, allocatable :: hopping_alpha(:,:),hopping_beta(:,:), interaction_mix(:,:,:,:),interaction_alpha(:,:,:,:),interaction_beta(:,:,:,:),hso_ab(:,:),hso_ba(:,:),vector(:),vector_new(:)
real, allocatable :: diagonal(:)
real :: start_time, end_time
real :: total_time
call cpu_time(start_time)

  !********* INPUT **********
relativistic=.true.
  
norb=6

n_alpha=3
n_beta=3
  
n_RAS_spaces_occ=1
n_RAS_spaces_virt=1

allocate(RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt))
allocate(excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt))
! number of orbitals cannot be equal 0

RAS_space_occ(1)=2
!RAS_space_occ(2)=0
active_space=2
RAS_space_virt(1)=2
!RAS_space_virt(2)=0
excit_array(:)=1
  
allocate(hopping_alpha(norb,norb))
allocate(hopping_beta(norb,norb))
allocate(interaction_mix(norb,norb,norb,norb))
allocate(interaction_alpha(norb,norb,norb,norb))
allocate(interaction_beta(norb,norb,norb,norb))
hopping_alpha(:,:) = 0
hopping_beta(:,:) = 0
interaction_alpha(:,:,:,:) = 0
interaction_beta(:,:,:,:) = 0
interaction_mix(:,:,:,:) = 0
if (relativistic .eqv. .true.) then
    allocate(hso_ab(norb,norb))
    allocate(hso_ba(norb,norb))
    hso_ab(:,:) = 0
    hso_ba(:,:) = 0
end if 
!***********************
!here we check if input is ok 

call check_orbital_space_declarations(norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)
call check_hamiltonian(relativistic,norb,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix,hso_ab,hso_ba)

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

! here we calculate sizea,sizeb,size_tot

call calculate_space_size(relativistic,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,size_tot)

write(*,*) "Space size alpha: ", sizea(1,2), sizea(2,2), sizea(3,2)
write(*,*) "Space size beta: ", sizeb(1,2), sizeb(2,2), sizeb(3,2)

write(*,*) "Total size:"
write(*,*) size_tot(1,:)
write(*,*) size_tot(2,:)
write(*,*) size_tot(3,:)
  
  !size_tot(1:3,1:2) - an array of pointers to the blocks of the vector (or Hamiltonian)
  
  !order of blocks: (n_alpha, n_beta); (n_alpha+1, n_beta-1); (n_alpha-1, n_beta+1)
  
  !size_tot(1,:) - for sector (n_alpha, n_beta)
  !size_tot(2,:) - for sector (n_alpha+1, n_beta-1)
  !size_tot(3,:) - for sector (n_alpha-1, n_beta+1) 

  !size_tot(:,1) - the beginning of such sector
  !size_tot(:,2) - the size (number of elements) of/in such sector

!Here we check if any crucial blocks are empty

if (size_tot(1,2) .eq. 0) then
    write(*,*) "MAIN BLOCK IS EMPTY"
    write(*,*) "QUITTING THE PROGRAM"
    !SHOULD WE INLCUDE VACUUM?
    stop
else if ((relativistic .eqv. .true.) .and. (size_tot(2,2)+size_tot(3,2) .eq. 0)) then
    write(*,*) "RELATIVISTIC BLOCKS ARE EMPTY"
    write(*,*) "QUITTING THE PROGRAM"
    stop
end if

! here we generate strings in first block (alpha and beta separately)

if (n_alpha .gt. 0)  then
  allocate(str_a(sizea(1,2),n_alpha))
  call fill_spin_strings(n_alpha,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(1,1),NRDa(1,2),sizea(1,2),str_a)
else if (n_alpha .eq. 0)  then
  allocate(str_a(sizea(1,2),1))
  str_a(sizea(1,2),:) = 0 !vacuum  
end if

if (n_beta .gt. 0) then
  allocate(str_b(sizeb(1,2),n_beta))
  call fill_spin_strings(n_beta,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(1,1),NRDb(1,2),sizeb(1,2),str_b)
else if  (n_beta .eq. 0)  then
  allocate(str_b(sizeb(1,2),1))
  str_b(sizeb(1,2),:) = 0 !vacuum  
end if

!here we generate creation-annihilation matrices for first block

allocate(alpha_annihilation_creation_matrix(norb,norb,sizea(1,2),2))  
call fill_annihilation_creation_matrix(n_alpha,sizea(1,2),str_a,norb,alpha_annihilation_creation_matrix)

allocate(beta_annihilation_creation_matrix(norb,norb,sizeb(1,2),2))
call fill_annihilation_creation_matrix(n_beta,sizeb(1,2),str_b,norb,beta_annihilation_creation_matrix)

!here we generate single spin part of hamiltonian for first block

allocate(alpha_hamiltonian(sizea(1,2),sizea(1,2)))
allocate(beta_hamiltonian(sizeb(1,2),sizeb(1,2)))

call fill_nr_single_spin_hamiltonian(str_a,sizea(1,2),n_alpha,norb,hopping_alpha,interaction_alpha,alpha_annihilation_creation_matrix,alpha_hamiltonian)
call fill_nr_single_spin_hamiltonian(str_b,sizeb(1,2),n_beta,norb,hopping_beta,interaction_beta,beta_annihilation_creation_matrix,beta_hamiltonian)

! here we prepare for matrix-vector multiplication

allocate(vector(size_tot(1,2)+size_tot(2,2)+size_tot(3,2)))
allocate(vector_new(size_tot(1,2)+size_tot(2,2)+size_tot(3,2)))
allocate(diagonal(size_tot(1,2)+size_tot(2,2)+size_tot(3,2)))

call generate_diagonal_elements(alpha_hamiltonian,str_a,sizea(1,2),n_alpha,beta_hamiltonian,str_b,sizeb(1,2),n_beta,interaction_mix,norb,diagonal(1:size_tot(1,2)))

!here we repeat for other blocks

if (relativistic .eqv. .true.) then

!here we generate strings
    
    if (n_alpha .gt. 1)  then
        allocate(str_a_m1(sizea(3,2),n_alpha-1))
        if (sizea(3,2) .ne. 0) then
            call fill_spin_strings(n_alpha-1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(3,1),NRDa(3,2),sizea(3,2),str_a_m1)
        end if
    else if (n_alpha .eq. 1) then
        allocate(str_a_m1(sizea(3,2),1))
        str_a_m1(:,:) = 0 !vacuum
    else
        allocate(str_a_m1(sizea(3,2),0))
    end if
      
    if (n_beta .gt. 1)  then
        allocate(str_b_m1(sizeb(3,2),n_beta-1))
        if (sizeb(3,2) .ne. 0) then
            call fill_spin_strings(n_beta-1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(3,1),NRDb(3,2),sizeb(3,2),str_b_m1)
        end if
    else if (n_beta .eq. 1) then
        allocate(str_b_m1(sizeb(3,2),1))
        str_b_m1(:,:) = 0 !vacuum
    else 
        allocate(str_b_m1(sizeb(3,2),0))
    end if
    
    allocate(str_a_p1(sizea(2,2),n_alpha+1))
    if (sizea(2,2) .ne. 0) then
     call fill_spin_strings(n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),str_a_p1)
    end if
      
    allocate(str_b_p1(sizeb(2,2),n_beta+1)) 
    if (sizeb(2,2) .ne. 0) then
        call fill_spin_strings(n_beta+1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(2,1),NRDb(2,2),sizeb(2,2),str_b_p1)
    end if

!here we generate annihilation, creation-annihilation and single spin hamiltonians
      
    allocate(alpha_annihilation_matrix(norb,sizea(1,2),2))
    allocate(alpha_annihilation_creation_matrix_m1(norb,norb,sizea(3,2),2))
    allocate(alpha_hamiltonian_m1(sizea(3,2),sizea(3,2)))

    if (sizea(3,2) .ne. 0) then      
        call fill_annihilation_results(sizea(3,2),sizea(1,2),n_alpha,str_a_m1,str_a,norb,alpha_annihilation_matrix)
        call fill_annihilation_creation_matrix(n_alpha-1,sizea(3,2),str_a_m1,norb,alpha_annihilation_creation_matrix_m1)
        call fill_nr_single_spin_hamiltonian(str_a_m1,sizea(3,2),n_alpha-1,norb,hopping_alpha,interaction_alpha,alpha_annihilation_creation_matrix_m1,alpha_hamiltonian_m1)
    end if
    
    allocate(beta_annihilation_matrix(norb,sizeb(1,2),2))
    allocate(beta_annihilation_creation_matrix_m1(norb,norb,sizeb(3,2),2))
    allocate(beta_hamiltonian_m1(sizeb(3,2),sizeb(3,2)))
    
    if (sizeb(3,2) .ne. 0) then      
        call fill_annihilation_results(sizeb(3,2),sizeb(1,2),n_beta,str_b_m1,str_b,norb,beta_annihilation_matrix)
        call fill_annihilation_creation_matrix(n_beta-1,sizeb(3,2),str_b_m1,norb,beta_annihilation_creation_matrix_m1)
        call fill_nr_single_spin_hamiltonian(str_b_m1,sizeb(3,2),n_beta-1,norb,hopping_beta,interaction_beta,beta_annihilation_creation_matrix_m1,beta_hamiltonian_m1)
    end if
 

    allocate(alpha_annihilation_matrix_p1(norb,sizea(2,2),2))
    allocate(alpha_annihilation_creation_matrix_p1(norb,norb,sizea(2,2),2))
    allocate(alpha_hamiltonian_p1(sizea(2,2),sizea(2,2)))
    
    if (sizea(2,2) .ne. 0) then
        call fill_spin_strings(n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),str_a_p1)
        call fill_annihilation_results(sizea(1,2),sizea(2,2),n_alpha+1,str_a,str_a_p1,norb,alpha_annihilation_matrix_p1)
        call fill_annihilation_creation_matrix(n_alpha+1,sizea(2,2),str_a_p1,norb,alpha_annihilation_creation_matrix_p1)
        call fill_nr_single_spin_hamiltonian(str_a_p1,sizea(2,2),n_alpha+1,norb,hopping_alpha,interaction_alpha,alpha_annihilation_creation_matrix_p1,alpha_hamiltonian_p1)
    end if

    allocate(beta_annihilation_matrix_p1(norb,sizeb(2,2),2))  
    allocate(beta_annihilation_creation_matrix_p1(norb,norb,sizeb(2,2),2))
    allocate(beta_hamiltonian_p1(sizeb(2,2),sizeb(2,2)))   

    if (sizeb(2,2) .ne. 0) then
        call fill_spin_strings(n_beta+1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(2,1),NRDb(2,2),sizeb(2,2),str_b_p1)
        call fill_annihilation_results(sizeb(1,2),sizeb(2,2),n_beta+1,str_b,str_b_p1,norb,beta_annihilation_matrix_p1)
        call fill_annihilation_creation_matrix(n_beta+1,sizeb(2,2),str_b_p1,norb,beta_annihilation_creation_matrix_p1)
        call fill_nr_single_spin_hamiltonian(str_b_p1,sizeb(2,2),n_beta+1,norb,hopping_beta,interaction_beta,beta_annihilation_creation_matrix_p1,beta_hamiltonian_p1)
    end if

! here we generate hamiltonian diagonal elements
    if (size_tot(2,2) .ne. 0) then
        call generate_diagonal_elements(alpha_hamiltonian_p1,str_a_p1,sizea(2,2),n_alpha+1,beta_hamiltonian_m1,str_b_m1,sizeb(3,2),n_beta-1,interaction_mix,norb,diagonal(size_tot(2,1):size_tot(2,2)+size_tot(1,2)))
    end if

    if (size_tot(3,2) .ne. 0) then 
        call generate_diagonal_elements(alpha_hamiltonian_m1,str_a_m1,sizea(3,2),n_alpha-1,beta_hamiltonian_p1,str_b_p1,sizeb(2,2),n_beta+1,interaction_mix,norb,diagonal(size_tot(3,1):size_tot(3,2)+size_tot(2,2)+size_tot(1,2)))
    end if

end if

vector = [0,0,0,0,1]
vector_new(:) = 0

if (relativistic .eqv. .false.) then
    call nr_matrix_vector_product(alpha_hamiltonian,str_a,sizea(1,2),n_alpha,beta_hamiltonian,str_b,sizeb(1,2),n_beta,interaction_mix,norb,alpha_annihilation_creation_matrix,beta_annihilation_creation_matrix,vector,vector_new)
else
    call rel_matrix_vector_product(alpha_hamiltonian,alpha_hamiltonian_p1,alpha_hamiltonian_m1,str_a,sizea(1,2),str_a_p1,sizea(2,2),str_a_m1,sizea(3,2),n_alpha,beta_hamiltonian,beta_hamiltonian_p1,beta_hamiltonian_m1,str_b,sizeb(1,2),str_b_p1,sizeb(2,2),str_b_m1,sizeb(3,2),n_beta,interaction_mix,hso_ab,hso_ba,norb,alpha_annihilation_creation_matrix,alpha_annihilation_creation_matrix_p1,alpha_annihilation_creation_matrix_m1,beta_annihilation_creation_matrix,beta_annihilation_creation_matrix_p1,beta_annihilation_creation_matrix_m1,alpha_annihilation_matrix,alpha_annihilation_matrix_p1,beta_annihilation_matrix,beta_annihilation_matrix_p1,size_tot,vector,vector_new)
end if
write(*,*) vector_new
write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time, size_tot(1,2),size_tot(2,2),size_tot(3,2),size_tot(1,2)+size_tot(2,2)+size_tot(3,2)
call flush(6)
end program relativistic_ed