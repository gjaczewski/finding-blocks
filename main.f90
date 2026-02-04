program new_rel

  implicit none
  integer:: n_alpha,n_beta,norb
  integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space, n_spaces
  integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
  integer, dimension(:,:), allocatable:: RAS_el_array_alpha,RAS_el_array_beta
  logical:: relativistic
  integer ::   NRD_spin_tmp, NRD_spin_alpha, NRD_spin_beta
  integer :: i, j,k,l,p,q, mu ,nu, temp, n_distributions_alpha,n_distributions_beta, n_distributions_alpha_p1,n_distributions_beta_p1, n_distributions_alpha_m1,n_distributions_beta_m1
  integer :: max,n_combinations
  integer, allocatable :: all_combinations(:,:)
  integer:: NRDa(3,2),NRDb(3,2),sizea(3,2),sizeb(3,2),size_tot(3,2)
  integer, dimension(:,:), allocatable :: str_a, str_b, str_a_p1, str_a_m1, str_b_p1, str_b_m1
  integer, allocatable :: alpha_annihilation_matrix(:,:,:), alpha_annihilation_matrix_p1(:,:,:), beta_annihilation_matrix(:,:,:), beta_annihilation_matrix_p1(:,:,:)
  integer, allocatable :: alpha_annihilation_creation_matrix(:,:,:,:), beta_annihilation_creation_matrix(:,:,:,:),alpha_annihilation_creation_matrix_p1(:,:,:,:), beta_annihilation_creation_matrix_p1(:,:,:,:),alpha_annihilation_creation_matrix_m1(:,:,:,:), beta_annihilation_creation_matrix_m1(:,:,:,:)
  integer :: verbose
  integer :: n_rows_alpha, n_rows_beta

  complex, allocatable :: hopping(:,:), interaction(:,:,:,:),hso(:,:),vector(:),vector_new(:)
  real :: start_time, end_time
  real :: total_time
  call cpu_time(start_time)
  !********* INPUT **********
  verbose = 3
  relativistic=.true.
  norb=10

  n_alpha=5
  n_beta=5
  
  n_RAS_spaces_occ=2
  n_RAS_spaces_virt=2

  allocate(RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt))
! number of orbitals cannot be equal 0

  RAS_space_occ(1)=2
  RAS_space_occ(2)=2
  active_space=2
  RAS_space_virt(1)=2
  RAS_space_virt(2)=2

  
  !***********************
  
  allocate(excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt))
  excit_array(:)=2

  call check_orbital_space_declarations(norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)

  !estimate maximal theoretical limit = product of el in occupied ras spaces*number of el in active_space*product of el in virtutal ras spaces

  NRD_spin_tmp=1
  do i=1,n_RAS_spaces_occ 
     NRD_spin_tmp=NRD_spin_tmp*RAS_space_occ(i)
  end do

  NRD_spin_tmp=NRD_spin_tmp*active_space
  do i=1, n_RAS_spaces_virt
     NRD_spin_tmp=NRD_spin_tmp*RAS_space_virt(i)
  end do
  !write(*,*)"NRD_spin_tmp",NRD_spin_tmp
  call flush(6)

  if (relativistic .eqv. .true. ) then 
     NRD_spin_tmp=3*NRD_spin_tmp !since we are considering n_spin, nspin+1, and n_spin-1 for a relativistic case
  !     write(*,*)"NRD_spin_tmp",NRD_spin_tmp
  call flush(6)
  end if


! calculate all possible combinations
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


! here we calculate number of possible distributions NRD_spin_alpha, NRD_spin_beta and number of distribution for fixed spin
call count_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta)


allocate(RAS_el_array_alpha(NRD_spin_alpha,n_spaces))
allocate(RAS_el_array_beta(NRD_spin_beta,n_spaces))

write(*,*) "NRD_spin_alpha: ",NRD_spin_alpha
write(*,*) "NRD_spin_beta: ",NRD_spin_beta
call flush(6)


! here we fill RAS_el_array_spin tables and calculate NRDa, NRDb
call find_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,NRDa,NRDb)

write(*,*) "NRDa:"
write(*,*) NRDa(1,:)
write(*,*) NRDa(2,:)
write(*,*) NRDa(3,:)

write(*,*) "NRDb: "
write(*,*) NRDa(1,:)
write(*,*) NRDa(2,:)
write(*,*) NRDa(3,:)

! here we calculate sizea(3,2),sizeb(3,2),size_tot(3,2) 
  call calculate_space_size(NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,size_tot,verbose)
write(*,*) "Space size alpha: ", sizea(1,2), sizea(2,2), sizea(3,2)
write(*,*) "Space size beta: ", sizeb(1,2), sizeb(2,2), sizeb(3,2)
write(*,*) "Total size:"
write(*,*) size_tot(1,:)
write(*,*) size_tot(2,:)
write(*,*) size_tot(3,:)
call flush(6)

  !size_tot(1:3,1:2) - an array of pointers to the blocks of the vector (or Hamiltonian)
  !order of blocks: (n_alpha, n_beta); (n_alpha+1, n_beta-1); (n_alpha-1, n_beta+1)
  !size_tot(1,:) - for sector (n_alpha, n_beta)
  !size_tot(2,:) - for sector (n_alpha+1, n_beta-1)
  !size_tot(3,:) - for sector (n_alpha-1, n_beta+1) 

  !size_tot(:,1) - the beginning of such sector
  !size_tot(:,2) - the size (number of elements) of/in such sector

  ! here we generate strings (alpha and beta separately)
  allocate(str_a(sizea(1,2),n_alpha))
  call fill_spin_strings(n_alpha,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(1,1),NRDa(1,2),sizea(1,2),str_a,verbose)


  allocate(str_b(sizeb(1,2),n_beta))
  call fill_spin_strings(n_beta,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(1,1),NRDb(1,2),sizeb(1,2),str_b,verbose)
  
  if (n_alpha .gt. 1) then
  allocate(str_a_m1(sizea(3,2),n_alpha-1))
  call fill_spin_strings(n_alpha-1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(3,1),NRDa(3,2),sizea(3,2),str_a_m1,verbose)
  else
  allocate(str_a_m1(sizea(3,2),1))
  str_a_m1(sizea(3,2),:) = 0 !vacuum
  end if

  if (n_beta .gt. 1) then
  allocate(str_b_m1(sizeb(3,2),n_beta-1))
  call fill_spin_strings(n_beta-1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(3,1),NRDb(3,2),sizeb(3,2),str_b_m1,verbose)
  else
  allocate(str_b_m1(sizeb(3,2),1))
  str_b_m1(sizeb(3,2),:) = 0 !vacuum
  end if
  
  if (relativistic .eqv. .true.) then
     allocate(str_a_p1(sizea(2,2),n_alpha+1))
     call fill_spin_strings(n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),str_a_p1,verbose)

     allocate(str_b_p1(sizeb(2,2),n_beta+1))
     call fill_spin_strings(n_beta,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(1,1),NRDb(1,2),sizeb(1,2),str_b_p1,verbose)
  end if

!here we generate how annihilation operator acts on states
allocate(alpha_annihilation_matrix(norb,sizea(1,2),2))
allocate(beta_annihilation_matrix(norb,sizeb(1,2),2))

call fill_annihilation_results(sizea(3,2),sizea(1,2),n_alpha,str_a_m1,str_a,norb,alpha_annihilation_matrix)
call fill_annihilation_results(sizeb(3,2),sizeb(1,2),n_beta,str_b_m1,str_b,norb,beta_annihilation_matrix)

if (relativistic .eqv. .true.) then
allocate(alpha_annihilation_matrix_p1(norb,sizea(2,2),2))
allocate(beta_annihilation_matrix_p1(norb,sizeb(2,2),2))
call fill_annihilation_results(sizea(1,2),sizea(2,2),n_alpha+1,str_a,str_a_p1,norb,alpha_annihilation_matrix_p1)
call fill_annihilation_results(sizeb(1,2),sizeb(2,2),n_beta+1,str_b,str_b_p1,norb,beta_annihilation_matrix_p1)
end if

!here we generate how pair of operators acts on states
allocate(alpha_annihilation_creation_matrix(norb,norb,sizea(1,2),2))
allocate(beta_annihilation_creation_matrix(norb,norb,sizeb(1,2),2))

call fill_annihilation_creation_matrix(n_alpha,sizea(1,2),str_a,norb,alpha_annihilation_creation_matrix)
call fill_annihilation_creation_matrix(n_beta,sizeb(1,2),str_b,norb,beta_annihilation_creation_matrix)

if (relativistic .eqv. .true.) then
allocate(alpha_annihilation_creation_matrix_p1(norb,norb,sizea(2,2),2))
allocate(beta_annihilation_creation_matrix_p1(norb,norb,sizeb(2,2),2))

allocate(alpha_annihilation_creation_matrix_m1(norb,norb,sizea(3,2),2))
allocate(beta_annihilation_creation_matrix_m1(norb,norb,sizeb(3,2),2))

call fill_annihilation_creation_matrix(n_alpha+1,sizea(2,2),str_a_p1,norb,alpha_annihilation_creation_matrix_p1)
call fill_annihilation_creation_matrix(n_beta+1,sizeb(2,2),str_b_p1,norb,beta_annihilation_creation_matrix_p1)

call fill_annihilation_creation_matrix(n_alpha-1,sizea(3,2),str_a_m1,norb,alpha_annihilation_creation_matrix_m1)
call fill_annihilation_creation_matrix(n_beta-1,sizeb(3,2),str_b_m1,norb,beta_annihilation_creation_matrix_m1)
end if



allocate(hopping(norb,norb))
allocate(interaction(norb,norb,norb,norb))
hopping(:,:) = 1
interaction(:,:,:,:) = 1
if (relativistic .eqv. .true.) then
allocate(hso(norb,norb))
hso(:,:) = 1
end if 

allocate(vector(sizea(1,2)*sizeb(1,2)))
allocate(vector_new(sizea(1,2)*sizeb(1,2)))
vector(:) = 1
call nr_matrix_vector_product(hopping,interaction,norb,str_a,sizea(1,2),n_alpha,str_b,sizeb(1,2),n_beta,alpha_annihilation_creation_matrix,beta_annihilation_creation_matrix,vector,vector_new)


write(*,*) "KONIEC"
call cpu_time(end_time)
total_time = end_time - start_time
write(*,*) "Execution time: ", total_time, size_tot(1,1),size_tot(2,1),size_tot(3,1)
call flush(6)
end program new_rel




















