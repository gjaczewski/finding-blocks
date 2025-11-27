program relatvistic_ed
  implicit none
  integer:: n_alpha,n_beta,norb
  integer:: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, dimension(:), allocatable:: RAS_space_occ,RAS_space_virt,excit_array
  integer:: NRD_spin_alpha,NRD_spin_beta,NRD_spin_tmp
  integer, dimension(:,:), allocatable:: RAS_el_array_alpha,RAS_el_array_beta
  integer:: orb,n_spaces
  logical:: relativistic
  double precision, dimension(:), allocatable:: scratch
  integer:: max_scratch
  integer, dimension(:), allocatable:: scratch_int
  integer:: max_scratch_int
  
  integer:: i,j,k,l
  integer:: count_unique


  !norb - total number of orbitals
  ! the principles for the occupied and virtual RAS spaces are the same
  !n_RAS_spaces_occ - the number of RAS spaces which are present for occupied orbitals (eg. 3 spaces with core, singly, and doubly occupied orbitals out of total number of occupied orbitals)
  !RAS_space_occ(number_RAS_spaces_occ) - dimension/number of orbitals in each of the occupied RAS space (eg. dimension_RAS_space_occ(1)=number of orbitals in the first RAS space
  !active_space - number of orbitals in the active space 
  !excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt) - an array describing the excitation level possible in each of the RAS spaces
  !relativistic - logical that is switching between relativistic and nonrelativistic calculations (relativistic = .true. is a relativistic calculation)

  max_scratch=10000000
  allocate(scratch(max_scratch))
  max_scratch_int=1000000
  allocate(scratch_int(max_scratch_int))
  !in a real program we will get these from the input
  
  relativistic=.true.
  norb=14

  n_alpha=6
  n_beta=6
  
  n_RAS_spaces_occ=2
  n_RAS_spaces_virt=2



  allocate(RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt))


  !in a real program we will get these from the input
  
  RAS_space_occ(1)=3
  RAS_space_occ(2)=2
  active_space=4
  RAS_space_virt(1)=3
  RAS_space_virt(2)=2

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

  n_spaces=n_RAS_spaces_occ+1+n_RAS_spaces_virt !total number of RAS spaces+active space
  allocate(RAS_el_array_alpha(NRD_spin_tmp,n_spaces))
  allocate(RAS_el_array_beta(NRD_spin_tmp,n_spaces))

  NRD_spin_alpha=0
  NRD_spin_beta=0
  
  call fill_RAS_electron_array(relativistic,norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,n_alpha,n_beta,NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta)

  k=1
  do i=1,NRD_spin_alpha
     do j=1,1 !n_RAS_spaces_virt+n_RAS_spaces_virt+1
        write(*,*)"RAS_el_array_alpha",i,j,RAS_el_array_alpha(i,j)
        scratch_int(k)=RAS_el_array_alpha(i,j)
!        write(*,*)"scratch(k)",k,scratch(k)
        k=k+1
     end do
  end do
  STOP
  !now try to find unique number of electrons in each of the spaces

!  call unique_elements(k-1,scratch_int(1:k-1),scratch_int(k),count_unique)
!  do i=1,count_unique
!     write(*,*)"unique el num",i,scratch_int(k+i)
!  end do
  

end program relatvistic_ed
