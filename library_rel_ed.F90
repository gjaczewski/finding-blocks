#include <slepc/finclude/slepceps.h>
module library_rel_ed
use slepceps
type :: parameters
    sequence
    integer :: n_strings_alpha, n_strings_beta, n_orb
    integer :: n_alpha, n_beta
    integer :: n_strings_alpha_p1, n_strings_beta_p1
    integer :: n_strings_alpha_m1, n_strings_beta_m1
    integer, allocatable :: alpha_annihilation_creation_matrix(:,:,:,:)
    integer, allocatable :: beta_annihilation_creation_matrix(:,:,:,:)
    integer, allocatable :: str_a(:,:), str_b(:,:)
    integer, allocatable :: alpha_annihilation_creation_matrix_p1(:,:,:,:)
    integer, allocatable :: beta_annihilation_creation_matrix_p1(:,:,:,:)
    integer, allocatable :: str_a_p1(:,:), str_b_p1(:,:)
    integer, allocatable :: alpha_annihilation_creation_matrix_m1(:,:,:,:)
    integer, allocatable :: beta_annihilation_creation_matrix_m1(:,:,:,:)
    integer, allocatable :: str_a_m1(:,:), str_b_m1(:,:)
    integer, allocatable :: alpha_annihilation_matrix(:,:,:)
    integer, allocatable :: beta_annihilation_matrix(:,:,:)
    integer, allocatable :: alpha_annihilation_matrix_p1(:,:,:)
    integer, allocatable :: beta_annihilation_matrix_p1(:,:,:)
    integer :: size_tot(3,2)
    complex(8), allocatable :: alpha_hamiltonian(:,:), beta_hamiltonian(:,:)
    complex(8), allocatable :: interaction_mix(:,:,:,:), hso_ab(:,:), hso_ba(:,:)
    complex(8), allocatable :: alpha_hamiltonian_p1(:,:), beta_hamiltonian_p1(:,:)
    complex(8), allocatable :: alpha_hamiltonian_m1(:,:), beta_hamiltonian_m1(:,:)
    complex(8), allocatable:: diag(:)
    real(8) :: nuclear_energy
    logical :: relativistic
end type parameters
  
  type(parameters), pointer :: ptr_dane
  complex(8), pointer       :: ptr_przekatna(:)
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
subroutine check_orbital_space_declarations(n_orb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)

  implicit none
  integer, intent(in) :: n_alpha, n_beta
  integer, intent(in) :: n_orb,n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer:: orb,i


  !n_orb - total number of orbitals
  ! the principles for the occupied and virtual RAS spaces are the same
  !number_RAS_spaces_occ - the number of RAS spaces which are present for occupied orbitals (eg. 3 spaces with core, singly, and doubly occupied orbitals out of total number of occupied orbitals)
  !dimension_RAS_space_occ(number_RAS_spaces_occ) - dimension/number of orbitals in each of the occupied RAS space (eg. dimension_RAS_space_occ(1)=number of orbitals in the first RAS space
  !dimension_active_space - number of orbitals in the active space 
  
  !check if number of all the orbitals is equal to n_orb

  orb=0
  do i=1, n_RAS_spaces_occ
     orb=orb+RAS_space_occ(i)
  end do

  orb=orb+active_space
  
  do i=1, n_RAS_spaces_virt
     orb=orb+RAS_space_virt(i)
  end do

  if (abs(n_RAS_spaces_occ-n_RAS_spaces_virt)>0) then
     write(*,*)"NUMBER of RAS SPACE OCC MUST BE EQUAL to the NUMBER of RAS SPACE VIRT",n_RAS_spaces_occ,n_RAS_spaces_virt
     call flush(6)
     STOP
  end if

  if (orb .ne. n_orb) then
     write(*,*)"SOMETHING WENT WRONG,number of orbitals declared in RAS spaces or active space is wrong (number of orbitals, number or orbitals in RAS+ACTIVE)", n_orb, orb
     call flush(6)
     do i=1, n_RAS_spaces_occ
        write(*,*)"RAS SPACE OCC",i,"number of orbitals", RAS_space_occ(i)
        call flush(6)
     end do

     write(*,*)"Active_space",active_space
     call flush(6)
     
     

     do i=1, n_RAS_spaces_virt
        write(*,*)"RAS SPACE VIRT",i,"number of orbitals", RAS_space_virt(i)
        call flush(6)
     end do

     write(*,*)"QUITING THE PROGRAM"
     call flush(6)
     STOP
     
  end if

  !check if number of all the n_alpha is smaller or equal to orb
  
  if (orb<n_alpha) then
     write(*,*)"SOMETHING WENT WRONG,number of declared alpha electrons", n_alpha, "is greater than total number of orbitals", orb
     call flush(6)
     write(*,*)"QUITING THE PROGRAM"
     call flush(6)
     STOP
  end if


  !check if number of all the n_beta is smaller or equal to orb
  
    if (orb<n_beta) then
     write(*,*)"SOMETHING WENT WRONG,number of declared alpha electrons", n_beta, "is greater than total number of orbitals", orb
     call flush(6)
     write(*,*)"QUITING THE PROGRAM"
     call flush(6)
     STOP
  end if
end subroutine check_orbital_space_declarations

subroutine check_hamiltonian(relativistic,n_orb,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix,hso_ab,hso_ba)
logical :: relativistic
integer, intent(in) :: n_orb
complex(8), intent(inout) :: hopping_alpha(n_orb,n_orb),hopping_beta(n_orb,n_orb), interaction_alpha(n_orb,n_orb,n_orb,n_orb), interaction_beta(n_orb,n_orb,n_orb,n_orb),interaction_mix(n_orb,n_orb,n_orb,n_orb),hso_ab(n_orb,n_orb),hso_ba(n_orb,n_orb)
integer :: i,j,k,l
do i=1,n_orb
   do j=1,n_orb
      if (abs(hopping_alpha(i,j)-conjg(hopping_alpha(j,i))) .gt. 1E-10) then
         write(*,*) "Hopping alpha is not hermitian"
         write(*,*) "QUITTING THE PROGRAM"
         STOP
      else
         hopping_alpha(i,j) = conjg(hopping_alpha(j,i))
      end if
      if (abs(hopping_beta(i,j)-conjg(hopping_beta(j,i))) .gt. 1E-10) then
         write(*,*) "Hopping beta is not hermitian"
         write(*,*) "QUITTING THE PROGRAM"
         STOP
      else
         hopping_beta(i,j) = conjg(hopping_beta(j,i))
      end if
      if (relativistic .eqv. .true.) then
         if (abs(hso_ab(i,j)-conjg(hso_ba(j,i))) .gt. 1E-10) then
            write(*,*) "Spin orbit hamiltonian is not hermitian"
            write(*,*) "QUITTING THE PROGRAM"
            STOP
         else
            hso_ab(i,j) = conjg(hso_ba(j,i))
         end if
      end if
      do k=1,n_orb
         do l=1,n_orb
            if ((abs(interaction_mix(i,j,k,l)-conjg(interaction_mix(k,l,i,j))) .gt. 1E-10 ) .or. (abs(interaction_mix(i,j,k,l)-interaction_mix(j,i,l,k)) .gt. 1E-10 )) then
               write(*,*) i,j,k,l
               write(*,*) abs(interaction_mix(i,j,k,l)-conjg(interaction_mix(k,l,i,j))), abs(interaction_mix(i,j,k,l)-interaction_mix(j,i,l,k)) 
               write(*,*) "Mixed interaction is not hermitian or not symmetric"
               write(*,*) "QUITTING THE PROGRAM"
               STOP
            else
               interaction_mix(i,j,k,l) = conjg(interaction_mix(k,l,i,j))
               interaction_mix(i,j,k,l) = interaction_mix(j,i,l,k)
            end if
            if ((abs(interaction_alpha(i,j,k,l)-conjg(interaction_alpha(k,l,i,j))) .gt. 1E-10) .or. (abs(interaction_alpha(i,j,k,l)-interaction_alpha(j,i,l,k)) .gt. 1E-10)) then
               write(*,*) "Alpha interaction is not hermitian or not symmetric"
               write(*,*) "QUITTING THE PROGRAM"
               STOP
            else
               interaction_alpha(i,j,k,l) = conjg(interaction_alpha(k,l,i,j))
               interaction_alpha(i,j,k,l) = interaction_alpha(j,i,l,k)
            end if
            if ((abs(interaction_beta(i,j,k,l)-conjg(interaction_beta(k,l,i,j))) .gt. 1E-10) .or. (abs(interaction_beta(i,j,k,l)-interaction_beta(j,i,l,k)) .gt. 1E-10)) then
               write(*,*) "Beta interaction is not hermitian or not symmetric"
               write(*,*) "QUITTING THE PROGRAM"
               STOP
            else
               interaction_beta(i,j,k,l) = conjg(interaction_beta(k,l,i,j))
               interaction_beta(i,j,k,l) = interaction_beta(j,i,l,k)
            end if
         end do
      end do
   end do
end do
end subroutine check_hamiltonian

subroutine count_distributions(n_s, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions)
implicit none
integer, intent(in) :: n_s, n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
integer ::i,j
integer, intent(in) :: n_combinations,n_spaces
integer, intent(in) :: all_combinations(n_combinations,n_spaces)
integer :: a,v,r, counter, sum, checker
integer, intent(inout) :: n_distributions
a = n_RAS_spaces_occ + 1
counter = 0
do i=1,n_combinations
   sum = 0
 
   do j=1,n_spaces
      sum = sum + all_combinations(i,j)
   end do
   if (sum .eq. n_s) then
      checker = n_spaces
      do r=1,n_RAS_spaces_occ
       
         if ((all_combinations(i,r) .lt. RAS_space_occ(r) - excit_array(r)) .or. (all_combinations(i,r) .gt. RAS_space_occ(r) )) then
            checker = checker - 1
         end if
      end do
 
 
      if  (all_combinations(i,a) .gt. active_space ) then
                     checker = checker - 1
      end if
      do v=1,n_RAS_spaces_virt
        if ((all_combinations(i,a+v) .gt.  excit_array(n_RAS_spaces_occ+v)) ) then
         
   
                     checker = checker - 1
        end if

      end do
   if (checker == n_spaces) then
      counter  = counter + 1
   end if


   end if
end do
n_distributions = counter
end subroutine count_distributions


subroutine find_distributions(n_s, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions,distributions)
implicit none
integer, intent(in) :: n_s, n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
integer :: n_combinations,n_spaces,i,j
integer,intent(inout) :: distributions(n_distributions,n_spaces)
integer, intent(in) ::  all_combinations(n_combinations,n_spaces)
integer :: a,v,r, counter, sum, checker
integer, intent(in) :: n_distributions
a = n_RAS_spaces_occ + 1
counter = 0
do i=1,n_combinations
   sum = 0
   do j=1,n_spaces
      sum = sum + all_combinations(i,j)
   end do
   if (sum .eq. n_s) then
      checker = n_spaces
      do r=1,n_RAS_spaces_occ
         
         if ((all_combinations(i,r) .lt. RAS_space_occ(r) - excit_array(r)) .or. (all_combinations(i,r) .gt. RAS_space_occ(r) )) then
            checker = checker - 1
         end if
      end do

 
      if  (all_combinations(i,a) .gt. active_space ) then
         checker = checker - 1
      end if
      do v=1,n_RAS_spaces_virt
        if (all_combinations(i,a+v) .gt. excit_array(n_RAS_spaces_occ+v))  then
         checker = checker - 1
        end if

      end do
   if (checker == n_spaces) then
      counter  = counter + 1

      distributions(counter,:) = all_combinations(i,:)
   end if

   end if
end do
end subroutine find_distributions


subroutine count_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta)
implicit none
logical, intent(in) :: relativistic
integer, intent(in) :: n_alpha,n_beta
integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
integer, intent(in) :: n_combinations,n_spaces
integer, intent(in) :: all_combinations(n_combinations,n_spaces)
integer, intent(out) :: n_distributions_alpha,n_distributions_beta
integer, intent(out) :: n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1
integer, intent(out) :: NRD_spin_alpha,NRD_spin_beta
n_distributions_alpha_p1 = 0    
n_distributions_beta_p1 = 0
n_distributions_alpha_m1 = 0    
n_distributions_beta_m1 = 0
call count_distributions(n_alpha, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha)
call count_distributions(n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta)

  call count_distributions(n_alpha-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_m1)
  call count_distributions(n_beta-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_m1)

  call count_distributions(n_alpha+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_p1)
  call count_distributions(n_beta+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_p1)


NRD_spin_alpha = n_distributions_alpha + n_distributions_alpha_m1 + n_distributions_alpha_p1
NRD_spin_beta= n_distributions_beta + n_distributions_beta_m1 + n_distributions_beta_p1
end subroutine count_spin_distributions


subroutine find_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,NRDa,NRDb)
implicit none
logical, intent(in) :: relativistic
integer, intent(in) :: n_alpha,n_beta
integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
integer, intent(in) :: n_combinations,n_spaces
integer, intent(in) :: all_combinations(n_combinations,n_spaces)
integer, intent(in) :: n_distributions_alpha,n_distributions_beta
integer, intent(in) :: n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1
integer, intent(in) :: NRD_spin_alpha,NRD_spin_beta
integer, intent(out) :: RAS_el_array_alpha(NRD_spin_alpha,n_spaces), RAS_el_array_beta(NRD_spin_beta,n_spaces)
integer, intent(out):: NRDa(3,2),NRDb(3,2)
NRDa(1,1) = 1
NRDa(1,2) = n_distributions_alpha

NRDa(2,1) = n_distributions_alpha+1
NRDa(2,2) = n_distributions_alpha + n_distributions_alpha_p1

NRDa(3,1) = n_distributions_alpha + n_distributions_alpha_p1 + 1
NRDa(3,2) = n_distributions_alpha + n_distributions_alpha_p1 + n_distributions_alpha_m1

NRDb(1,1) = 1
NRDb(1,2) = n_distributions_beta

NRDb(2,1) = n_distributions_beta+1
NRDb(2,2) = n_distributions_beta + n_distributions_beta_p1

NRDb(3,1) = n_distributions_beta + n_distributions_beta_p1 + 1
NRDb(3,2) = n_distributions_beta + n_distributions_beta_p1 + n_distributions_beta_m1

call find_distributions(n_alpha, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,RAS_el_array_alpha(1:n_distributions_alpha,:))
call find_distributions(n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta,RAS_el_array_beta(1:n_distributions_beta,:))


call find_distributions(n_alpha+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_p1,RAS_el_array_alpha(n_distributions_alpha+1:n_distributions_alpha+n_distributions_alpha_p1,:))
call find_distributions(n_beta+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_p1,RAS_el_array_beta(n_distributions_beta+1:n_distributions_beta+n_distributions_beta_p1,:))
call find_distributions(n_alpha-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_m1,RAS_el_array_alpha(n_distributions_alpha+n_distributions_alpha_p1+1:n_distributions_alpha_m1+n_distributions_alpha+n_distributions_alpha_p1,:))
call find_distributions(n_beta-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_m1,RAS_el_array_beta(n_distributions_beta+n_distributions_beta_p1+1:n_distributions_beta_m1+n_distributions_beta+n_distributions_beta_p1,:))


end subroutine find_spin_distributions


subroutine calculate_space_size(relativistic,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,size_tot)
  implicit none
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in)::NRDa(3,2),NRDb(3,2)
  integer, intent(in):: NRD_spin_alpha,NRD_spin_beta
  integer, intent(in):: RAS_el_array_alpha(NRD_spin_alpha,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_el_array_beta(NRD_spin_beta,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(out)::sizea(3,2),sizeb(3,2),size_tot(3,2) 
  logical, intent(in) :: relativistic

  
  ! n_alpha n_beta block
  sizea(1,1)=1

  call calc_size_block(NRDa(1,1),NRDa(1,2),NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizea(1,2))

  sizeb(1,1)=1 
  call calc_size_block(NRDb(1,1),NRDb(1,2),NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizeb(1,2))

  size_tot(1,1)=1
  size_tot(1,2)=sizea(1,2)*sizeb(1,2)
  

  
  ! n_alpha+1 n_beta-1 block
  sizea(2,1)=1 
  call calc_size_block(NRDa(2,1),NRDa(2,2),NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizea(2,2))

  sizeb(3,1)=1 
  call calc_size_block(NRDb(3,1),NRDb(3,2),NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizeb(3,2))
  
  size_tot(2,1)=size_tot(1,2)+1
  size_tot(2,2)=sizea(2,2)*sizeb(3,2)


  ! n_alpha-1 n_beta+1 block

  sizea(3,1)=1 
  call calc_size_block(NRDa(3,1),NRDa(3,2),NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizea(3,2))
  sizeb(2,1)=1 
  call calc_size_block(NRDb(2,1),NRDb(2,2),NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizeb(2,2))
  
  size_tot(3,1)=size_tot(2,1)+sizea(2,2)*sizeb(3,2)  
  size_tot(3,2)=sizea(3,2)*sizeb(2,2)

  if (relativistic .eqv. .false.) then
   size_tot(2,2) = 0
   size_tot(3,2) = 0
  end if
end subroutine calculate_space_size


subroutine calc_size_block(range1,range2,NRD_spin,RAS_el_array_spin,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,vec_length)
  implicit none
  integer, intent(in):: range1,range2
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in):: NRD_spin
  integer, intent(in):: RAS_el_array_spin(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(out)::vec_length

  integer:: i,r,v,a
  integer:: siz
  !real(8) :: nCr_dp
  real(8) :: c
  vec_length=0
 
  do i=range1,range2
     a=n_RAS_spaces_occ+1 
     siz=1
     do r=1,n_RAS_spaces_occ
        call bin_coeff(RAS_space_occ(r),RAS_el_array_spin(i,r),c)
        siz=siz*int(c)
     end do
     do v=1,n_RAS_spaces_virt
        call bin_coeff(RAS_space_virt(v),RAS_el_array_spin(i,a+v),c)
        siz=siz*int(c)
     end do
     call bin_coeff(active_space,RAS_el_array_spin(i,a),c)

     siz=siz*int(c)
     vec_length=vec_length+siz
  end do

end subroutine calc_size_block

subroutine bin_coeff(n,r,c)
implicit none
integer, intent(in) :: n,r
real(8), intent(inout) :: c
integer :: i,k
    c = 0
    if (r < 0 .or. r > n) then
        c = 0.0d0
        return
    end if

    k = min(r, n - r)
    c = 1.0d0

    do i = 1, k
        c = c * dble(n - k + i) / dble(i)
    end do
end subroutine bin_coeff


  subroutine fill_spin_strings(orbital_energies,n_orb,n_s,NRD_spin,RAS_el_array_spin,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,range1,range2,sizes,str_s)
  
  implicit none
  integer, intent(in):: n_s, n_orb
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in):: range1,range2
  integer, intent(in):: NRD_spin
  integer, intent(in):: RAS_el_array_spin(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(in)::sizes
  integer, intent(inout)::str_s(sizes,n_s)
  real(8), intent(in) :: orbital_energies(n_orb)
  integer :: space_sizes(n_RAS_spaces_occ+1+n_RAS_spaces_virt)
  integer, allocatable :: str_temp(:,:)
  integer, allocatable :: another_str_temp(:,:)
  integer :: siz_space, tot_siz_space
  integer :: n_el_temp
  integer :: orbital_index
  integer :: n_str_temp(n_RAS_spaces_occ+1+n_RAS_spaces_virt)
  integer :: j
  integer :: a, i , r ,v
  real(8) :: c
  !real(8) :: nCr_dp

  a=n_RAS_spaces_occ+1
  tot_siz_space = 0
  do i=range1,range2
      siz_space = 1
     orbital_index = 1
     j = 1
     n_el_temp = 0

   do r=1,n_RAS_spaces_occ
      call bin_coeff(RAS_space_occ(r),RAS_el_array_spin(i,r),c)
      space_sizes(j) = int(c)
      siz_space=siz_space*space_sizes(j)
      n_str_temp(j) = siz_space
      j = j + 1

   end do
   call bin_coeff(active_space,RAS_el_array_spin(i,a),c)
   space_sizes(j) = int(c)
   siz_space=siz_space*space_sizes(j)
   n_str_temp(j) = siz_space
   j = j + 1

   do v=1,n_RAS_spaces_virt
   call bin_coeff(RAS_space_virt(v),RAS_el_array_spin(i,a+v),c)
   space_sizes(j) = int(c)
   siz_space=siz_space*space_sizes(j)
   n_str_temp(j) = siz_space
   j = j + 1
   end do

   allocate(str_temp(siz_space,n_s))
   allocate(another_str_temp(siz_space,n_s))

   j = 0
   do r=1,n_RAS_spaces_occ
      call truegenerate(RAS_space_occ(r),RAS_el_array_spin(i,r),str_temp(1:space_sizes(r),1+n_el_temp:RAS_el_array_spin(i,r)+n_el_temp),orbital_index,space_sizes(r),RAS_el_array_spin(i,r))

      if (r .gt. 1) then
         call string_direct_product(str_temp(1:space_sizes(r),1+n_el_temp:RAS_el_array_spin(i,r)+n_el_temp),space_sizes(r),RAS_el_array_spin(i,r),str_temp(1:n_str_temp(j),1:n_el_temp),n_str_temp(j),n_el_temp,another_str_temp(1:n_str_temp(j+1),1:n_el_temp + RAS_el_array_spin(i,r)))
         str_temp = another_str_temp 

         
      end if 

      n_el_temp = n_el_temp + RAS_el_array_spin(i,r)
      j = j + 1
      orbital_index=orbital_index+RAS_space_occ(r)

   end do


      call truegenerate(active_space,RAS_el_array_spin(i,a),str_temp(1:space_sizes(a),1+n_el_temp:RAS_el_array_spin(i,a)+n_el_temp),orbital_index,space_sizes(a),RAS_el_array_spin(i,a))

      if (n_RAS_spaces_occ .ne. 0) then
      call string_direct_product(str_temp(1:space_sizes(a),1+n_el_temp:RAS_el_array_spin(i,a)+n_el_temp),space_sizes(a),RAS_el_array_spin(i,a),str_temp(1:n_str_temp(j),1:n_el_temp),n_str_temp(j),n_el_temp,another_str_temp(1:n_str_temp(j+1),1:n_el_temp+RAS_el_array_spin(i,a)))
      str_temp = another_str_temp

      end if


      n_el_temp = n_el_temp + RAS_el_array_spin(i,a)
      orbital_index=orbital_index+active_space
      j = j + 1



   do v=1,n_RAS_spaces_virt
      call truegenerate(RAS_space_virt(v),RAS_el_array_spin(i,a+v),str_temp(1:space_sizes(a+v),1+n_el_temp:RAS_el_array_spin(i,a+v)+n_el_temp),orbital_index,space_sizes(a+v),RAS_el_array_spin(i,a+v))

      
      call string_direct_product(str_temp(1:space_sizes(a+v),1+n_el_temp:RAS_el_array_spin(i,a+v)+n_el_temp),space_sizes(a+v),RAS_el_array_spin(i,a+v),str_temp(1:n_str_temp(j),1:n_el_temp),n_str_temp(j),n_el_temp,another_str_temp(1:n_str_temp(j+1),1:n_el_temp + RAS_el_array_spin(i,a+v)))
      
      str_temp = another_str_temp
      
      n_el_temp = n_el_temp + RAS_el_array_spin(i,a+v)
      orbital_index=orbital_index+RAS_space_virt(v)

      j = j + 1
   end do

   str_s(tot_siz_space+1:tot_siz_space+siz_space,:) = str_temp
 
   tot_siz_space = tot_siz_space + siz_space

   deallocate(str_temp)
   
   deallocate(another_str_temp)

   end do
   call sort_states(orbital_energies,n_orb,str_s,sizes,n_s)
end subroutine fill_spin_strings


subroutine truegenerate(n,k,options,orbital_index,s1,s2)
    implicit none
    integer,intent(in) :: n, k
    integer :: arr(n)
    integer, intent(inout) :: options(s1,s2)
    integer :: index
    integer, intent(in) :: orbital_index
    integer, intent(in) :: s1,s2
   
    ! Example values (you may change these or read from input)

    index = 0


    !options = 0
    call generate(arr, n, k, 1,options,index,s1,s2)
    options = options + orbital_index - 1
contains
recursive subroutine generate(arr, n, k, pos,options, index,s1,s2)
        implicit none
        integer, intent(inout) :: arr(n)
        integer, intent(in)    :: n, k, pos
        integer, intent(inout) :: options(s1,s2)
        integer, intent(inout) :: index
        integer, intent(in) :: s1,s2
        integer :: i,j
        
        ! If we filled all positions:
        if (pos > n) then
            if (count(arr .ne. 0) == k) then
                index = index + 1
                j = 1
                do i = 1, n
                    if (arr(i) .ne. 0) then
                        options(index,j) = arr(i)
                        j = j+1
                    end if
                end do
            end if
            return
        end if

        ! Option 1: place 0 here
        arr(pos) = 0
        call generate(arr, n, k, pos + 1,options,index,s1,s2)

        ! Option 2: place orbital number here
        arr(pos) = pos
    
        call generate(arr, n, k, pos + 1,options,index,s1,s2)

end subroutine generate
end subroutine truegenerate


subroutine string_direct_product(string_array1,nstr_1,n_el_1,string_array2,nstr_2,n_el_2,string_array_combined)
   implicit none
   integer, intent(in) ::string_array1(nstr_1,n_el_1)
   integer, intent(in):: string_array2(nstr_2,n_el_2)
   integer, intent(in) ::nstr_1,n_el_1,nstr_2,n_el_2
   integer, intent(inout):: string_array_combined(nstr_1*nstr_2,n_el_1+n_el_2)

   integer :: i,j
   integer :: n_s
   n_s = 0
   !allocate(string_array_combined(nstr_1*nstr_2,n_el_1+n_el_2))
   do j = 1,nstr_2
      do i = 1,nstr_1
         n_s = n_s + 1
         !string_array_combined(n_s,:) = [string_array1(nstr_1,:),string_array2(nstr_2,:)] do sprawdzenia
         string_array_combined(n_s,1:n_el_2) = string_array2(j,1:n_el_2)
         string_array_combined(n_s,n_el_2+1:n_el_1+n_el_2) = string_array1(i,1:n_el_1)
      end do
   end do
end subroutine string_direct_product


subroutine sort_states(orbital_energies,n_orb,spin_strings,n_spin_strings,n_spin)
implicit none
integer, intent(in) :: n_spin_strings,n_spin,n_orb
integer, intent(inout) :: spin_strings(n_spin_strings,n_spin)
real(8), intent(in) :: orbital_energies(n_orb)
real(8) :: temp_array(n_spin_strings), temp_energy
integer :: i,j
integer :: temp_state(n_spin)
if (n_spin_strings .ge. 2) then
temp_array(:) = 0
do i=1,n_spin_strings
   do j=1,n_spin
      temp_array(i) = temp_array(i) + orbital_energies(spin_strings(i,j))
   end do
end do
do i = 1, n_spin_strings-1
        do j = 1, n_spin_strings-i
            if (temp_array(j) > temp_array(j+1)) then
                temp_state   =  spin_strings(j,:)
                temp_energy = temp_array(j)

                spin_strings(j,:) = spin_strings(j+1,:)
                temp_array(j) = temp_array(j+1)

                spin_strings(j+1,:) = temp_state
                temp_array(j+1) = temp_energy
            end if
        end do
end do
end if
end subroutine sort_states


subroutine IS_EQUAL(vector1,vector2,length,indicator)
implicit none
    integer, intent(in) :: length
    logical, intent(inout) :: indicator
    integer, intent(in) :: vector1(length), vector2(length)
    integer :: i
    integer :: temp(length)
    indicator = .true.
    do i=1,length
        temp(i) = vector1(i) - vector2(i)
        if (temp(i) .ne. 0) then 
         indicator = .false.
         exit
        end if
    end do
end subroutine IS_EQUAL





subroutine annihilation(orbital,string_spin,n_spin,new_string_spin,sign)
implicit none
integer, intent(in) :: orbital,n_spin
integer, intent(in) :: string_spin(n_spin)
integer, intent(out) :: new_string_spin(n_spin-1)
integer, intent(out) :: sign
integer :: i,j,k
logical :: check
check = .false.
k = 1
do i=1,n_spin
   if (string_spin(i) .eq. orbital) then
      check = .true.
      exit
   end if 
end do
if (check .eqv. .true.) then
   do i=1,n_spin
   if (string_spin(i) .ne. orbital ) then
   new_string_spin(k) = string_spin(i)
   k = k + 1
   end if
   end do
   do j=1,n_spin
      if (string_spin(j) .eq. orbital) then
      sign = (-1)**(j-1)
      exit
      end if
   end do
else
   sign = 1
   new_string_spin(:) = 0
end if
end subroutine annihilation


subroutine fill_annihilation_results(n_spin_strings_m1,n_spin_strings,n_spin,spin_strings_m1,spin_strings,n_orb,spin_annihilation_matrix)
implicit none
integer, intent(in) :: n_spin_strings_m1,n_spin_strings,n_spin,n_orb
integer, intent(in) :: spin_strings_m1(n_spin_strings_m1,n_spin-1), spin_strings(n_spin_strings,n_spin)
integer, intent(inout) :: spin_annihilation_matrix(n_orb,n_spin_strings,2)
integer :: temp_string(n_spin-1)
integer :: sign, i, j, k
logical :: indicator
spin_annihilation_matrix(:,:,1) = 0
spin_annihilation_matrix(:,:,2) = 1
if (n_spin .gt. 1) then
do i=1,n_spin_strings
   do j=1,n_spin
      call annihilation(spin_strings(i,j),spin_strings(i,:),n_spin,temp_string,sign)
      do k=1,n_spin_strings_m1
         call IS_EQUAL(temp_string,spin_strings_m1(k,:),n_spin-1,indicator)
         if (indicator .eqv. .true.) then
            spin_annihilation_matrix(spin_strings(i,j),i,1) = k
            spin_annihilation_matrix(spin_strings(i,j),i,2) = sign
            exit
         end if
      end do
   end do
end do
else if (n_spin .eq. 1) then
do i=1,n_spin_strings
   spin_annihilation_matrix(spin_strings(i,1),i,1) = 1
end do
end if
end subroutine fill_annihilation_results





!***********************THIS WHOLE PART MIGHT BE OPTIMISED USING SMALL MATRICES**********************


subroutine fill_annihilation_creation_matrix(n_spin,n_spin_strings,spin_strings,n_orb,annihilation_creation_matrix)
implicit none
integer, intent(in) ::  n_spin, n_spin_strings, n_orb
integer, intent(in) :: spin_strings(n_spin_strings,n_spin)
integer, intent(inout) :: annihilation_creation_matrix(n_orb,n_orb,n_spin_strings,2)
integer :: i,j,k,l
integer :: temp_string1(n_spin-1), temp_string2(n_spin-1)
integer :: sign1,sign2
logical :: indicator
annihilation_creation_matrix(:,:,:,1) = 0
annihilation_creation_matrix(:,:,:,2) = 1
if (n_spin .gt. 1) then
do i=1,n_spin_strings
   do j=1,n_spin
      call annihilation(spin_strings(i,j),spin_strings(i,:),n_spin,temp_string1,sign1)
      do k=1,n_spin_strings
         do l=1,n_spin
            call annihilation(spin_strings(k,l),spin_strings(k,:),n_spin,temp_string2,sign2)
            call IS_EQUAL(temp_string1,temp_string2,n_spin-1,indicator)
            if (indicator .eqv. .true.) then
               annihilation_creation_matrix(spin_strings(i,j),spin_strings(k,l),k,1) = i
               annihilation_creation_matrix(spin_strings(i,j),spin_strings(k,l),k,2) = sign1*sign2
               exit
            end if
         end do
      end do
   end do
end do
else if (n_spin .eq. 1) then
do i=1,n_spin_strings
      do k=1,n_spin_strings
         annihilation_creation_matrix(spin_strings(i,1),spin_strings(k,1),k,1) = i
   end do
end do
end if
end subroutine fill_annihilation_creation_matrix


subroutine fill_nr_single_spin_hamiltonian(spin_strings,n_spin_strings,n_spin,n_orb,hopping,interaction,spin_annihilation_creation_matrix,spin_nr_hamiltonian)
integer, intent(in) :: n_spin_strings,n_spin,n_orb
integer, intent(in) :: spin_strings(n_spin_strings,n_spin), spin_annihilation_creation_matrix(n_orb,n_orb,n_spin_strings,2)
complex(8), intent(in) :: hopping(n_orb,n_orb), interaction(n_orb,n_orb,n_orb,n_orb)
complex(8), intent(out) :: spin_nr_hamiltonian(n_spin_strings,n_spin_strings)
integer :: j,p,q,r,s, temp_state1,temp_state2,temp_sign1,temp_sign2
spin_nr_hamiltonian(:,:) = 0
if (n_spin .gt. 0) then
do j=1,n_spin_strings
   do q=1,n_spin
      do p=1,n_orb
         temp_state1 = spin_annihilation_creation_matrix(p,spin_strings(j,q),j,1)
         temp_sign1 = spin_annihilation_creation_matrix(p,spin_strings(j,q),j,2)
         if (temp_state1 .ne. 0) then
            spin_nr_hamiltonian(temp_state1,j) = spin_nr_hamiltonian(temp_state1,j) + temp_sign1*hopping(p,spin_strings(j,q))
            do r=1,n_orb
               spin_nr_hamiltonian(temp_state1,j) = spin_nr_hamiltonian(temp_state1,j) + 0.5*temp_sign1*interaction(p,r,spin_strings(j,q),r)
               do s=1,n_spin
                  temp_state2 = spin_annihilation_creation_matrix(r,spin_strings(temp_state1,s),temp_state1,1)
                  temp_sign2 = temp_sign1*spin_annihilation_creation_matrix(r,spin_strings(temp_state1,s),temp_state1,2)
                  if (temp_state2 .ne. 0) then
                     spin_nr_hamiltonian(temp_state2,j) = spin_nr_hamiltonian(temp_state2,j) - 0.5*temp_sign2*interaction(r,p,spin_strings(j,q),spin_strings(temp_state1,s)) 
                  end if
               end do
            end do
         end if
      end do
   end do
end do
else
spin_nr_hamiltonian(:,:) = 0
end if
end subroutine fill_nr_single_spin_hamiltonian


subroutine nr_matrix_vector_product(alpha_hamiltonian,strings_alpha,n_strings_alpha,n_alpha,beta_hamiltonian,strings_beta,n_strings_beta,n_beta,interaction_mix,n_orb,alpha_annihilation_creation_matrix,beta_annihilation_creation_matrix,vector,vector_new)
integer, intent (in) :: n_strings_alpha,n_strings_beta,n_orb, n_alpha, n_beta
integer, intent (in) :: alpha_annihilation_creation_matrix(n_orb,n_orb,n_strings_alpha,2), beta_annihilation_creation_matrix(n_orb,n_orb,n_strings_beta,2), strings_alpha(n_strings_alpha,n_alpha), strings_beta(n_strings_beta,n_beta)
complex(8), intent (in) :: alpha_hamiltonian(n_strings_alpha,n_strings_alpha), beta_hamiltonian(n_strings_beta,n_strings_beta), interaction_mix(n_orb,n_orb,n_orb,n_orb)
complex(8), intent (inout) :: vector_new(n_strings_alpha*n_strings_beta),vector(n_strings_alpha*n_strings_beta)
complex(8) :: temp_table(n_strings_alpha,n_strings_beta), new_table(n_strings_alpha,n_strings_beta)
integer :: mu, i,j,k,l, p,q,r,s,temp_state1,temp_state2,temp_sign1,temp_sign2
new_table(:,:) = 0
temp_table(:,:) = 0
vector_new(:) = 0
do i=1,n_strings_alpha
   mu = (i-1)*n_strings_beta+1
   do j=1,n_strings_beta
      temp_table(i,j) = vector(mu)
      mu = mu + 1
   end do
end do
do i=1,n_strings_alpha
   do j=1,n_strings_beta
      do k=1,n_strings_alpha
         new_table(i,j) = new_table(i,j) + temp_table(k,j)*alpha_hamiltonian(i,k)
         if (new_table(i,j) /= new_table(i,j)) then
          print *, 'NaN pojawil sie dla', i,j,k 
          stop
         end if
      end do
      do l=1,n_strings_beta
         new_table(i,j) = new_table(i,j) + temp_table(i,l)*beta_hamiltonian(j,l)
         if (new_table(i,j) /= new_table(i,j)) then
          print *, 'NaN pojawil sie dla', i,j,k,l 
          write(*,*) temp_table(i,l)
          write(*,*) beta_hamiltonian(j,l)
          stop
         end if
      end do
      do p=1,n_alpha
         do r=1,n_orb
            temp_state1 = alpha_annihilation_creation_matrix(r,strings_alpha(i,p),i,1)
            if (temp_state1 .ne. 0) then
               temp_sign1 = alpha_annihilation_creation_matrix(r,strings_alpha(i,p),i,2)
               do q=1,n_beta
                  do s=1,n_orb
                     temp_state2 = beta_annihilation_creation_matrix(s,strings_beta(j,q),j,1)
                     if (temp_state2 .ne. 0) then
                        temp_sign2 = temp_sign1*beta_annihilation_creation_matrix(s,strings_beta(j,q),j,2)
                        new_table(i,j) = new_table(i,j) + temp_sign2*temp_table(temp_state1,temp_state2)*interaction_mix(strings_alpha(i,p),strings_beta(j,q),r,s)
                        if (new_table(i,j) /= new_table(i,j)) then
                           print *, 'NaN pojawil sie dla', i,j,p,r,q,s
                           stop
                        end if
                     end if 
                  end do
               end do   
            end if
         end do
      end do
   end do
end do
do i=1,n_strings_alpha
   mu = (i-1)*n_strings_beta+1
   do j=1,n_strings_beta
      vector_new(mu) = new_table(i,j)
      mu = mu + 1
   end do
end do
end subroutine nr_matrix_vector_product


subroutine rel_matrix_vector_product(dane,vector,vector_new)
type(parameters), intent(inout) :: dane
complex(8), intent(inout) :: vector(dane%n_strings_alpha*dane%n_strings_beta+dane%n_strings_alpha_p1*dane%n_strings_beta_m1+dane%n_strings_alpha_m1*dane%n_strings_beta_p1)
complex(8), intent (inout) :: vector_new(dane%n_strings_alpha*dane%n_strings_beta+dane%n_strings_alpha_p1*dane%n_strings_beta_m1+dane%n_strings_alpha_m1*dane%n_strings_beta_p1)
complex(8) ::  vector_new_1(dane%n_strings_alpha*dane%n_strings_beta),vector_new_2(dane%n_strings_alpha_p1*dane%n_strings_beta_m1),vector_new_3(dane%n_strings_alpha_m1*dane%n_strings_beta_p1)
complex(8) :: temp_table1(dane%n_strings_alpha,dane%n_strings_beta), new_table1(dane%n_strings_alpha,dane%n_strings_beta),temp_table2(dane%n_strings_alpha_p1,dane%n_strings_beta_m1), new_table2(dane%n_strings_alpha_p1,dane%n_strings_beta_m1),temp_table3(dane%n_strings_alpha_m1,dane%n_strings_beta_p1), new_table3(dane%n_strings_alpha_m1,dane%n_strings_beta_p1)
integer :: i,j,k,l,p,q,mu, temp_state1, temp_state2, temp_sign1, temp_sign2
temp_table1(:,:) = 0
temp_table2(:,:) = 0
temp_table3(:,:) = 0 
vector_new(:) = 0
vector_new_1(:) = 0
vector_new_2(:) = 0
vector_new_3(:) = 0
call nr_matrix_vector_product(dane%alpha_hamiltonian,dane%str_a,dane%n_strings_alpha,dane%n_alpha,dane%beta_hamiltonian,dane%str_b,dane%n_strings_beta,dane%n_beta,dane%interaction_mix,dane%n_orb,dane%alpha_annihilation_creation_matrix,dane%beta_annihilation_creation_matrix,vector(dane%size_tot(1,1):dane%size_tot(1,2)),vector_new_1)
call nr_matrix_vector_product(dane%alpha_hamiltonian_p1,dane%str_a_p1,dane%n_strings_alpha_p1,dane%n_alpha+1,dane%beta_hamiltonian_m1,dane%str_b_m1,dane%n_strings_beta_m1,dane%n_beta-1,dane%interaction_mix,dane%n_orb,dane%alpha_annihilation_creation_matrix_p1,dane%beta_annihilation_creation_matrix_m1,vector(dane%size_tot(2,1):dane%size_tot(1,2)+dane%size_tot(2,2)),vector_new_2)
call nr_matrix_vector_product(dane%alpha_hamiltonian_m1,dane%str_a_m1,dane%n_strings_alpha_m1,dane%n_alpha-1,dane%beta_hamiltonian_p1,dane%str_b_p1,dane%n_strings_beta_p1,dane%n_beta+1,dane%interaction_mix,dane%n_orb,dane%alpha_annihilation_creation_matrix_m1,dane%beta_annihilation_creation_matrix_p1,vector(dane%size_tot(3,1):dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)),vector_new_3)
do i = 1, dane%n_strings_alpha*dane%n_strings_beta
   if (vector_new_1(i) /= vector_new_1(i)) then
      print *, "NaN pojawił się w pętli 1 dla i =", i
      stop
   end if
end do
do i = 1, dane%n_strings_alpha_p1*dane%n_strings_beta_m1
   if (vector_new_2(i) /= vector_new_2(i)) then
      print *, "NaN pojawił się w pętli 2 dla i =", i
      stop
   end if
end do
do i = 1, dane%n_strings_alpha_m1*dane%n_strings_beta_p1
   if (vector_new_3(i) /= vector_new_3(i)) then
      print *, "NaN pojawił się w pętli 3 dla i =", i
      stop
   end if
end do
new_table1(:,:) = 0
new_table2(:,:) = 0
new_table3(:,:) = 0

do i=1,dane%n_strings_alpha
   mu = (i-1)*dane%n_strings_beta+1
   do j=1,dane%n_strings_beta
      temp_table1(i,j) = vector(mu)
      mu = mu + 1
   end do
end do


do i=1,dane%n_strings_alpha_p1
   mu = (i-1)*dane%n_strings_beta_m1+dane%size_tot(2,1)
   do j=1,dane%n_strings_beta_m1
      temp_table2(i,j) = vector(mu)
      mu = mu + 1
   end do
end do



do i=1,dane%n_strings_alpha_m1
   mu = (i-1)*dane%n_strings_beta_p1+dane%size_tot(3,1)
   do j=1,dane%n_strings_beta_p1
      temp_table3(i,j) = vector(mu)
      mu = mu + 1
   end do
end do


do j=1,dane%n_strings_beta
      do p=1,dane%n_beta
         temp_state1 = dane%beta_annihilation_matrix(dane%str_b(j,p),j,1)
         if (temp_state1 .ne. 0) then
            temp_sign1 = dane%beta_annihilation_matrix(dane%str_b(j,p),j,2)
            do k=1,dane%n_strings_alpha_p1
               do q=1,dane%n_alpha+1
                  temp_state2 = dane%alpha_annihilation_matrix_p1(dane%str_a_p1(k,q),k,1)
                  if (temp_state2 .ne. 0) then
                     temp_sign2 = temp_sign1 * dane%alpha_annihilation_matrix_p1(dane%str_a_p1(k,q),k,2)
                     new_table1(temp_state2,j) = new_table1(temp_state2,j) +  temp_sign2*conjg(dane%hso_ab(dane%str_a_p1(k,q),dane%str_b(j,p)))*temp_table2(k,temp_state1)*(-1)**dane%n_alpha
                  end if
               end do
            end do
         end if
      end do
end do

do i=1,dane%n_strings_alpha
   do q=1,dane%n_alpha
      temp_state1 = dane%alpha_annihilation_matrix(dane%str_a(i,q),i,1) 
      if (temp_state1 .ne. 0) then
         temp_sign1 = dane%alpha_annihilation_matrix(dane%str_a(i,q),i,2)
         do l=1,dane%n_strings_beta_p1
            do p=1,dane%n_beta+1
               temp_state2 = dane%beta_annihilation_matrix_p1(dane%str_b_p1(l,p),l,1)
               if (temp_state2 .ne. 0) then
                  temp_sign2 = temp_sign1*dane%beta_annihilation_matrix_p1(dane%str_b_p1(l,p),l,2)
                  new_table1(i,temp_state2) = new_table1(i,temp_state2) + temp_sign2*dane%hso_ab(dane%str_a(i,q),dane%str_b_p1(l,p))*temp_table3(temp_state1,l)*(-1)**(dane%n_alpha-1)
               end if
            end do
         end do 
      end if     
   end do
end do

do i=1,dane%n_strings_alpha_p1
   do q=1,dane%n_alpha+1
      temp_state1 = dane%alpha_annihilation_matrix_p1(dane%str_a_p1(i,q),i,1)
      if (temp_state1 .ne. 0) then
        temp_sign1 = dane%alpha_annihilation_matrix_p1(dane%str_a_p1(i,q),i,2)
        do l=1,dane%n_strings_beta
         do p=1,dane%n_beta
            temp_state2 = dane%beta_annihilation_matrix(dane%str_b(l,p),l,1)
            if (temp_state2 .ne. 0) then
               temp_sign2 = temp_sign1*dane%beta_annihilation_matrix(dane%str_b(l,p),l,2)
               new_table2(i,temp_state2) = new_table2(i,temp_state2) + temp_sign2*dane%hso_ab(dane%str_a_p1(i,q),dane%str_b(l,p))*temp_table1(temp_state1,l)*(-1)**dane%n_alpha 
            end if
         end do
        end do        
      end if
   end do
end do

do j=1,dane%n_strings_beta_p1
   do p=1,dane%n_beta+1
      temp_state1 = dane%beta_annihilation_matrix_p1(dane%str_b_p1(j,p),j,1)
      if (temp_state1 .ne. 0) then
         temp_sign1 = dane%beta_annihilation_matrix_p1(dane%str_b_p1(j,p),j,2)
         do k=1,dane%n_strings_alpha
            do q=1,dane%n_alpha
               temp_state2 = dane%alpha_annihilation_matrix(dane%str_a(k,q),k,1)
               if (temp_state2 .ne. 0) then
                  temp_sign2 = temp_sign1*dane%alpha_annihilation_matrix(dane%str_a(k,q),k,2)
                  new_table3(temp_state2,j) = new_table3(temp_state2,j) + temp_sign2*conjg(dane%hso_ab(dane%str_a(k,q),dane%str_b_p1(j,p)))*temp_table1(k,temp_state1)*(-1)**(dane%n_alpha-1)
               end if
            end do
         end do
      end if
   end do
end do

do i=1,dane%n_strings_alpha
   mu = (i-1)*dane%n_strings_beta+1
   do j=1,dane%n_strings_beta
      vector_new_1(mu) = vector_new_1(mu) + new_table1(i,j)
      mu = mu + 1
   end do
end do

do i=1,dane%n_strings_alpha_p1
   mu = (i-1)*dane%n_strings_beta_m1+1
   do j=1,dane%n_strings_beta_m1
      vector_new_2(mu) = vector_new_2(mu) + new_table2(i,j)
      mu = mu + 1
   end do
end do

do i=1,dane%n_strings_alpha_m1
   mu = (i-1)*dane%n_strings_beta_p1+1
   do j=1,dane%n_strings_beta_p1
      vector_new_3(mu) = vector_new_3(mu) + new_table3(i,j)
      mu = mu + 1
   end do
end do
vector_new(1:dane%size_tot(1,2)) = vector_new_1(:)
vector_new(dane%size_tot(2,1):dane%size_tot(3,1)-1) = vector_new_2(:)
vector_new(dane%size_tot(3,1):dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)) = vector_new_3(:)

end subroutine rel_matrix_vector_product

  subroutine matrix_vector_product(dane, vector, vector_new)
    implicit none
    type(parameters), intent(inout) :: dane
    complex(8), intent(inout) :: vector(:)
    complex(8), intent(inout) :: vector_new(:)
    !integer, intent(in) :: n
    !n = dane%size_tot(1,2) + dane%size_tot(2,2) + dane%size_tot(3,2)
    vector_new(:) = 0
    if (dane%relativistic) then
      call rel_matrix_vector_product(dane, vector, vector_new)
    else

      call nr_matrix_vector_product(dane%alpha_hamiltonian, dane%str_a, dane%n_strings_alpha, dane%n_alpha, dane%beta_hamiltonian, dane%str_b, dane%n_strings_beta, dane%n_beta, dane%interaction_mix, dane%n_orb, dane%alpha_annihilation_creation_matrix, dane%beta_annihilation_creation_matrix, vector, vector_new)
    end if
  end subroutine matrix_vector_product




subroutine generate_diagonal_elements(alpha_hamiltonian,strings_alpha,n_strings_alpha,n_alpha,beta_hamiltonian,strings_beta,n_strings_beta,n_beta,interaction,n_orb,diagonal)
integer, intent (in) :: n_strings_alpha,n_strings_beta,n_orb, n_alpha, n_beta
integer, intent (in) :: strings_alpha(n_strings_alpha,n_alpha), strings_beta(n_strings_beta,n_beta)
complex(8), intent (in) :: alpha_hamiltonian(n_strings_alpha,n_strings_alpha), beta_hamiltonian(n_strings_beta,n_strings_beta), interaction(n_orb,n_orb,n_orb,n_orb)
complex(8), intent(out) :: diagonal(n_strings_alpha*n_strings_beta)
integer :: mu,i,j,p,q
complex(8) :: temp_diagonal(n_strings_alpha*n_strings_beta)
do mu=1,n_strings_alpha*n_strings_beta
   i = (mu-1)/n_strings_beta+1
   j = mod(mu-1,n_strings_beta)+1
   temp_diagonal(mu) = alpha_hamiltonian(i,i)+beta_hamiltonian(j,j)
   do p=1,n_alpha
      do q=1,n_beta
         temp_diagonal(mu) = temp_diagonal(mu) + interaction(strings_alpha(i,p),strings_beta(j,q),strings_alpha(i,p),strings_beta(j,q))
      end do
   end do
end do
diagonal = temp_diagonal
end subroutine generate_diagonal_elements

subroutine generate_params(dane,relativistic,n_orb,n_alpha,n_beta,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,RAS_space_virt,active_space,excit_array,orbital_energies,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix, hso_ab, hso_ba)
   implicit none
   type(parameters), intent(inout) :: dane
   logical, intent(in) :: relativistic
   integer, intent(in) :: n_orb, n_alpha, n_beta, n_RAS_spaces_occ, n_RAS_spaces_virt, active_space
   integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt),excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
   real(8), intent(in) :: orbital_energies(n_orb)
   complex(8), intent(inout) :: hopping_alpha(n_orb,n_orb),hopping_beta(n_orb,n_orb),interaction_alpha(n_orb,n_orb,n_orb,n_orb),interaction_beta(n_orb,n_orb,n_orb,n_orb), interaction_mix(n_orb,n_orb,n_orb,n_orb), hso_ab(n_orb,n_orb), hso_ba(n_orb,n_orb)
   integer :: max,n_spaces,n_combinations, temp, n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta
   integer :: NRDa(3,2),NRDb(3,2),sizea(3,2),sizeb(3,2)
   integer :: i,j
   integer, allocatable :: all_combinations(:,:),RAS_el_array_alpha(:,:),RAS_el_array_beta(:,:)
   call check_orbital_space_declarations(n_orb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)

   call check_hamiltonian(relativistic,n_orb,hopping_alpha,hopping_beta,interaction_alpha,interaction_beta,interaction_mix,hso_ab,hso_ba)

   allocate(dane%interaction_mix(n_orb,n_orb,n_orb,n_orb))
   allocate(dane%hso_ab(n_orb,n_orb))
   allocate(dane%hso_ba(n_orb,n_orb))

   dane%relativistic = relativistic
   dane%interaction_mix = interaction_mix
   dane%hso_ab = hso_ab
   dane%hso_ba = hso_ba
   dane%n_alpha = n_alpha
   dane%n_beta = n_beta
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
   call count_spin_distributions(dane%relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta)
   allocate(RAS_el_array_alpha(NRD_spin_alpha,n_spaces))
   allocate(RAS_el_array_beta(NRD_spin_beta,n_spaces))

   write(*,*) "NRD_spin_alpha: ",NRD_spin_alpha
   write(*,*) "NRD_spin_beta: ",NRD_spin_beta
   ! here we fill RAS_el_array_spin tables and calculate NRDa, NRDb
   !NRD(:,1) denotes beginning of spin block, NRD(:,2) denotes end
   call find_spin_distributions(dane%relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,NRDa,NRDb)
   write(*,*) "NRDa:"
   write(*,*) NRDa(1,:)
   write(*,*) NRDa(2,:)
   write(*,*) NRDa(3,:)

   write(*,*) "NRDb: "
   write(*,*) NRDb(1,:)
   write(*,*) NRDb(2,:)
   write(*,*) NRDb(3,:)
   ! here we calculate sizea,sizeb,dane%size_tot
   call calculate_space_size(dane%relativistic,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,dane%size_tot)

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
   else if ((dane%relativistic .eqv. .true.) .and. (dane%size_tot(2,2)+dane%size_tot(3,2) .eq. 0)) then
    write(*,*) "RELATIVISTIC BLOCKS ARE EMPTY"
    write(*,*) "QUITTING THE PROGRAM"
    stop
   end if
! here we generate strings in first block (alpha and beta separately)

if (n_alpha .gt. 0)  then
  allocate(dane%str_a(sizea(1,2),n_alpha))
  call fill_spin_strings(orbital_energies,n_orb,n_alpha,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(1,1),NRDa(1,2),sizea(1,2),dane%str_a)
else if (n_alpha .eq. 0)  then
  allocate(dane%str_a(sizea(1,2),1))
  dane%str_a(sizea(1,2),:) = 0 !vacuum  
end if

if (n_beta .gt. 0) then
  allocate(dane%str_b(sizeb(1,2),n_beta))
  call fill_spin_strings(orbital_energies,n_orb,n_beta,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(1,1),NRDb(1,2),sizeb(1,2),dane%str_b)
else if  (n_beta .eq. 0)  then
  allocate(dane%str_b(sizeb(1,2),1))
  dane%str_b(sizeb(1,2),:) = 0 !vacuum  
end if

!here we generate creation-annihilation matrices for first block

allocate(dane%alpha_annihilation_creation_matrix(n_orb,n_orb,sizea(1,2),2))  
call fill_annihilation_creation_matrix(n_alpha,sizea(1,2),dane%str_a,n_orb,dane%alpha_annihilation_creation_matrix)

allocate(dane%beta_annihilation_creation_matrix(n_orb,n_orb,sizeb(1,2),2))
call fill_annihilation_creation_matrix(n_beta,sizeb(1,2),dane%str_b,n_orb,dane%beta_annihilation_creation_matrix)

!here we generate single spin part of hamiltonian for first block

allocate(dane%alpha_hamiltonian(sizea(1,2),sizea(1,2)))
allocate(dane%beta_hamiltonian(sizeb(1,2),sizeb(1,2)))

call fill_nr_single_spin_hamiltonian(dane%str_a,sizea(1,2),n_alpha,n_orb,hopping_alpha,interaction_alpha,dane%alpha_annihilation_creation_matrix,dane%alpha_hamiltonian)
call fill_nr_single_spin_hamiltonian(dane%str_b,sizeb(1,2),n_beta,n_orb,hopping_beta,interaction_beta,dane%beta_annihilation_creation_matrix,dane%beta_hamiltonian)

! here we prepare for matrix-vector multiplication


allocate(dane%diag(dane%size_tot(1,2)+dane%size_tot(2,2)+dane%size_tot(3,2)))

call generate_diagonal_elements(dane%alpha_hamiltonian,dane%str_a,sizea(1,2),n_alpha,dane%beta_hamiltonian,dane%str_b,sizeb(1,2),n_beta,dane%interaction_mix,n_orb,dane%diag(1:dane%size_tot(1,2)))
!here we repeat for other blocks



!here we generate strings

    if (n_alpha .gt. 1)  then
        allocate(dane%str_a_m1(sizea(3,2),n_alpha-1))
        if (sizea(3,2) .ne. 0) then
            call fill_spin_strings(orbital_energies,n_orb,n_alpha-1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(3,1),NRDa(3,2),sizea(3,2),dane%str_a_m1)
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
            call fill_spin_strings(orbital_energies,n_orb,n_beta-1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(3,1),NRDb(3,2),sizeb(3,2),dane%str_b_m1)
        end if
    else if (n_beta .eq. 1) then
        allocate(dane%str_b_m1(sizeb(3,2),1))
        dane%str_b_m1(:,:) = 0 !vacuum
    else 
        allocate(dane%str_b_m1(sizeb(3,2),0))
    end if
    
    allocate(dane%str_a_p1(sizea(2,2),n_alpha+1))
    if (sizea(2,2) .ne. 0) then
     call fill_spin_strings(orbital_energies,n_orb,n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),dane%str_a_p1)
    end if
      
    allocate(dane%str_b_p1(sizeb(2,2),n_beta+1)) 
    if (sizeb(2,2) .ne. 0) then
        call fill_spin_strings(orbital_energies,n_orb,n_beta+1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(2,1),NRDb(2,2),sizeb(2,2),dane%str_b_p1)
    end if

!here we generate annihilation, creation-annihilation and single spin hamiltonians
    allocate(dane%alpha_annihilation_matrix(n_orb,sizea(1,2),2))
    allocate(dane%alpha_annihilation_creation_matrix_m1(n_orb,n_orb,sizea(3,2),2))
    allocate(dane%alpha_hamiltonian_m1(sizea(3,2),sizea(3,2)))
    if (sizea(3,2) .ne. 0) then      
        call fill_annihilation_results(sizea(3,2),sizea(1,2),n_alpha,dane%str_a_m1,dane%str_a,n_orb,dane%alpha_annihilation_matrix)
        call fill_annihilation_creation_matrix(n_alpha-1,sizea(3,2),dane%str_a_m1,n_orb,dane%alpha_annihilation_creation_matrix_m1)
        call fill_nr_single_spin_hamiltonian(dane%str_a_m1,sizea(3,2),n_alpha-1,n_orb,hopping_alpha,interaction_alpha,dane%alpha_annihilation_creation_matrix_m1,dane%alpha_hamiltonian_m1)
    end if
  
    allocate(dane%beta_annihilation_matrix(n_orb,sizeb(1,2),2))
    allocate(dane%beta_annihilation_creation_matrix_m1(n_orb,n_orb,sizeb(3,2),2))
    allocate(dane%beta_hamiltonian_m1(sizeb(3,2),sizeb(3,2)))
    
    if (sizeb(3,2) .ne. 0) then      
        call fill_annihilation_results(sizeb(3,2),sizeb(1,2),n_beta,dane%str_b_m1,dane%str_b,n_orb,dane%beta_annihilation_matrix)
        call fill_annihilation_creation_matrix(n_beta-1,sizeb(3,2),dane%str_b_m1,n_orb,dane%beta_annihilation_creation_matrix_m1)
        call fill_nr_single_spin_hamiltonian(dane%str_b_m1,sizeb(3,2),n_beta-1,n_orb,hopping_beta,interaction_beta,dane%beta_annihilation_creation_matrix_m1,dane%beta_hamiltonian_m1)
    end if
 

    allocate(dane%alpha_annihilation_matrix_p1(n_orb,sizea(2,2),2))
    allocate(dane%alpha_annihilation_creation_matrix_p1(n_orb,n_orb,sizea(2,2),2))
    allocate(dane%alpha_hamiltonian_p1(sizea(2,2),sizea(2,2)))
    
    if (sizea(2,2) .ne. 0) then
        call fill_spin_strings(orbital_energies,n_orb,n_alpha+1,NRD_spin_alpha,RAS_el_array_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa(2,1),NRDa(2,2),sizea(2,2),dane%str_a_p1)
        call fill_annihilation_results(sizea(1,2),sizea(2,2),n_alpha+1,dane%str_a,dane%str_a_p1,n_orb,dane%alpha_annihilation_matrix_p1)
        call fill_annihilation_creation_matrix(n_alpha+1,sizea(2,2),dane%str_a_p1,n_orb,dane%alpha_annihilation_creation_matrix_p1)
        call fill_nr_single_spin_hamiltonian(dane%str_a_p1,sizea(2,2),n_alpha+1,n_orb,hopping_alpha,interaction_alpha,dane%alpha_annihilation_creation_matrix_p1,dane%alpha_hamiltonian_p1)
    end if

    allocate(dane%beta_annihilation_matrix_p1(n_orb,sizeb(2,2),2))  
    allocate(dane%beta_annihilation_creation_matrix_p1(n_orb,n_orb,sizeb(2,2),2))
    allocate(dane%beta_hamiltonian_p1(sizeb(2,2),sizeb(2,2)))   

    if (sizeb(2,2) .ne. 0) then
        call fill_spin_strings(orbital_energies,n_orb,n_beta+1,NRD_spin_beta,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDb(2,1),NRDb(2,2),sizeb(2,2),dane%str_b_p1)
        call fill_annihilation_results(sizeb(1,2),sizeb(2,2),n_beta+1,dane%str_b,dane%str_b_p1,n_orb,dane%beta_annihilation_matrix_p1)
        call fill_annihilation_creation_matrix(n_beta+1,sizeb(2,2),dane%str_b_p1,n_orb,dane%beta_annihilation_creation_matrix_p1)
        call fill_nr_single_spin_hamiltonian(dane%str_b_p1,sizeb(2,2),n_beta+1,n_orb,hopping_beta,interaction_beta,dane%beta_annihilation_creation_matrix_p1,dane%beta_hamiltonian_p1)
    end if
! here we generate hamiltonian diagonal elements
    if (dane%size_tot(2,2) .ne. 0) then
        call generate_diagonal_elements(dane%alpha_hamiltonian_p1,dane%str_a_p1,sizea(2,2),n_alpha+1,dane%beta_hamiltonian_m1,dane%str_b_m1,sizeb(3,2),n_beta-1,dane%interaction_mix,n_orb,dane%diag(dane%size_tot(2,1):dane%size_tot(2,2)+dane%size_tot(1,2)))
    end if

    if (dane%size_tot(3,2) .ne. 0) then 
        call generate_diagonal_elements(dane%alpha_hamiltonian_m1,dane%str_a_m1,sizea(3,2),n_alpha-1,dane%beta_hamiltonian_p1,dane%str_b_p1,sizeb(2,2),n_beta+1,dane%interaction_mix,n_orb,dane%diag(dane%size_tot(3,1):dane%size_tot(3,2)+dane%size_tot(2,2)+dane%size_tot(1,2)))
    end if
dane%n_strings_alpha = sizea(1,2)
dane%n_strings_beta = sizeb(1,2)
dane%n_orb = n_orb
dane%n_alpha = n_alpha
dane%n_beta = n_beta
dane%n_strings_alpha_p1 = sizea(2,2)
dane%n_strings_beta_p1 = sizeb(2,2)
dane%n_strings_alpha_m1 = sizea(3,2)
dane%n_strings_beta_m1 = sizeb(3,2)
end subroutine generate_params
!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_electron(orbital,state,new_state,dane,new_dane,spin)
type(parameters), intent(in) :: dane, new_dane
integer, intent(in) :: orbital
integer,intent(in) :: spin
complex(8), intent(in) :: state(dane%size_tot(1,2) + dane%size_tot(2,2) + dane%size_tot(3,2))
complex(8), intent(out) :: new_state(new_dane%size_tot(1,2) + new_dane%size_tot(2,2) + new_dane%size_tot(3,2))
complex(8) :: temp_table1(dane%n_strings_alpha,dane%n_strings_beta), temp_table2(dane%n_strings_alpha_p1,dane%n_strings_beta_m1), temp_table3(dane%n_strings_alpha_m1,dane%n_strings_beta_p1)
complex(8) :: new_temp_table1(new_dane%n_strings_alpha,new_dane%n_strings_beta), new_temp_table2(new_dane%n_strings_alpha_p1,new_dane%n_strings_beta_m1), new_temp_table3(new_dane%n_strings_alpha_m1,new_dane%n_strings_beta_p1)
integer :: i,j,mu
   new_state(:) = 0
   new_temp_table1(:,:)=0
   new_temp_table2(:,:)=0
   new_temp_table3(:,:)=0
   if (((dane%n_alpha-new_dane%n_alpha .eq. -1) .and. ((dane%n_beta .eq. new_dane%n_beta) .and. (spin .eq. 1)))) then
      continue
   else if (((dane%n_beta-new_dane%n_beta .eq. -1) .and. ((dane%n_alpha .eq. new_dane%n_alpha) .and. (spin .eq. -1)))) then
      continue
   else
      write(*,*) "WRONG NUMBER OF ELECTRONS"
      stop
   end if
! TBD: ZABEZPIECZENIA JESLI JAKIES BLOKI PUSTE
   do i=1,dane%n_strings_alpha
      mu = (i-1)*dane%n_strings_beta+1
      do j=1,dane%n_strings_beta
         temp_table1(i,j) = state(mu)
         mu = mu + 1
      end do
   end do
 if (dane%relativistic .eqv. .true.) then
   do i=1,dane%n_strings_alpha_p1
      mu = (i-1)*dane%n_strings_beta_m1+1
      do j=1,dane%n_strings_beta_m1
         temp_table2(i,j) = state(mu+dane%size_tot(1,2))
         mu = mu + 1
      end do
   end do
   
   do i=1,dane%n_strings_alpha_m1
      mu = (i-1)*dane%n_strings_beta_p1+1
      do j=1,dane%n_strings_beta_p1
         temp_table3(i,j) = state(mu+dane%size_tot(1,2)+dane%size_tot(2,2))
         mu = mu + 1
      end do
   end do
  end if
   if (spin .eq. 1) then
      do i=1,new_dane%n_strings_alpha
          if (new_dane%alpha_annihilation_matrix(orbital,i,1) .ne. 0) then
            do j=1,new_dane%n_strings_beta
               new_temp_table1(i,j) = new_dane%alpha_annihilation_matrix(orbital,i,2)*temp_table1(new_dane%alpha_annihilation_matrix(orbital,i,1),j)
            end do
          end if
      end do
      
      if (dane%relativistic .eqv. .true.) then
      do i=1,new_dane%n_strings_alpha_p1
           if (new_dane%alpha_annihilation_matrix_p1(orbital,i,1) .ne. 0) then
             do j=1,new_dane%n_strings_beta_m1
               new_temp_table2(i,j) = new_dane%alpha_annihilation_matrix_p1(orbital,i,2)*temp_table2(new_dane%alpha_annihilation_matrix_p1(orbital,i,1),j)
             end do
           end if
      end do  

      do i=1,new_dane%n_strings_alpha_m1
           if (dane%alpha_annihilation_matrix(orbital,i,1) .ne. 0) then
            do j=1,new_dane%n_strings_beta_p1
               new_temp_table3(i,j) = dane%alpha_annihilation_matrix(orbital,i,2)*temp_table3(dane%alpha_annihilation_matrix(orbital,i,1),j)
            end do
           end if
      end do 
      end if
   else
         do j=1,new_dane%n_strings_beta
           if (new_dane%beta_annihilation_matrix(orbital,j,1) .ne. 0) then
            do i=1,new_dane%n_strings_alpha
               new_temp_table1(i,j) = new_dane%beta_annihilation_matrix(orbital,j,2)*temp_table1(i,new_dane%beta_annihilation_matrix(orbital,j,1))*(-1)**new_dane%n_alpha
            end do
           end if
         end do
      
      
       if (dane%relativistic .eqv. .true.) then
         do j=1,new_dane%n_strings_beta_m1
           if (dane%beta_annihilation_matrix(orbital,j,1) .ne. 0) then
            do i=1,new_dane%n_strings_alpha_p1
               new_temp_table2(i,j) = dane%beta_annihilation_matrix(orbital,j,2)*temp_table2(i,dane%beta_annihilation_matrix(orbital,j,1))*(-1)**(new_dane%n_alpha+1)
            end do
           end if
         end do


         do j=1,new_dane%n_strings_beta_p1
           if (new_dane%beta_annihilation_matrix_p1(orbital,j,1) .ne. 0) then
            do i=1,new_dane%n_strings_alpha_m1
               new_temp_table3(i,j) = new_dane%beta_annihilation_matrix_p1(orbital,j,2)*temp_table3(i,new_dane%beta_annihilation_matrix_p1(orbital,j,1))*(-1)**(new_dane%n_alpha-1)
            end do
           end if
         end do
         end if

   end if
do i=1,new_dane%n_strings_alpha
   mu = (i-1)*new_dane%n_strings_beta+1
   do j=1,new_dane%n_strings_beta
      new_state(mu) =  new_temp_table1(i,j)
      mu = mu + 1
   end do
end do
if (dane%relativistic .eqv. .true.) then
do i=1,new_dane%n_strings_alpha_p1
   mu = (i-1)*new_dane%n_strings_beta_m1+1
   do j=1,new_dane%n_strings_beta_m1
      new_state(mu+new_dane%size_tot(1,2)) =  new_temp_table2(i,j)
      mu = mu + 1
   end do
end do

do i=1,new_dane%n_strings_alpha_m1
   mu = (i-1)*new_dane%n_strings_beta_p1+1
   do j=1,new_dane%n_strings_beta_p1
      new_state(mu+new_dane%size_tot(1,2)+new_dane%size_tot(2,2)) =  new_temp_table3(i,j)
      mu = mu + 1
   end do
end do
end if
end subroutine create_electron

subroutine create_hole(orbital,state,new_state,dane,new_dane,spin)
type(parameters), intent(in) :: dane, new_dane
integer, intent(in) :: orbital
integer, intent(in) :: spin
complex(8), intent(in) :: state(dane%size_tot(1,2) + dane%size_tot(2,2) + dane%size_tot(3,2))
complex(8), intent(out) :: new_state(new_dane%size_tot(1,2) + new_dane%size_tot(2,2) + new_dane%size_tot(3,2))
complex(8) :: temp_table1(dane%n_strings_alpha,dane%n_strings_beta), temp_table2(dane%n_strings_alpha_p1,dane%n_strings_beta_m1), temp_table3(dane%n_strings_alpha_m1,dane%n_strings_beta_p1)
complex(8) :: new_temp_table1(new_dane%n_strings_alpha,new_dane%n_strings_beta), new_temp_table2(new_dane%n_strings_alpha_p1,new_dane%n_strings_beta_m1), new_temp_table3(new_dane%n_strings_alpha_m1,new_dane%n_strings_beta_p1)
integer :: i,j,mu
   new_state(:) = 0
   new_temp_table1(:,:)=0
   new_temp_table2(:,:)=0
   new_temp_table3(:,:)=0
   if (((dane%n_alpha-new_dane%n_alpha .eq. 1) .and. ((dane%n_beta .eq. new_dane%n_beta) .and. (spin .eq. 1)))) then
      continue
   else if (((dane%n_beta-new_dane%n_beta .eq. 1) .and. ((dane%n_alpha .eq. new_dane%n_alpha) .and. (spin .eq. -1)))) then
      continue
   else
      write(*,*) "WRONG NUMBER OF ELECTRONS"
      stop
   end if
   do i=1,dane%n_strings_alpha
      mu = (i-1)*dane%n_strings_beta+1
      do j=1,dane%n_strings_beta
         temp_table1(i,j) = state(mu)
         mu = mu + 1
      end do
   end do
if (dane%relativistic .eqv. .true.) then
   do i=1,dane%n_strings_alpha_p1
      mu = (i-1)*dane%n_strings_beta_m1+1
      do j=1,dane%n_strings_beta_m1
         temp_table2(i,j) = state(mu+dane%size_tot(1,2))
         mu = mu + 1
      end do
   end do
   
   do i=1,dane%n_strings_alpha_m1
      mu = (i-1)*dane%n_strings_beta_p1+1
      do j=1,dane%n_strings_beta_p1
         temp_table3(i,j) = state(mu+dane%size_tot(1,2)+dane%size_tot(2,2))
         mu = mu + 1
      end do
   end do
end if
   if (spin .eq. 1) then
      do i=1,dane%n_strings_alpha
         if (dane%alpha_annihilation_matrix(orbital,i,1) .ne. 0) then
            do j=1,dane%n_strings_beta
               new_temp_table1(dane%alpha_annihilation_matrix(orbital,i,1),j) = dane%alpha_annihilation_matrix(orbital,i,2)*temp_table1(i,j)
            end do
         end if
      end do
if (dane%relativistic .eqv. .true.) then
      do i=1,dane%n_strings_alpha_p1
         if (dane%alpha_annihilation_matrix_p1(orbital,i,1) .ne. 0) then
            do j=1,dane%n_strings_beta_m1
               new_temp_table2(dane%alpha_annihilation_matrix_p1(orbital,i,1),j) = dane%alpha_annihilation_matrix_p1(orbital,i,2)*temp_table2(i,j)
            end do
         end if
      end do
   
      do i=1,dane%n_strings_alpha_m1
         if (new_dane%alpha_annihilation_matrix(orbital,i,1) .ne. 0) then
            do j=1,dane%n_strings_beta_p1
               new_temp_table3(new_dane%alpha_annihilation_matrix(orbital,i,1),j) = new_dane%alpha_annihilation_matrix(orbital,i,2)*temp_table3(i,j)
            end do
         end if
      end do
end if 
   else

      do j=1,dane%n_strings_beta
         if (dane%beta_annihilation_matrix(orbital,j,1) .ne. 0) then
            do i=1,dane%n_strings_alpha
               new_temp_table1(i,dane%beta_annihilation_matrix(orbital,j,1)) = dane%beta_annihilation_matrix(orbital,j,1)*temp_table1(i,j)*(-1)**dane%n_alpha
            end do
         end if
      end do
if (dane%relativistic .eqv. .true.) then
      do j=1,dane%n_strings_beta_m1
         if (new_dane%beta_annihilation_matrix(orbital,j,1) .ne. 0) then
            do i=1,dane%n_strings_alpha_p1
               new_temp_table2(i,new_dane%beta_annihilation_matrix(orbital,j,1)) = new_dane%beta_annihilation_matrix(orbital,j,1)*temp_table2(i,j)*(-1)**(dane%n_alpha+1)
            end do
         end if
      end do

      do j=1,dane%n_strings_beta_p1
         if (dane%beta_annihilation_matrix_p1(orbital,j,1) .ne. 0) then
            do i=1,dane%n_strings_alpha_m1
               new_temp_table3(i,dane%beta_annihilation_matrix_p1(orbital,j,1)) = dane%beta_annihilation_matrix_p1(orbital,j,1)*temp_table3(i,j)*(-1)**(dane%n_alpha-1)
            end do
         end if
      end do
   end if
   end if
do i=1,new_dane%n_strings_alpha
   mu = (i-1)*new_dane%n_strings_beta+1
   do j=1,new_dane%n_strings_beta
      new_state(mu) =  new_temp_table1(i,j)
      mu = mu + 1
   end do
end do
if (dane%relativistic .eqv. .true.) then
do i=1,new_dane%n_strings_alpha_p1
   mu = (i-1)*new_dane%n_strings_beta_m1+1
   do j=1,new_dane%n_strings_beta_m1
      new_state(mu+new_dane%size_tot(1,2)) =  new_temp_table2(i,j)
      mu = mu + 1
   end do
end do

do i=1,new_dane%n_strings_alpha_m1
   mu = (i-1)*new_dane%n_strings_beta_p1+1
   do j=1,new_dane%n_strings_beta_p1
      new_state(mu+new_dane%size_tot(1,2)+new_dane%size_tot(2,2)) =  new_temp_table3(i,j)
      mu = mu + 1
   end do
end do
end if
end subroutine create_hole

subroutine diff_spin_product(new_state_alpha,new_dane_alpha,new_state_beta,new_dane_beta,scalar_product,sign)
type(parameters), intent(in) :: new_dane_alpha, new_dane_beta
complex(8), intent(in) :: new_state_alpha(new_dane_alpha%size_tot(1,2)+new_dane_alpha%size_tot(2,2)+new_dane_alpha%size_tot(3,2)),new_state_beta(new_dane_beta%size_tot(1,2)+new_dane_beta%size_tot(2,2)+new_dane_beta%size_tot(3,2))
integer, intent(in) :: sign
complex(8), intent(out) :: scalar_product
integer :: i
scalar_product = 0
!ALPHA ALWAYS ON THE LEFT
if ((new_dane_alpha%relativistic .eqv. .true.) .and. (new_dane_beta%relativistic .eqv. .true.)) then
if (sign .eq. 1) then
do i=1,new_dane_alpha%size_tot(1,2)
   scalar_product = scalar_product + conjg(new_state_alpha(i))*new_state_beta(i+new_dane_beta%size_tot(1,2))
end do

do i=1,new_dane_alpha%size_tot(3,2)
   scalar_product = scalar_product + conjg(new_state_alpha(i+new_dane_alpha%size_tot(1,2)+new_dane_alpha%size_tot(2,2)))*new_state_beta(i)
end do
else if (sign .eq. -1) then
do i=1,new_dane_alpha%size_tot(1,2)
   scalar_product = scalar_product + conjg(new_state_alpha(i))*new_state_beta(i+new_dane_beta%size_tot(1,2)+new_dane_beta%size_tot(2,2))
end do

do i=1,new_dane_alpha%size_tot(2,2)
   scalar_product = scalar_product + conjg(new_state_alpha(i+new_dane_alpha%size_tot(1,2)))*new_state_beta(i)
end do
end if
end if
end subroutine diff_spin_product


subroutine continued_fraction(params,state,z,krylov_size,sign,fraction)
type(parameters), intent(inout) :: params
complex(8), intent(in) :: state(params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2))
complex(8), intent(in) :: z
integer, intent(in) :: krylov_size, sign
complex(8), intent(out) :: fraction
complex(8) :: a(krylov_size), b(krylov_size)
complex(8) :: state_1(params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)),state_2(params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)),state_3(params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2))
complex(8) :: temp_state(params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2))
integer :: i

state_1 = state
b(1) = 0
do i=1,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)
   b(1) = b(1) + abs(state_1(i))**2
end do
b(1) = sqrt(b(1))


state_1 = state_1/b(1)


call matrix_vector_product(params, state_1, temp_state)

a(1) = 0
do i=1,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)
   a(1) = a(1) + conjg(state_1(i))*temp_state(i)
end do

state_2 = temp_state - a(1)*state_1

b(2) = 0
do i=1,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)
   b(2) = b(2) + abs(state_2(i))**2
end do
b(2) = sqrt(b(2))

state_2 = state_2/b(2)
call matrix_vector_product(params, state_2, temp_state)
a(2) = 0
do i=1,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)
   a(2) = a(2) + conjg(state_2(i))*temp_state(i)
end do
do i=3,krylov_size
   state_3 = temp_state - a(i-1)*state_2 - b(i-1)*state_1
   b(i) = 0
   do j=1,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)
      b(i) = b(i) + abs(state_3(j))**2
   end do
   b(i) = sqrt(b(i))


   state_3 = state_3/b(i)

   call matrix_vector_product(params, state_3, temp_state)

   a(i) = 0
   do j=1,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)
      a(i) = a(i) + conjg(state_3(j))*temp_state(j)
   end do
   state_1 = state_2
   state_2 = state_3
end do

a = sign*a
fraction = 0
do i=krylov_size,1,-1
   fraction = b(i)**2/(z - a(i) - fraction)
end do
end subroutine continued_fraction

subroutine mixed_continued_fraction(params1,state_1,params2,state_2,z,krylov_size,sign,fraction)
type(parameters), intent(inout) :: params1, params2
complex(8), intent(in):: state_1(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)), state_2(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2))
complex(8), intent(in) :: z
integer, intent(in) :: sign
complex(8), intent(out) :: fraction
complex(8) :: state_1_1(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)),state_1_2(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)),state_1_3(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2))
complex(8) :: state_2_1(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)),state_2_2(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)),state_2_3(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2))
complex(8) :: temp_state_1(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2))
complex(8) :: temp_state_2(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2))
complex(8) :: a(krylov_size), b(krylov_size)
integer :: i

state_1_1 = state_1
state_2_1 = state_2
call diff_spin_product(state_1_1,params1,state_2_1,params2,b(1),sign)
b(1) = b(1) + conjg(b(1))
do i=1,params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)
   b(1) = b(1) + conjg(state_1_1(i))*state_1_1(i)
end do

do i=1,params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)
   b(1) = b(1) + conjg(state_2_1(i))*state_2_1(i)
end do
b(1) = sqrt(b(1))
state_1_1 = state_1_1/b(1)
state_2_1 = state_2_1/b(1)

call matrix_vector_product(params1, state_1_1, temp_state_1)
call matrix_vector_product(params2, state_2_1, temp_state_2)

call diff_spin_product(state_1_1,params1,temp_state_2,params2,a(1),sign)
a(1) = a(1) + conjg(a(1))
do i=1,params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)
   a(1) = a(1) + conjg(state_1_1(i))*temp_state_1(i)
end do

do i=1,params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)
   a(1) = a(1) + conjg(state_2_1(i))*temp_state_2(i)
end do


state_1_2 = temp_state_1 - a(1)*state_1_1
state_2_2 = temp_state_2 - a(1)*state_2_1

call diff_spin_product(state_1_2,params1,state_2_2,params2,b(2),sign)
b(2) = b(2) + conjg(b(2))
do i=1,params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)
   b(2) = b(2) + conjg(state_1_2(i))*state_1_2(i)
end do

do i=1,params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)
   b(2) = b(2) + conjg(state_2_2(i))*state_2_2(i)
end do
b(2) = sqrt(b(2))
state_1_2 = state_1_2/b(2)
state_2_2 = state_2_2/b(2)

call matrix_vector_product(params1, state_1_2, temp_state_1)
call matrix_vector_product(params2, state_2_2, temp_state_2)

call diff_spin_product(state_1_2,params1,temp_state_2,params2,a(2),sign)
a(2) = a(2) + conjg(a(2))
do i=1,params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)
   a(2) = a(2) + conjg(state_1_2(i))*temp_state_1(i)
end do

do i=1,params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)
   a(2) = a(2) + conjg(state_2_2(i))*temp_state_2(i)
end do

if (krylov_size .gt. 2) then
   do i=3,krylov_size
      state_1_3 = temp_state_1 - a(i-1)*state_1_2 - b(i-1)*state_1_1
      state_2_3 = temp_state_2 - a(i-1)*state_2_2 - b(i-1)*state_2_1
      call diff_spin_product(state_1_3,params1,state_2_3,params2,b(i),sign)
      b(i) = b(i) + conjg(b(i))
      
      do j=1,params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)
         b(i) = b(i) + conjg(state_1_3(j))*state_1_3(j)
      end do

      do j=1,params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)
         b(i) = b(i) + conjg(state_2_3(j))*state_2_3(j)
      end do
      b(i) = sqrt(b(i))
      state_1_3 = state_1_3/b(i)
      state_2_3 = state_2_3/b(i)

      call matrix_vector_product(params1, state_1_3, temp_state_1)
      call matrix_vector_product(params2, state_2_3, temp_state_2)

      call diff_spin_product(state_1_3,params1,temp_state_2,params2,a(i),sign)
      a(i) = a(i) + conjg(a(i))
      do j=1,params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)
         a(i) = a(i) + conjg(state_1_3(j))*temp_state_1(j)
      end do

      do j=1,params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2)
         a(i) = a(i) + conjg(state_2_3(j))*temp_state_2(j)
      end do

      state_1_1 = state_1_2
      state_2_1 = state_2_2
      state_1_2 = state_1_3
      state_2_2 = state_2_3
   end do
end if
a = a * sign
fraction = 0
do i=krylov_size,1,-1
   fraction = b(i)**2/(z - a(i) - fraction)
end do
end subroutine mixed_continued_fraction

subroutine calc_e_gf(orbital1,orbital2,spin1,spin2,params1,params2,gs_params,ground_state,gs_energy,omega,krylov_size,gf)
type(parameters), intent(inout) :: params1, params2, gs_params
integer, intent(in) :: orbital1,orbital2,spin1,spin2, krylov_size
complex(8), intent(in) :: ground_state(gs_params%size_tot(1,2)+gs_params%size_tot(2,2)+gs_params%size_tot(3,2))
real(8), intent(in) :: gs_energy
complex(8), intent(in) :: omega
complex(8), intent(out) :: gf
complex(8) :: fraction_plus, fraction1, fraction2, fraction_i
complex(8) :: state_1(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)),state_2(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2))
call create_electron(orbital1,ground_state,state_1,gs_params,params1,spin1)
call continued_fraction(params1,state_1,omega+gs_energy,krylov_size,1,fraction1)


if ((orbital1 .eq. orbital2) .and. (spin1 .eq. spin2)) then
   gf = fraction1
else
   call create_electron(orbital2,ground_state,state_2,gs_params,params2,spin2)
   call continued_fraction(params2,state_2,omega+gs_energy,krylov_size,1,fraction2)
   if (spin1 .eq. spin2) then
      call continued_fraction(params1,state_1+state_2,omega+gs_energy,krylov_size,1,fraction_plus)
      call continued_fraction(params1,(0.0d0,-1.0d0)*state_1+state_2,omega+gs_energy,krylov_size,1,fraction_i)
      gf = 0.5*(fraction_plus-(0.0d0,1.0d0)*fraction_i-(1.0d0,-1.0d0)*(fraction1+fraction2))
   else
      call mixed_continued_fraction(params1,state_1,params2,state_2,omega+gs_energy,krylov_size,1,fraction_plus)
      call mixed_continued_fraction(params1,(0.0d0,-1.0d0)*state_1,params2,state_2,omega+gs_energy,krylov_size,1,fraction_i)
      gf = 0.5*(fraction_plus-(0.0d0,1.0d0)*fraction_i-(1.0d0,-1.0d0)*(fraction1+fraction2))
   end if
end if
end subroutine calc_e_gf



subroutine calc_h_gf(orbital1,orbital2,spin1,spin2,params1,params2,gs_params,ground_state,gs_energy,omega,krylov_size,gf)
type(parameters), intent(inout) :: params1, params2, gs_params
integer, intent(in) :: orbital1,orbital2,spin1,spin2, krylov_size
complex(8), intent(in) :: ground_state(gs_params%size_tot(1,2)+gs_params%size_tot(2,2)+gs_params%size_tot(3,2))
real(8), intent(in) :: gs_energy
complex(8), intent(in) :: omega
complex(8), intent(out) :: gf
complex(8) :: fraction_plus, fraction1, fraction2, fraction_i
complex(8) :: state_1(params1%size_tot(1,2)+params1%size_tot(2,2)+params1%size_tot(3,2)),state_2(params2%size_tot(1,2)+params2%size_tot(2,2)+params2%size_tot(3,2))
!INDEX CONVENTION NEEDS TO BE ESTABLISHED
call create_hole(orbital1,ground_state,state_1,gs_params,params1,spin1)
call continued_fraction(params1,state_1,omega-gs_energy,krylov_size,-1,fraction1)

if ((orbital1 .eq. orbital2) .and. (spin1 .eq. spin2)) then
   gf = fraction1
else
   call create_hole(orbital2,ground_state,state_2,gs_params,params2,spin2)
   call continued_fraction(params2,state_2,omega-gs_energy,krylov_size,-1,fraction2)
   if (spin1 .eq. spin2) then
      call continued_fraction(params1,state_1+state_2,omega-gs_energy,krylov_size,-1,fraction_plus)
      call continued_fraction(params1,(0.0d0,1.0d0)*state_1+state_2,omega-gs_energy,krylov_size,-1,fraction_i)
      gf = 0.5*(fraction_plus-(0.0d0,1.0d0)*fraction_i-(1.0d0,-1.0d0)*(fraction1+fraction2))
   else
      call mixed_continued_fraction(params1,state_1,params2,state_2,omega-gs_energy,krylov_size,-1,fraction_plus)
      call mixed_continued_fraction(params1,(0.0d0,1.0d0)*state_1,params2,state_2,omega-gs_energy,krylov_size,-1,fraction_i)
      gf = 0.5*(fraction_plus-(0.0d0,1.0d0)*fraction_i-(1.0d0,-1.0d0)*(fraction1+fraction2))
   end if
end if
end subroutine calc_h_gf









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eigensystem(params,n_pairs,eigenenergies,eigenstates,ierr)
type(parameters),target, intent(in) :: params
integer, intent(in) :: n_pairs
real(8), intent(inout) :: eigenenergies(n_pairs)
complex(8), intent(inout) :: eigenstates(n_pairs,params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2)) 
integer :: N
PetscScalar          :: real_part, imaginary_part
Vec                  :: vector_petsc, im_vector_petsc
PetscScalar, pointer :: tablica_wynikowa(:)
PetscErrorCode, intent(inout) :: ierr
Mat            :: A
EPS            :: eps
integer              :: nconv
complex(8),target :: diagonal(params%size_tot(1,2)+params%size_tot(2,2)+params%size_tot(3,2))

diagonal = params%diag
ptr_dane => params
ptr_przekatna => diagonal
N = params%size_tot(1,2) + params%size_tot(2,2) + params%size_tot(3,2)


call MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, &
                      N, N, PETSC_NULL_INTEGER, A, ierr)

call MatShellSetOperation(A, MATOP_MULT, WrapperMatMult, ierr)
call MatShellSetOperation(A, MATOP_GET_DIAGONAL, WrapperMatGetDiagonal, ierr)

call EPSCreate(PETSC_COMM_WORLD, eps, ierr)
call EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr)

call EPSSetProblemType(eps, EPS_HEP, ierr)

call EPSSetType(eps,EPSJD, ierr)
  
  
call EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr)
call EPSSetDimensions(eps, n_pairs, PETSC_DECIDE, PETSC_DECIDE, ierr)
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
    
    call VecGetArrayRead(vector_petsc, tablica_wynikowa, ierr)
    print *, "Eigenvalue no.", i + 1, "=", real(real_part) + params%nuclear_energy
    if (i .lt. n_pairs) then
    eigenenergies(i+1) = real(real_part) + params%nuclear_energy
    eigenstates(i+1,:) = tablica_wynikowa
    end if
  end do
    

call VecRestoreArrayRead(vector_petsc, tablica_wynikowa, ierr)

call VecDestroy(vector_petsc, ierr)
call VecDestroy(im_vector_petsc, ierr)



call EPSDestroy(eps, ierr)
call MatDestroy(A, ierr)



nullify(ptr_dane)
nullify(ptr_przekatna)
end subroutine eigensystem
end module library_rel_ed

