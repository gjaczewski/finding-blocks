subroutine check_orbital_space_declarations(norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,n_alpha,n_beta)

  implicit none
  integer, intent(in) :: n_alpha, n_beta
  integer, intent(in) :: norb,n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer:: orb,i


  !norb - total number of orbitals
  ! the principles for the occupied and virtual RAS spaces are the same
  !number_RAS_spaces_occ - the number of RAS spaces which are present for occupied orbitals (eg. 3 spaces with core, singly, and doubly occupied orbitals out of total number of occupied orbitals)
  !dimension_RAS_space_occ(number_RAS_spaces_occ) - dimension/number of orbitals in each of the occupied RAS space (eg. dimension_RAS_space_occ(1)=number of orbitals in the first RAS space
  !dimension_active_space - number of orbitals in the active space 
  
  !check if number of all the orbitals is equal to norb

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

  if (orb .ne. norb) then
     write(*,*)"SOMETHING WENT WRONG,number of orbitals declared in RAS spaces or active space is wrong (number of orbitals, number or orbitals in RAS+ACTIVE)", norb, orb
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


subroutine count_distributions(n_s, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions)
implicit none
integer, intent(in) :: n_s, n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
integer ::i,j
integer, intent(in) :: n_combinations,n_spaces
integer, intent(in) :: all_combinations(n_combinations,n_spaces)
integer :: a,v,r, counter, sum, checker
integer, intent(out) :: n_distributions
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


subroutine count_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta)
implicit none
logical, intent(in) :: relativistic
integer, intent(in) :: n_alpha,n_beta
integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt), RAS_space_virt(n_RAS_spaces_virt)
integer, intent(in) :: n_combinations,n_spaces
integer, intent(in) :: all_combinations(n_combinations,n_spaces)
integer, intent(out) :: n_distributions_alpha,n_distributions_beta
integer, intent(out) :: n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1
integer, intent(out) :: NRD_spin_alpha,NRD_spin_beta

call count_distributions(n_alpha, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha)
call count_distributions(n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta)
  if (n_alpha .le. 0) then
    n_distributions_alpha_m1 = 0

  else
    call count_distributions(n_alpha-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_m1)
  end if
    if (n_beta .le. 0) then
    n_distributions_beta_m1 = 0

    else
    call count_distributions(n_beta-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_m1)

  end if
if (relativistic .eqv. .true.) then
  

      call count_distributions(n_alpha+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_p1)
      call count_distributions(n_beta+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_p1)
      
else 
    n_distributions_alpha_p1 = 0    
    n_distributions_beta_p1 = 0
end if 


NRD_spin_alpha = n_distributions_alpha + n_distributions_alpha_m1 + n_distributions_alpha_p1
NRD_spin_beta= n_distributions_beta + n_distributions_beta_m1 + n_distributions_beta_p1
end subroutine count_spin_distributions


subroutine find_spin_distributions(relativistic,n_alpha,n_beta, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,RAS_space_virt,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha,n_distributions_beta,n_distributions_alpha_p1,n_distributions_beta_p1,n_distributions_alpha_m1,n_distributions_beta_m1,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,NRDa,NRDb)
implicit none
logical, intent(in) :: relativistic
integer, intent(in) :: n_alpha,n_beta
integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt, active_space
integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ), excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt), RAS_space_virt(n_RAS_spaces_virt)
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
if (n_alpha .gt. 0) then
call find_distributions(n_alpha-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_m1,RAS_el_array_alpha(n_distributions_alpha+n_distributions_alpha_p1+1:n_distributions_alpha_m1+n_distributions_alpha+n_distributions_alpha_p1,:))
end if
if (n_beta .gt. 0) then
call find_distributions(n_beta-1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_m1,RAS_el_array_beta(n_distributions_beta+n_distributions_beta_p1+1:n_distributions_beta_m1+n_distributions_beta+n_distributions_beta_p1,:))
end if
if (relativistic .eqv. .true.) then
call find_distributions(n_alpha+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_alpha_p1,RAS_el_array_alpha(n_distributions_alpha+1:n_distributions_alpha+n_distributions_alpha_p1,:))
call find_distributions(n_beta+1, n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,active_space,excit_array,n_combinations,n_spaces,all_combinations,n_distributions_beta_p1,RAS_el_array_beta(n_distributions_beta+1:n_distributions_beta+n_distributions_beta_p1,:))
end if

end subroutine find_spin_distributions


subroutine calculate_space_size(NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,size_tot,verbose)
  implicit none
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in)::NRDa(3,2),NRDb(3,2)
  integer, intent(in):: NRD_spin_alpha,NRD_spin_beta
  integer, intent(in):: RAS_el_array_alpha(NRD_spin_alpha,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_el_array_beta(NRD_spin_beta,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(out)::sizea(3,2),sizeb(3,2),size_tot(3,2) 
  integer, intent(in):: verbose
  

  integer:: vector_length
  
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

  if (verbose==2) then
     write(*,*)"SIZES na, nb", size_tot(1,1),size_tot(1,2),sizea(1,2),sizeb(1,2)
     write(*,*)"SIZES na+1,nb-1", size_tot(2,1),size_tot(2,2),sizea(2,2),sizeb(3,2)
     write(*,*)"SIZES na-1,nb+1", size_tot(3,1),size_tot(3,2),sizea(3,2),sizeb(2,2)
     
     vector_length=size_tot(3,1)+size_tot(3,2)-1
     write(*,*)"vector_length",vector_length
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
  real(8) :: nCr_dp
  !write(*,*)"range1,range2",range1,range2
  vec_length=0
  do i=range1,range2
     a=n_RAS_spaces_occ+1 
     siz=1
     do r=1,n_RAS_spaces_occ
        siz=siz*int(nCr_dp(RAS_space_occ(r),RAS_el_array_spin(i,r)))
     end do
     do v=1,n_RAS_spaces_virt
        siz=siz*int(nCr_dp(RAS_space_virt(v),RAS_el_array_spin(i,a+v)))
     end do
     siz=siz*int(nCr_dp(active_space,RAS_el_array_spin(i,a)))
     vec_length=vec_length+siz
  end do
  !write(*,*)"vec_length",vec_length
  call flush(6)
end subroutine calc_size_block


function nCr_dp(n, r) result(c)
    implicit none
    integer, intent(in) :: n, r
    real(8) :: c
    integer :: i, k

    if (r < 0 .or. r > n) then
        c = 0.0d0
        return
    end if

    k = min(r, n - r)
    c = 1.0d0

    do i = 1, k
        c = c * dble(n - k + i) / dble(i)
    end do
  end function nCr_dp


  subroutine fill_spin_strings(n_s,NRD_spin,RAS_el_array_spin,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,range1,range2,sizes,str_s)
  
  implicit none
  integer, intent(in):: n_s
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in):: range1,range2
  integer, intent(in):: NRD_spin
  integer, intent(in):: RAS_el_array_spin(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(in)::sizes
  integer, intent(inout)::str_s(sizes,n_s)
  integer :: space_sizes(n_RAS_spaces_occ+1+n_RAS_spaces_virt)
  integer, allocatable :: str_temp(:,:)
  integer, allocatable :: another_str_temp(:,:)
  integer :: siz_space, tot_siz_space
  integer :: n_el_temp
  integer :: orbital_index
  integer :: n_str_temp(n_RAS_spaces_occ+1+n_RAS_spaces_virt)
  integer :: j
  integer :: a, i , r ,v
  real(8) :: nCr_dp

  a=n_RAS_spaces_occ+1
  tot_siz_space = 0
  do i=range1,range2
      siz_space = 1
     orbital_index = 1
     j = 1
     n_el_temp = 0

   do r=1,n_RAS_spaces_occ
      space_sizes(j) = int(nCr_dp(RAS_space_occ(r),RAS_el_array_spin(i,r)))
      siz_space=siz_space*space_sizes(j)
      n_str_temp(j) = siz_space
      j = j + 1

   end do

   space_sizes(j) = int(nCr_dp(active_space,RAS_el_array_spin(i,a)))
   siz_space=siz_space*space_sizes(j)
   n_str_temp(j) = siz_space
   j = j + 1

   do v=1,n_RAS_spaces_virt
   space_sizes(j) = int(nCr_dp(RAS_space_virt(v),RAS_el_array_spin(i,a+v)))

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
   call sort_states(str_s,sizes,n_s)
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


subroutine sort_states(spin_strings,n_spin_strings,n_spin)
implicit none
integer, intent(in) :: n_spin_strings,n_spin
integer, intent(inout) :: spin_strings(n_spin_strings,n_spin)
real :: temp_array(n_spin_strings), temp_energy
integer :: i,j
integer :: temp_state(n_spin)
if (n_spin_strings .ge. 2) then
temp_array(:) = 0
do i=1,n_spin_strings
   do j=1,n_spin
      temp_array(i) = temp_array(i) - 1/(real(spin_strings(i,j)))**2  
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


subroutine ISEQUAL(vector1,vector2,length,indicator)
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
end subroutine ISEQUAL





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


subroutine fill_annihilation_results(n_spin_strings_m1,n_spin_strings,n_spin,spin_strings_m1,spin_strings,norb,spin_annihilation_matrix)
implicit none
integer, intent(in) :: n_spin_strings_m1,n_spin_strings,n_spin,norb
integer, intent(in) :: spin_strings_m1(n_spin_strings_m1,n_spin-1), spin_strings(n_spin_strings,n_spin)
integer, intent(out) :: spin_annihilation_matrix(norb,n_spin_strings,2)
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
         call ISEQUAL(temp_string,spin_strings_m1(k,:),n_spin-1,indicator)
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





!***********************THIS WHOLE PART MIGHT BE OPTIMISED USING SMALL MATRICE**********************


subroutine fill_annihilation_creation_matrix(n_spin,n_spin_strings,spin_strings,norb,annihilation_creation_matrix)
implicit none
integer, intent(in) ::  n_spin, n_spin_strings, norb
integer, intent(in) :: spin_strings(n_spin_strings,n_spin)
integer, intent(inout) :: annihilation_creation_matrix(norb,norb,n_spin_strings,2)
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
            call ISEQUAL(temp_string1,temp_string2,n_spin-1,indicator)
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


subroutine fill_nr_single_spin_hamiltonian(spin_strings,n_spin_strings,n_spin,norb,hopping,interaction,spin_annihilation_creation_matrix,spin_nr_hamiltonian)
integer, intent(in) :: n_spin_strings,n_spin,norb
integer, intent(in) :: spin_strings(n_spin_strings,n_spin), spin_annihilation_creation_matrix(norb,norb,n_spin_strings,2)
complex, intent(in) :: hopping(norb,norb), interaction(norb,norb,norb,norb)
complex, intent(out) :: spin_nr_hamiltonian(n_spin_strings,n_spin_strings)
integer :: j,p,q,r,s, temp_state1,temp_state2,temp_sign1,temp_sign2

do j=1,n_spin_strings
   do q=1,n_spin
      do p=1,norb
         temp_state1 = spin_annihilation_creation_matrix(p,spin_strings(j,q),j,1)
         temp_sign1 = spin_annihilation_creation_matrix(p,spin_strings(j,q),j,2)
         if (temp_state1 .ne. 0) then
            spin_nr_hamiltonian(temp_state1,j) = spin_nr_hamiltonian(temp_state1,j) + temp_sign1*hopping(p,spin_strings(j,q))
            do r=1,norb
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
end subroutine fill_nr_single_spin_hamiltonian


subroutine nr_matrix_vector_product(alpha_hamiltonian,strings_alpha,n_strings_alpha,n_alpha,beta_hamiltonian,strings_beta,n_strings_beta,n_beta,interaction_mix,norb,alpha_annihilation_creation_matrix,beta_annihilation_creation_matrix,vector,vector_new)
integer, intent (in) :: n_strings_alpha,n_strings_beta,norb, n_alpha, n_beta
integer, intent (in) :: alpha_annihilation_creation_matrix(norb,norb,n_strings_alpha,2), beta_annihilation_creation_matrix(norb,norb,n_strings_beta,2), strings_alpha(n_strings_alpha,n_alpha), strings_beta(n_strings_beta,n_beta)
complex, intent (in) :: alpha_hamiltonian(n_strings_alpha,n_strings_alpha), beta_hamiltonian(n_strings_beta,n_strings_beta), interaction_mix(norb,norb,norb,norb), vector(n_strings_alpha*n_strings_beta)
complex, intent (inout) :: vector_new(n_strings_alpha*n_strings_beta)
complex :: temp_table(n_strings_alpha,n_strings_beta), new_table(n_strings_alpha,n_strings_beta)
integer :: mu, i,j,k,l, p,q,r,s,temp_state1,temp_state2,temp_sign1,temp_sign2
new_table(:,:) = 0
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
      end do
      do l=1,n_strings_beta
         new_table(i,j) = new_table(i,j) + temp_table(i,l)*beta_hamiltonian(j,l)
      end do
      do p=1,n_alpha
         do r=1,norb
            temp_state1 = alpha_annihilation_creation_matrix(r,strings_alpha(i,p),i,1)
            if (temp_state1 .ne. 0) then
               temp_sign1 = alpha_annihilation_creation_matrix(r,strings_alpha(i,p),i,2)
               do q=1,n_beta
                  do s=1,norb
                     temp_state2 = beta_annihilation_creation_matrix(s,strings_beta(j,q),j,1)
                     if (temp_state2 .ne. 0) then
                        temp_sign2 = temp_sign1*beta_annihilation_creation_matrix(s,strings_beta(j,q),j,2)
                        new_table(i,j) = new_table(i,j) + temp_sign2*temp_table(temp_state1,temp_state2)*interaction_mix(strings_alpha(i,p),strings_beta(j,q),r,s)
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


subroutine rel_matrix_vector_product(alpha_hamiltonian,alpha_hamiltonian_p1,alpha_hamiltonian_m1,strings_alpha,n_strings_alpha,strings_alpha_p1,n_strings_alpha_p1,strings_alpha_m1,n_strings_alpha_m1,n_alpha,beta_hamiltonian,beta_hamiltonian_p1,beta_hamiltonian_m1,strings_beta,n_strings_beta,strings_beta_p1,n_strings_beta_p1,strings_beta_m1,n_strings_beta_m1,n_beta,interaction_mix,hso,norb,alpha_annihilation_creation_matrix,alpha_annihilation_creation_matrix_p1,alpha_annihilation_creation_matrix_m1,beta_annihilation_creation_matrix,beta_annihilation_creation_matrix_p1,beta_annihilation_creation_matrix_m1,alpha_annihilation_matrix,alpha_annihilation_matrix_p1,beta_annihilation_matrix,beta_annihilation_matrix_p1,size_tot,vector,vector_new)
integer, intent (in) :: n_strings_alpha,n_strings_beta,norb, n_alpha, n_beta, n_strings_alpha_p1,n_strings_beta_p1,n_strings_alpha_m1,n_strings_beta_m1
integer, intent (in) :: alpha_annihilation_creation_matrix(norb,norb,n_strings_alpha,2), beta_annihilation_creation_matrix(norb,norb,n_strings_beta,2), strings_alpha(n_strings_alpha,n_alpha), strings_beta(n_strings_beta,n_beta)
integer, intent (in) :: alpha_annihilation_creation_matrix_p1(norb,norb,n_strings_alpha_p1,2), beta_annihilation_creation_matrix_p1(norb,norb,n_strings_beta_p1,2), strings_alpha_p1(n_strings_alpha_p1,n_alpha+1), strings_beta_p1(n_strings_beta_p1,n_beta+1)
integer, intent (in) :: alpha_annihilation_creation_matrix_m1(norb,norb,n_strings_alpha_m1,2), beta_annihilation_creation_matrix_m1(norb,norb,n_strings_beta_m1,2), strings_alpha_m1(n_strings_alpha_m1,n_alpha-1), strings_beta_m1(n_strings_beta_m1,n_beta-1)
integer, intent (in) :: alpha_annihilation_matrix(norb,n_strings_alpha,2),beta_annihilation_matrix(norb,n_strings_beta,2),alpha_annihilation_matrix_p1(norb,n_strings_alpha_p1,2),beta_annihilation_matrix_p1(norb,n_strings_beta_p1,2)
integer, intent (in) :: size_tot(3,2)
complex, intent (in) :: alpha_hamiltonian(n_strings_alpha,n_strings_alpha), beta_hamiltonian(n_strings_beta,n_strings_beta), interaction_mix(norb,norb,norb,norb),hso(norb,norb),vector(n_strings_alpha*n_strings_beta+n_strings_alpha_p1*n_strings_beta_m1+n_strings_alpha_m1*n_strings_beta_p1)
complex, intent (in) :: alpha_hamiltonian_p1(n_strings_alpha_p1,n_strings_alpha_p1), beta_hamiltonian_p1(n_strings_beta_p1,n_strings_beta_p1)
complex, intent (in) :: alpha_hamiltonian_m1(n_strings_alpha_m1,n_strings_alpha_m1), beta_hamiltonian_m1(n_strings_beta_m1,n_strings_beta_m1)
complex, intent (inout) :: vector_new(n_strings_alpha*n_strings_beta+n_strings_alpha_p1*n_strings_beta_m1+n_strings_alpha_m1*n_strings_beta_p1)
complex ::  vector_new_1(n_strings_alpha*n_strings_beta),vector_new_2(n_strings_alpha_p1*n_strings_beta_m1),vector_new_3(n_strings_alpha_m1*n_strings_beta_p1)
complex :: temp_table1(n_strings_alpha,n_strings_beta), new_table1(n_strings_alpha,n_strings_beta),temp_table2(n_strings_alpha_p1,n_strings_beta_m1), new_table2(n_strings_alpha_p1,n_strings_beta_m1),temp_table3(n_strings_alpha_m1,n_strings_beta_p1), new_table3(n_strings_alpha_m1,n_strings_beta_p1)
integer :: checker, i,j,k,l,p,q,mu, temp_state1, temp_state2, temp_sign1, temp_sign2
checker = 0
if (n_alpha-1 .eq. 0) then
checker = checker + 1
end if
if (n_beta-1 .eq. 0) then
checker = checker + 2
end if
!DO WDROZENIA
call nr_matrix_vector_product(alpha_hamiltonian,strings_alpha,n_strings_alpha,n_alpha,beta_hamiltonian,strings_beta,n_strings_beta,n_beta,interaction_mix,norb,alpha_annihilation_creation_matrix,beta_annihilation_creation_matrix,vector(size_tot(1,1):size_tot(1,2)),vector_new_1)
call nr_matrix_vector_product(alpha_hamiltonian_p1,strings_alpha_p1,n_strings_alpha_p1,n_alpha+1,beta_hamiltonian_m1,strings_beta_m1,n_strings_beta_m1,n_beta-1,interaction_mix,norb,alpha_annihilation_creation_matrix_p1,beta_annihilation_creation_matrix_m1,vector(size_tot(2,1):size_tot(2,2)),vector_new_2)
call nr_matrix_vector_product(alpha_hamiltonian_m1,strings_alpha_m1,n_strings_alpha_m1,n_alpha-1,beta_hamiltonian_p1,strings_beta_p1,n_strings_beta_p1,n_beta+1,interaction_mix,norb,alpha_annihilation_creation_matrix_m1,beta_annihilation_creation_matrix_p1,vector(size_tot(3,1):size_tot(3,2)),vector_new_3)

new_table1(:,:) = 0
new_table2(:,:) = 0
new_table3(:,:) = 0
do i=1,n_strings_alpha
   mu = (i-1)*n_strings_beta+1
   do j=1,n_strings_beta
      temp_table1(i,j) = vector(mu)
      mu = mu + 1
   end do
end do

do i=1,n_strings_alpha_p1
   mu = (i-1)*n_strings_beta_m1+size_tot(2,1)
   do j=1,n_strings_beta_m1
      temp_table2(i,j) = vector(mu)
      mu = mu + 1
   end do
end do

do i=1,n_strings_alpha_m1
   mu = (i-1)*n_strings_beta_p1+size_tot(3,1)
   do j=1,n_strings_beta_p1
      temp_table3(i,j) = vector(mu)
      mu = mu + 1
   end do
end do

do j=1,n_strings_beta
      do p=1,n_beta
         temp_state1 = beta_annihilation_matrix(strings_beta(j,p),j,1)
         if (temp_state1 .ne. 0) then
            temp_sign1 = beta_annihilation_matrix(strings_beta(j,p),j,2)
            do k=1,n_strings_alpha_p1
               do q=1,n_alpha+1
                  temp_state2 = alpha_annihilation_matrix_p1(strings_alpha_p1(k,q),k,1)
                  if (temp_state2 .ne. 0) then
                     temp_sign2 = temp_sign1 * alpha_annihilation_matrix_p1(strings_alpha_p1(k,q),k,2)
                     new_table1(temp_state2,j) = new_table1(temp_state2,j) +  temp_sign2*hso(strings_beta(j,p),strings_alpha_p1(k,q))*temp_table2(k,temp_state1)*(-1)**n_alpha
                  end if
               end do
            end do
         end if
      end do
end do

do i=1,n_strings_alpha
   do q=1,n_alpha
      temp_state1 = alpha_annihilation_matrix(strings_alpha(i,q),i,1) 
      if (temp_state1 .ne. 0) then
         temp_sign1 = alpha_annihilation_matrix(strings_alpha(i,q),i,2)
         do l=1,n_strings_beta_p1
            do p=1,n_beta+1
               temp_state2 = beta_annihilation_matrix_p1(strings_beta_p1(l,p),l,1)
               if (temp_state2 .ne. 0) then
                  temp_sign2 = temp_sign1*beta_annihilation_matrix_p1(strings_beta_p1(l,p),l,2)
                  new_table1(i,temp_state2) = new_table1(i,temp_state2) + temp_sign2*hso(strings_alpha(i,q),strings_beta_p1(l,p))*temp_table3(temp_state1,l)*(-1)**(n_alpha-1)
               end if
            end do
         end do 
      end if     
   end do
end do

do i=1,n_strings_alpha_p1
   do q=1,n_alpha+1
      temp_state1 = alpha_annihilation_matrix_p1(strings_alpha_p1(i,q),i,1)
      if (temp_state1 .ne. 0) then
        temp_sign1 = alpha_annihilation_matrix_p1(strings_alpha_p1(i,q),i,2)
        do l=1,n_strings_beta
         do p=1,n_beta
            temp_state2 = beta_annihilation_matrix(strings_beta(l,p),l,1)
            if (temp_state2 .ne. 0) then
               temp_sign2 = temp_sign1*beta_annihilation_matrix(strings_beta(l,p),l,2)
               new_table2(i,temp_state2) = new_table2(i,temp_state2) + temp_sign2*hso(strings_alpha_p1(i,q),strings_beta(l,p))*temp_table1(temp_state1,l)*(-1)**n_alpha 
            end if
         end do
        end do        
      end if
   end do
end do

do j=1,n_strings_beta_p1
   do p=1,n_beta+1
      temp_state1 = beta_annihilation_matrix_p1(strings_beta_p1(j,p),j,1)
      if (temp_state1 .ne. 0) then
         temp_sign1 = beta_annihilation_matrix_p1(strings_beta_p1(j,p),j,2)
         do k=1,n_strings_alpha
            do q=1,n_alpha
               temp_state2 = alpha_annihilation_matrix(strings_alpha(k,q),k,1)
               if (temp_state2 .ne. 0) then
                  temp_sign2 = temp_sign1*alpha_annihilation_matrix(strings_alpha(k,q),k,2)
                  new_table3(temp_state2,j) = new_table3(temp_state2,j) + temp_sign2*hso(strings_beta_p1(j,p),strings_alpha(k,q))*temp_table1(k,temp_state1)*(-1)**(n_alpha-1)
               end if
            end do
         end do
      end if
   end do
end do

do i=1,n_strings_alpha
   mu = (i-1)*n_strings_beta+1
   do j=1,n_strings_beta
      vector_new_1(mu) = vector_new_1(mu) + new_table1(i,j)
      mu = mu + 1
   end do
end do

do i=1,n_strings_alpha_p1
   mu = (i-1)*n_strings_beta_m1+1
   do j=1,n_strings_beta_m1
      vector_new_2(mu) = vector_new_2(mu) + new_table2(i,j)
      mu = mu + 1
   end do
end do

do i=1,n_strings_alpha_m1
   mu = (i-1)*n_strings_beta_p1+1
   do j=1,n_strings_beta_p1
      vector_new_3(mu) = vector_new_3(mu) + new_table3(i,j)
      mu = mu + 1
   end do
end do
vector_new(1:size_tot(1,2)) = vector_new_1(:)
vector_new(size_tot(2,1):size_tot(3,1)-1) = vector_new_2(:)
vector_new(size_tot(3,1):size_tot(1,2)+size_tot(2,2)+size_tot(3,2)) = vector_new_3(:)

end subroutine rel_matrix_vector_product


subroutine generate_diagonal_elements(alpha_hamiltonian,strings_alpha,n_strings_alpha,n_alpha,beta_hamiltonian,strings_beta,n_strings_beta,n_beta,interaction,norb,diagonal)
integer, intent (in) :: n_strings_alpha,n_strings_beta,norb, n_alpha, n_beta
integer, intent (in) :: strings_alpha(n_strings_alpha,n_alpha), strings_beta(n_strings_beta,n_beta)
complex, intent (in) :: alpha_hamiltonian(n_strings_alpha,n_strings_beta), beta_hamiltonian(n_strings_beta,n_strings_beta), interaction(norb,norb,norb,norb)
real, intent(out) :: diagonal(n_strings_alpha*n_strings_beta)
integer :: mu,i,j,p,q
complex :: temp_diagonal(n_strings_alpha*n_strings_beta)
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
!**********UNUSED, BUT POTENTIALLY USEFUL ROUTINES**************


subroutine sort(a, n)
    implicit none
    integer, intent(in)    :: n
    integer, intent(inout) :: a(n)
    integer :: i, j, temp

    do i = 1, n-1
        do j = 1, n-i
            if (a(j) > a(j+1)) then
                temp   = a(j)
                a(j)   = a(j+1)
                a(j+1) = temp
            end if
        end do
    end do
end subroutine sort


subroutine find_number(strings_alpha,strings_beta,string_alpha,string_beta,n_strings_alpha,n_strings_beta,n_alpha,n_beta,string_number)
implicit none
integer, intent(in) :: n_strings_alpha,n_strings_beta,n_alpha,n_beta
integer, intent(in) :: strings_alpha(n_strings_alpha,n_alpha), strings_beta(n_strings_beta,n_beta)
integer, intent(in) :: string_alpha(n_alpha), string_beta(n_beta)
integer, intent(inout) :: string_number
integer :: i, j, k
logical :: indicator
do i=1,n_strings_alpha
   call ISEQUAL(string_alpha,strings_alpha(i,:),n_alpha,indicator)
   if (indicator .eqv. .true.) then
      j = i
      exit
   end if
end do

do i=1,n_strings_beta
   call ISEQUAL(string_beta,strings_beta(i,:),n_beta,indicator)
   if (indicator .eqv. .true.) then
      k = i
      exit
   end if
end do
string_number = (j-1)*n_strings_beta+k
end subroutine find_number


subroutine creation(orbital,string_spin,n_spin,new_string_spin,sign)
implicit none
integer, intent(in) :: orbital,n_spin
integer, intent(in) :: string_spin(n_spin)
integer, intent(out) :: new_string_spin(n_spin+1)
integer, intent(out) :: sign
integer :: i,j
logical :: check
check = .true.
do i=1,n_spin
   if (string_spin(i) .eq. orbital) then
      check = .false.
      exit
   end if 
end do
if (check .eqv. .true.) then
   do i=1,n_spin
   new_string_spin(i+1) = string_spin(i)
   end do
   new_string_spin(1) = orbital
   call sort(new_string_spin,n_spin+1)
   do j=1,n_spin+1
      if (new_string_spin(j) .eq. orbital) then
      sign = (-1)**(j-1)
      exit
      end if
   end do
else
   sign = 1
   new_string_spin(:) = 0
end if
end subroutine creation


subroutine fill_creation_results(n_spin_strings_p1,n_spin_strings,n_spin,spin_strings_p1,spin_strings,norb,spin_creation_matrix)
implicit none 
integer, intent(in) :: n_spin_strings_p1,n_spin_strings,n_spin,norb
integer, intent(in) :: spin_strings_p1(n_spin_strings_p1,n_spin+1), spin_strings(n_spin_strings,n_spin)
integer, intent(out) :: spin_creation_matrix(norb,n_spin_strings,2)
integer :: temp_string(n_spin+1)
integer :: sign, i, j, k
logical :: indicator
spin_creation_matrix(:,:,1) = 0
spin_creation_matrix(:,:,2) = 1
do i=1,n_spin_strings
   do j=1,norb
      call creation(j,spin_strings(i,:),n_spin,temp_string,sign)
      do k=1,n_spin_strings_p1
         call ISEQUAL(temp_string,spin_strings_p1(k,:),n_spin+1,indicator)
         if (indicator .eqv. .true.) then
            spin_creation_matrix(j,i,1) = k
            spin_creation_matrix(j,i,2) = sign
            exit
         end if
      end do
   end do
end do
end subroutine fill_creation_results








