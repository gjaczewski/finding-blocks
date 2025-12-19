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




subroutine fill_RAS_electron_array(relativistic,norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,n_alpha,n_beta,NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,NRDa,NRDb)
  
  implicit none
  logical, intent(in) :: relativistic
  integer, intent(in) :: n_alpha, n_beta
  integer, intent(in) :: norb,n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in) :: excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
  integer, intent(inout):: NRD_spin_alpha,NRD_spin_beta
  integer, intent(in) :: NRD_spin_tmp
  integer, intent(inout):: RAS_el_array_alpha(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_el_array_beta(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(out):: NRDa(3,2),NRDb(3,2)
  integer:: orb,n_spaces
  

  integer:: i,j,k,l
  
  !excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt) - an array describing the excitation level possible in each of the RAS spaces
  k=1
  do i=1,n_RAS_spaces_occ
     if (excit_array(i)>RAS_space_occ(i)) then
        write(*,*)"Excitation level of RAS_space_occ",i,"greater than the size of the space", excit_array(i),RAS_space_occ(i)
        call flush(6)
        STOP
     end if
     k=k+1
  end do
  
  do i=1, n_RAS_spaces_virt
     if (excit_array(k)>RAS_space_virt(i)) then
        write(*,*)"Excitation level of RAS_space_virt",i,"greater than the size of the space", excit_array(k),RAS_space_virt(i)
        call flush(6)
        STOP
     end if
     k=k+1
  end do

  !NRDa(1:3,1:2) - an array of pointers to the RAS_el_array_alpha
  !NRDa(1,:) - for n_alpha electrons
  !NRDa(2,:) - for n_alpha+1 electrons
  !NRDa(3,:) - for n_alpha-1 electrons 
  !NRDa(:,1) - the beginning of such sector
  !NRDa(:,2) - the size (number of elements) of/in such sector
  

  
  ! for alpha electrons
  NRD_spin_alpha=0
  NRDa(1,1)=1
  call fill_RAS_spaces_per_spin(n_alpha,NRD_spin_alpha,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_alpha)
  NRDa(1,2)=NRD_spin_alpha
  
  if (relativistic .eqv. .true.) then
     
     ! for alpha+1 electrons
     if (n_alpha+1 .le.norb) then
        NRDa(2,1)=NRD_spin_alpha+1
        call fill_RAS_spaces_per_spin(n_alpha+1,NRD_spin_alpha,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_alpha)
        NRDa(2,2)=NRD_spin_alpha
     else
        NRDa(2,1)=0
        NRDa(2,2)=0
     end if
     
     ! for alpha-1 electrons
     if (n_alpha-1 .ge.0) then
        NRDa(3,1)=NRD_spin_alpha+1
        call fill_RAS_spaces_per_spin(n_alpha-1,NRD_spin_alpha,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_alpha)
        NRDa(3,2)=NRD_spin_alpha
     else
        NRDa(3,1)=0
        NRDa(3,2)=0
     end if
     
  end if

    
  ! for beta electrons
  NRD_spin_beta=0
  NRDb(1,1)=1
  call fill_RAS_spaces_per_spin(n_beta,NRD_spin_beta,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_beta)
  NRDb(1,2)=NRD_spin_beta

  if (relativistic .eqv. .true.) then

    
     ! for beta+1 electrons
     if (n_beta+1 .le.norb) then
        NRDb(2,1)=NRD_spin_beta+1
        call fill_RAS_spaces_per_spin(n_beta+1,NRD_spin_beta,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_beta)
        NRDb(2,2)=NRD_spin_beta
     else
        NRDb(2,1)=0
        NRDb(2,2)=0
     end if



     ! for beta-1 electrons
     if (n_beta-1 .ge.0) then
        NRDb(3,1)=NRD_spin_beta+1
        call fill_RAS_spaces_per_spin(n_beta-1,NRD_spin_beta,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_beta)
        NRDb(3,2)=NRD_spin_beta
     else
        NRDb(3,1)=0
        NRDb(3,2)=0
     end if


  end if

write(*,*)"Number of distributions alpha",NRD_spin_alpha,"beta",NRD_spin_beta
call flush(6)
end subroutine fill_RAS_electron_array
  
  
subroutine fill_RAS_spaces_per_spin(n_spin,NRD_good,NRD_spin_tmp,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_dis)
  implicit none
  integer, intent(in) :: n_spin
  integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt),excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
  integer, intent(in) :: NRD_spin_tmp
  integer, intent(inout):: NRD_good
  integer, intent(inout):: RAS_dis(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer:: verbose
  

  
  integer, dimension(:,:), allocatable::occ_RAS_el,virt_RAS_el,RAS_dis_temp
  integer, dimension(:), allocatable::counter
  integer:: n_el_temp,active_el,NRD_spin
  integer:: i,j,k,l,m,n,p,q
  logical:: count_possibilities,fill
  integer, dimension(:), allocatable::scratch
  !NRD_spin - number of possible electron distributions for a given alpha or beta spin
  verbose=0

  allocate(scratch(10000))
  
  allocate(counter(n_RAS_spaces_occ+n_RAS_spaces_virt))
  count_possibilities=.true.
  counter(:)=0
  call  count_el_nums_in_ras_spaces(count_possibilities,fill,counter,scratch(1),scratch(1),MAXVAL(counter),n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,excit_array)
  count_possibilities=.false.
  
  allocate(occ_RAS_el(n_RAS_spaces_occ,MAXVAL(counter)))
  allocate(virt_RAS_el(n_RAS_spaces_virt,MAXVAL(counter)))
  write(*,*)"MAXVAL(counter)",MAXVAL(counter)

  write(*,*)"counter",counter(:)
  write(*,*)"n_RAS_spaces_occ,n_RAS_spaces_virt",n_RAS_spaces_occ,n_RAS_spaces_virt
  write(*,*)"RAS_space_occ",RAS_space_occ(:)
   write(*,*)"excit_array",excit_array(:)
  
  fill=.true.
  call  count_el_nums_in_ras_spaces(count_possibilities,fill,counter,occ_RAS_el,virt_RAS_el,MAXVAL(counter),n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,excit_array)
  fill=.false.

    do i=1,n_RAS_spaces_virt
       do j=1,counter(n_RAS_spaces_occ+i)
          write(*,*)"virt_RAS_el(p,l)",virt_RAS_el(i,j)
       end do
    end do


  
  NRD_spin=1
  do i=1,n_RAS_spaces_occ 
     NRD_spin=NRD_spin*counter(i)
  end do
  do k=1,n_RAS_spaces_virt
     NRD_spin=NRD_spin*counter(n_RAS_spaces_occ+k)
  end do
 
  allocate(RAS_dis_temp(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1))
  
  call fill_distribution_arrays(n_spin,NRD_spin_tmp,NRD_spin,NRD_good,RAS_dis_temp,RAS_dis,n_RAS_spaces_occ,n_RAS_spaces_virt,occ_RAS_el,virt_RAS_el,active_space,counter,verbose)

  deallocate(scratch)
  deallocate(counter)
  deallocate(occ_RAS_el)
  deallocate(virt_RAS_el)
  deallocate(RAS_dis_temp)
end subroutine fill_RAS_spaces_per_spin
  
subroutine fill_distribution_arrays(n_spin,NRD_spin_tmp,NRD_spin,NRD_good,RAS_dis_temp,RAS_dis,n_RAS_spaces_occ,n_RAS_spaces_virt,occ_RAS_el,virt_RAS_el,active_space,counter,verbose)
  implicit none
  integer, intent(in)::n_spin,NRD_spin,NRD_spin_tmp
  integer, intent(inout)::NRD_good
  integer, intent(in)::n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
    integer, intent(inout)::RAS_dis_temp(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_dis(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
    integer::counter(n_RAS_spaces_occ+n_RAS_spaces_virt)
    integer::occ_RAS_el(n_RAS_spaces_occ,MAXVAL(counter)),virt_RAS_el(n_RAS_spaces_virt,MAXVAL(counter))
    integer::verbose

    integer:: i,j,k,l,m,n,p,q,NRD_s_t,NRD_g_t
    integer:: n_el_temp,active_el

   
     
    NRD_s_t=0
    do i=1,n_RAS_spaces_occ
       do j =1,counter(i)
          do m=1,n_RAS_spaces_occ
             if ( m .gt. i) then
                do n =1,counter(m)
                   do k=1,n_RAS_spaces_virt
                      do l=1,counter(n_RAS_spaces_occ+k)
                         do p=1,n_RAS_spaces_virt
                            if ( p .gt. k) then
                               do q =1,counter(p)
                                  NRD_s_t=NRD_s_t+1
                                  RAS_dis_temp(NRD_s_t,i)=occ_RAS_el(i,j)
                                  RAS_dis_temp(NRD_s_t,m)=occ_RAS_el(m,n)
                                  RAS_dis_temp(NRD_s_t,n_RAS_spaces_occ+k+1)=virt_RAS_el(k,l)
                                  RAS_dis_temp(NRD_s_t,n_RAS_spaces_occ+p+1)=virt_RAS_el(p,q)
                               end do
                            end if
                         end do
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do

        
    call flush(6)

    if (verbose .eq. 2) then
       do i=1,NRD_s_t
          write(*,*)"NRD_spin",i,RAS_dis_temp(i,1),RAS_dis_temp(i,2),RAS_dis_temp(i,4),RAS_dis_temp(i,5)
       end do
       call flush(6)
    end if
    
    
   NRD_g_t=NRD_good
    

    do i=1,NRD_s_t
       n_el_temp=0
       do j=1,n_RAS_spaces_occ
          n_el_temp=n_el_temp+RAS_dis_temp(i,j)
       end do
       do k=1,n_RAS_spaces_virt
          n_el_temp=n_el_temp+RAS_dis_temp(i,n_RAS_spaces_occ+1+k)
       end do
       active_el=n_spin-n_el_temp
       if (active_el.ge.0 .AND. active_el.le.active_space) then
          NRD_g_t=NRD_g_t+1
          RAS_dis_temp(i,n_RAS_spaces_occ+1)=active_el
          do j=1,n_RAS_spaces_occ
             RAS_dis(NRD_g_t,j)=RAS_dis_temp(i,j)
          end do
          RAS_dis(NRD_g_t,n_RAS_spaces_occ+1)=RAS_dis_temp(i,n_RAS_spaces_occ+1)
          do k=1,n_RAS_spaces_virt
             RAS_dis(NRD_g_t,n_RAS_spaces_occ+1+k)=RAS_dis_temp(i,n_RAS_spaces_occ+1+k)
          end do
       else
          if (verbose .eq. 2) then
             write(*,*)i,"NIE",i,n_spin-n_el_temp,RAS_dis_temp(i,1),RAS_dis_temp(i,2),RAS_dis_temp(i,4),RAS_dis_temp(i,5)
          end if
       end if
    end do



!    write(*,*)"NRD_spin,NRD_good******",NRD_good,NRD_g_t
!    call flush(6)

    
    if (verbose .eq. 2) then
       do i=NRD_good+1,NRD_g_t
          write(*,*)"NRD_good",i,RAS_dis(i,1),RAS_dis(i,2),RAS_dis(i,3),RAS_dis(i,4),RAS_dis(i,5)
       end do
       call flush(6)
    end if

    NRD_good=NRD_g_t
    

    call flush(6)
    
    
  end subroutine fill_distribution_arrays
  
  


subroutine count_el_nums_in_ras_spaces(count_possibilities,fill,counter,occ_RAS_el,virt_RAS_el,maxcounter,n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,excit_array)
  implicit none
  logical, intent(in)::count_possibilities,fill
  integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
  integer, intent(inout) :: counter(n_RAS_spaces_occ+n_RAS_spaces_virt)
  integer, intent(in)::maxcounter
  integer, intent(inout)::occ_RAS_el(n_RAS_spaces_occ,maxcounter),virt_RAS_el(n_RAS_spaces_virt,maxcounter)
  integer:: i,j,k


  write(*,*)"count_possibilities,fill",count_possibilities,fill

  do i=1,n_RAS_spaces_occ   
     k=0
     do j=excit_array(i),0,-1
        k=k+1
        if (fill .eqv. .true.) then
           occ_RAS_el(i,k+1)=RAS_space_occ(i)-j
        end if
     end do
     if (count_possibilities .eqv. .true.) then
        counter(i)=k
     end if
  end do


  do i=1,n_RAS_spaces_virt
     k=0
     do j=0,excit_array(n_RAS_spaces_occ+i),+1
        k=k+1
        if (fill .eqv. .true.) then 
           virt_RAS_el(i,k)=j
           write(*,*)"VIRT",n_RAS_spaces_occ+i,j
           write(*,*)"IK",i,k,virt_RAS_el(i,k)
           call flush(6)
        end if
     end do
     if (count_possibilities .eqv. .true.) then
        counter(n_RAS_spaces_occ+i)=k
     end if
  end do


    do i=1,n_RAS_spaces_virt
       do j=1,counter(n_RAS_spaces_occ+i)
          write(*,*)"virt_RAS_el(p,l)****",i,j,virt_RAS_el(i,j)
       end do
    end do


  
end subroutine count_el_nums_in_ras_spaces


subroutine calculate_space_size(NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,NRDa,NRDb,sizea,sizeb,size_tot,verbose)
  implicit none
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in)::NRDa(3,2),NRDb(3,2)
  integer, intent(in):: NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta
  integer, intent(in):: RAS_el_array_alpha(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_el_array_beta(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(out)::sizea(3,2),sizeb(3,2),size_tot(3,2) 
  integer, intent(in):: verbose
  
  integer:: i,r,v,a
  integer:: vector_length
  
  ! n_alpha n_beta block
  sizea(1,1)=1
  call flush(6)
  call calc_size_block(NRDa(1,1),NRDa(1,2),NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizea(1,2))
  
  sizeb(1,1)=1 
  call calc_size_block(NRDb(1,1),NRDb(1,2),NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizeb(1,2))

  size_tot(1,1)=1
  size_tot(1,2)=sizea(1,2)*sizeb(1,2)
  

  
  ! n_alpha+1 n_beta-1 block
  sizea(2,1)=1 
  call calc_size_block(NRDa(2,1),NRDa(2,2),NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizea(2,2))

  sizeb(3,1)=1 
  call calc_size_block(NRDb(3,1),NRDb(3,2),NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizeb(3,2))
  
  size_tot(2,1)=size_tot(1,2)+1
  size_tot(2,2)=sizea(2,2)*sizeb(3,2)


  ! n_alpha-1 n_beta+1 block
  sizea(3,1)=1 
  call calc_size_block(NRDa(3,1),NRDa(3,2),NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizea(3,2))

  sizeb(2,1)=1 
  call calc_size_block(NRDb(2,1),NRDb(2,2),NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,sizeb(2,2))
  
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



subroutine calc_size_block(range1,range2,NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,vec_length)
  implicit none
  integer, intent(in):: range1,range2
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in):: NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta
  integer, intent(in):: RAS_el_array_alpha(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_el_array_beta(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(out)::vec_length

  integer:: i,r,v,a
  integer:: siz
  real(8) :: nCr_dp
  write(*,*)"range1,range2",range1,range2
  vec_length=0
  do i=range1,range2
     a=n_RAS_spaces_occ+1 
     siz=1
     do r=1,n_RAS_spaces_occ
        siz=siz*int(nCr_dp(RAS_space_occ(r),RAS_el_array_alpha(i,r)))
     end do
     do v=1,n_RAS_spaces_virt
        siz=siz*int(nCr_dp(RAS_space_virt(v),RAS_el_array_alpha(i,a+v)))
     end do
     siz=siz*int(nCr_dp(active_space,RAS_el_array_alpha(i,a)))
     vec_length=vec_length+siz
  end do
  write(*,*)"vec_length",vec_length
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

  




subroutine fill_spin_strings(n_s,NRD_spin,RAS_el_array_spin,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,range1,range2,sizes,str_s,verbose)
  implicit none
  integer, intent(in):: n_s
  integer, intent(in):: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in):: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in):: range1,range2
  integer, intent(in):: NRD_spin
  integer, intent(in):: RAS_el_array_spin(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer, intent(in)::sizes
  integer, intent(inout)::str_s(sizes,n_s)
  integer, intent(in):: verbose
  integer :: space_sizes(n_RAS_spaces_occ+1+n_RAS_spaces_virt)
  integer, allocatable :: str_temp(:,:)
  integer, allocatable :: another_str_temp(:,:)
  integer :: siz_space, tot_siz_space
  integer :: n_el_temp
  integer:: orbital_index
  integer :: n_str_temp(n_RAS_spaces_occ+1+n_RAS_spaces_virt)
  integer :: j
  integer :: a, i , r ,v
  real(8) :: nCr_dp
      write(*,*) "tutaj",RAS_el_array_spin(61,:),NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1
call flush(6)
stop
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
   if (i .eq. 61) then
   write(*,*) v,int(nCr_dp(RAS_space_virt(v),RAS_el_array_spin(i,a+v))),RAS_space_virt(v),RAS_el_array_spin(i,a+v),nCr_dp(RAS_space_virt(v),RAS_el_array_spin(i,a+v))
   write(*,*) "tuuuu",RAS_el_array_spin(61,:)
call flush(6)
stop
   end if
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
      call string_direct_product(str_temp(1:space_sizes(a),1+n_el_temp:RAS_el_array_spin(i,a)+n_el_temp),space_sizes(a),RAS_el_array_spin(i,a),str_temp(1:n_str_temp(j),1:n_el_temp),n_str_temp(j),n_el_temp,another_str_temp(1:n_str_temp(j+1),1:n_el_temp+RAS_el_array_spin(i,a)))
      str_temp = another_str_temp

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
   if (i .eq. 61) then
   write(*,*) "i", space_sizes
   call flush(6)
   stop
   end if
   end do
write(*,*) "koniec2"
call flush(6)
stop 



   

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
   integer, intent(in):: string_array2(nstr_1,nstr_2)
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

