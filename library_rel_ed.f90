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




subroutine fill_RAS_electron_array(relativistic,norb,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,n_alpha,n_beta,NRD_spin_tmp,NRD_spin_alpha,NRD_spin_beta,RAS_el_array_alpha,RAS_el_array_beta)
  
  implicit none
  logical, intent(in) :: relativistic
  integer, intent(in) :: n_alpha, n_beta
  integer, intent(in) :: norb,n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt)
  integer, intent(in) :: excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
  integer, intent(inout):: NRD_spin_alpha,NRD_spin_beta
  integer, intent(in) :: NRD_spin_tmp
  integer, intent(inout):: RAS_el_array_alpha(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_el_array_beta(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
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



  
  ! for alpha electrons
  call fill_RAS_spaces_per_spin(n_alpha,NRD_spin_tmp,NRD_spin_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_alpha)

  STOP
  ! for beta electrons
!  call fill_RAS_spaces_per_spin(n_beta,NRD_spin_tmp,NRD_spin_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_beta)

  if (relativistic == .true.) then
     
  ! for alpha+1 electrons
     call fill_RAS_spaces_per_spin(n_alpha+1,NRD_spin_tmp,NRD_spin_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_alpha)
    ! for alpha-1 electrons   
       call fill_RAS_spaces_per_spin(n_alpha-1,NRD_spin_tmp,NRD_spin_alpha,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_alpha)

  
  ! for beta+1 electrons
!  call fill_RAS_spaces_per_spin(n_beta+1,NRD_spin_tmp,NRD_spin_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_beta)


    ! for beta-1 electrons
!  call fill_RAS_spaces_per_spin(n_beta-1,NRD_spin_tmp,NRD_spin_beta,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_el_array_beta)

end if


end subroutine fill_RAS_electron_array
  
  
subroutine fill_RAS_spaces_per_spin(n_spin,NRD_spin_tmp,NRD_spin,n_RAS_spaces_occ,RAS_space_occ,n_RAS_spaces_virt,RAS_space_virt,active_space,excit_array,RAS_dis_temp,RAS_dis)
  implicit none
  integer, intent(in) :: n_spin
  integer, intent(in) :: n_RAS_spaces_occ,n_RAS_spaces_virt,active_space
  integer, intent(in) :: RAS_space_occ(n_RAS_spaces_occ),RAS_space_virt(n_RAS_spaces_virt),excit_array(n_RAS_spaces_occ+n_RAS_spaces_virt)
  integer, intent(in) :: NRD_spin_tmp
  integer, intent(inout):: NRD_spin
  integer, intent(inout):: RAS_dis(NRD_spin_tmp,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
  integer:: verbose
  

  
  integer, dimension(:,:), allocatable::occ_RAS_el,virt_RAS_el,RAS_dis_temp
  integer, dimension(:), allocatable::counter
  integer:: n_el_temp,active_el,NRD_good
  integer:: i,j,k,l,m,n,p,q
  logical:: count_possibilities,fill
  double precision, dimension(:), allocatable::scratch
  !NRD_spin - number of possible electron distributions for a given alpha or beta spin
  verbose=2

  allocate(scratch(10000))
  
  allocate(counter(n_RAS_spaces_occ+n_RAS_spaces_virt))
  count_possibilities=.true.
  counter(:)=0
  call  count_el_nums_in_ras_spaces(count_possibilities,fill,counter,scratch(1),scratch(1),MAXVAL(counter),n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,excit_array)
  
  allocate(occ_RAS_el(n_RAS_spaces_occ,MAXVAL(counter)))
  allocate(virt_RAS_el(n_RAS_spaces_virt,MAXVAL(counter)))
  write(*,*)"MAXVAL(counter)",MAXVAL(counter)
  
  fill=.true.
  call  count_el_nums_in_ras_spaces(count_possibilities,fill,counter,occ_RAS_el,virt_RAS_el,MAXVAL(counter),n_RAS_spaces_occ,n_RAS_spaces_virt,RAS_space_occ,excit_array)

  NRD_spin=1
  do i=1,n_RAS_spaces_occ 
     NRD_spin=NRD_spin*counter(i)
  end do
  do k=1,n_RAS_spaces_virt
     NRD_spin=NRD_spin*counter(n_RAS_spaces_occ+k)
  end do

  allocate(RAS_dis_temp(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1))
  
  call fill_distribution_arrays(n_spin,NRD_spin,NRD_good,RAS_dis_temp,RAS_dis,n_RAS_spaces_occ,n_RAS_spaces_virt,occ_RAS_el,virt_RAS_el,counter,verbose)

  deallocate(RAS_dis_temp)
end subroutine fill_RAS_spaces_per_spin
  
subroutine fill_distribution_arrays(n_spin,NRD_spin,NRD_good,RAS_dis_temp,RAS_dis,n_RAS_spaces_occ,n_RAS_spaces_virt,occ_RAS_el,virt_RAS_el,counter,verbose)
  implicit none
    integer::n_spin,NRD_spin,NRD_good,n_RAS_spaces_occ,n_RAS_spaces_virt
    integer::RAS_dis_temp(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1),RAS_dis(NRD_spin,n_RAS_spaces_occ+n_RAS_spaces_virt+1)
    integer::counter(n_RAS_spaces_occ+n_RAS_spaces_virt)
    integer::occ_RAS_el(n_RAS_spaces_occ,MAXVAL(counter)),virt_RAS_el(n_RAS_spaces_virt,MAXVAL(counter))
    integer::verbose

    integer:: i,j,k,l,m,n,p,q
    integer:: n_el_temp,active_el
    
    write(*,*)"przed gigantyczna petla"
    call flush(6)
    RAS_dis_temp(:,:)=-1000
    NRD_spin=0
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
                                  NRD_spin=NRD_spin+1
                                  RAS_dis_temp(NRD_spin,i)=occ_RAS_el(i,j)
                                  RAS_dis_temp(NRD_spin,m)=occ_RAS_el(m,n)
                                  RAS_dis_temp(NRD_spin,n_RAS_spaces_occ+k+1)=virt_RAS_el(k,l)
                                  RAS_dis_temp(NRD_spin,n_RAS_spaces_occ+p+1)=virt_RAS_el(p,q)
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
       do i=1,NRD_spin
          write(*,*)"NRD_spin",i,RAS_dis_temp(i,1),RAS_dis_temp(i,2),RAS_dis_temp(i,4),RAS_dis_temp(i,5)
       end do
       call flush(6)
    end if
    
    
    NRD_good=0
    do i=1,NRD_spin
       n_el_temp=0
       do j=1,n_RAS_spaces_occ
          n_el_temp=n_el_temp+RAS_dis_temp(i,j)
          write(*,*)i,"o",j,RAS_dis_temp(i,j)
       end do
       do k=1,n_RAS_spaces_virt
          n_el_temp=n_el_temp+RAS_dis_temp(i,n_RAS_spaces_occ+1+k)
          write(*,*)i,"v",n_RAS_spaces_occ+1+k,RAS_dis_temp(i,n_RAS_spaces_occ+1+k)
       end do
       if ((n_spin-n_el_temp).ge.0) then
          active_el=n_spin-n_el_temp
          NRD_good=NRD_good+1
          RAS_dis_temp(i,n_RAS_spaces_occ+1)=active_el
          write(*,*)i,"a",active_el
          do j=1,n_RAS_spaces_occ
             RAS_dis(NRD_good,j)=RAS_dis_temp(i,j)
          end do
          RAS_dis(NRD_good,n_RAS_spaces_occ+1)=RAS_dis_temp(i,n_RAS_spaces_occ+1)
          do k=1,n_RAS_spaces_virt
             RAS_dis(NRD_good,n_RAS_spaces_occ+1+k)=RAS_dis_temp(i,n_RAS_spaces_occ+1+k)
          end do
          call flush(6)
          !do j=1,n_RAS_spaces_occ
          !   RAS_el_array_good(NRD_good,i)=RAS_dis_temp(i,j)
          !end do
          !RAS_el_array_good(NRD_good,n_RAS_spaces_occ+1)=active_el
          !do j=1,n_RAS_spaces_virt
          !   RAS_el_array_good(NRD_good,n_RAS_spaces_occ+1+j)=RAS_dis_temp(i,n_RAS_spaces_occ+1+j)
          !end do
       else
          write(*,*)i,"NIE",n_spin-n_el_temp
       end if
    end do
    write(*,*)"NRD_spin,NRD_good",NRD_spin,NRD_good


    
    if (verbose .eq. 2) then
       do i=1,NRD_good
          write(*,*)"NRD_good",i,RAS_dis(i,1),RAS_dis(i,2),RAS_dis(i,4),RAS_dis(i,5)
       end do
       call flush(6)
    end if
    
    

    call flush(6)
    STOP

    
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


 
  do i=1,n_RAS_spaces_occ   
     k=0
     do j=excit_array(i),0,-1
        if (fill ==.true.) then
           write(*,*)"Jestem tu",i,k+1
           call flush(6)           
           occ_RAS_el(i,k+1)=RAS_space_occ(i)-j
           write(*,*)"OCC",i,RAS_space_occ(i)-j
           call flush(6)
        end if
        k=k+1
     end do
     if (count_possibilities == .true.) then
        counter(i)=k
        write(*,*)"OCCcounter",i,k
     end if
  end do


  do i=1,n_RAS_spaces_virt
     k=0
     do j=0,excit_array(n_RAS_spaces_occ+i),+1
        if (fill ==.true.) then 
           virt_RAS_el(i,k)=j
           write(*,*)"VIRT",n_RAS_spaces_occ+i,j
           call flush(6)
        end if
        k=k+1
     end do
     if (count_possibilities == .true.) then
        counter(n_RAS_spaces_occ+i)=k
        write(*,*)"VIRTcounter",n_RAS_spaces_occ+i,k
     end if
  end do
  
end subroutine count_el_nums_in_ras_spaces



subroutine unique_elements(size,arr,uniq,count)
    implicit none
    integer, intent(in):: size
    integer, intent(out)::count
    integer, intent(in):: arr(size)
    integer, intent(out)::uniq(size)
    integer :: i, j, n
    logical :: is_new

    ! Example input array

    count = 0

    
    do i = 1, size
       write(*,*)"ARR",i,arr(i)
       is_new = .true.
        do j = 1, count
            if (abs(arr(i)-arr(j))==0) then
                is_new = .false.
                exit
            end if
        end do

        if (is_new) then
            count = count + 1
            uniq(count) = arr(i)
        end if
    end do

    ! Print results
    print *, 'Original array: ', arr
    print *, 'Unique elements:', uniq(1:count)
  end subroutine unique_elements
