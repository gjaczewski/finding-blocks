! Module for calculation of Green function blocks
module gf_mod 
contains 
function compute_part_gf(n_alpha,n_beta) result(partitions)
    implicit none
    integer, intent(in) :: n_alpha, n_beta !Numbers of electrons
    integer, allocatable :: partitions(:,:)
    integer :: i,j,k !Loop variables
    k = 0
    ! Loop below calculates number of blocks
        do i = -2,2
            do j = -2,2
                if ((abs(i+j) == 1) .and. (n_alpha + i >= 0) .and. (n_beta + j >=0)) then
                    k = k + 1
                end if
            end do
        end do
    allocate(partitions(k,2))
    k = 1
    !Loop below saves blocks to array
    do i = -2,2
        do j = -2,2
            if ((abs(i+j) == 1) .and. (n_alpha + i >= 0) .and. (n_beta + j >=0)) then
                partitions(k,1) = n_alpha + i 
                partitions(k,2) = n_beta  + j
                k = k + 1
            end if 
        end do
    end do
end function compute_part_gf
end module gf_mod




module wf_mod
contains

function compute_part_wf(n_alpha, n_beta) result(partitions)
    implicit none
    integer, intent(in) :: n_alpha, n_beta
    integer, allocatable :: partitions(:,:)
    integer :: i, j
    i = 1

    if (n_alpha >= 1) i = i + 1
    if (n_beta >= 1) i = i + 1
    allocate(partitions(i,2))
    if (n_alpha==0) then
     do j = 0, 1
        partitions(j+1,1) = n_alpha + j
        partitions(j+1,2) = n_beta - j
     end do
    end if
    if (n_beta==0) then
     do j = 0, 1
        partitions(j+1,1) = n_alpha - j
        partitions(j+1,2) = n_beta + j
     end do
    end if
    if (n_beta*n_alpha>0) then
        do j=-1,1
            partitions(j+2,1) = n_alpha - j
            partitions(j+2,2) = n_beta + j
        end do
    end if
end function compute_part_wf
end module wf_mod


program prototype
    use gf_mod
    use wf_mod
    implicit none
    integer :: n_alpha, n_beta
    integer :: k
    integer, allocatable :: gf_partitions(:,:), wf_partitions(:,:)

    print *, 'Please enter n_\alpha and n_\beta. '
    read(*,*) n_alpha, n_beta
    gf_partitions = compute_part_gf(n_alpha, n_beta)
    wf_partitions = compute_part_wf(n_alpha,n_beta)
    print*, 'Wave function function blocks:'
    do k=1,size(wf_partitions,dim=1)
        print*,wf_partitions(k,:)
    end do
    print*, 'Green function blocks:'
    do k=1,size(gf_partitions,dim=1)
        print*,gf_partitions(k,:)
    end do
end program prototype
