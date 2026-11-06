program probe
 use elems
 use nodes
 use thermos
 implicit none
 integer :: n,i,j,k,l
 real(8),external :: gtemp0,gtemp2
 read(*,*) n
 elem(1)%node=(/1,2,3,4/)
 matprop(1)%kx=.01d0;matprop(1)%ky=.01d0;matprop(1)%kz=.008d0
 matprop(1)%rho=1d-10/8.854187d-12;matprop(1)%cv=0
 do i=1,n
  do j=1,4
   read(*,*) node(j)%xyz
  end do
  call cdecal()
  do j=1,4
   do k=1,4
    write(*,'(2es25.16)')gtemp0(1,j,k,1),gtemp2(1,j,k,1)
   end do
  end do
 end do
end program
