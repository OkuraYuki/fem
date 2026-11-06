module elems
 implicit none
 type element
  integer :: node(4)
  real(8) :: cc(4),dd(4),ee(4),de
 end type
 integer :: nelem
 type(element),allocatable :: elem(:)
end module
module nodes
 implicit none
 type nodetype
  real(8) :: xyz(3)
 end type
 type(nodetype),allocatable :: node(:)
end module
module thermos
 implicit none
 type material
  real(8) :: kx,ky,kz,rho,cv
 end type
 real(8) :: dt=1d-9
 type(material) :: matprop(3)
end module
program full_residual
 use elems
 use nodes
 use thermos
 implicit none
 integer :: n,unused1,unused2,i,j,k,l,nb,nr,nm,n2,nz,lo,hi,mat
 integer,allocatable :: material_id(:)
 logical,allocatable :: fixed(:)
 real(8) :: fact,dummy,vals(6),kval,cval
 real(8),external :: gtemp0,gtemp2
 real(8),allocatable :: p0(:),p24(:),p25(:),ks0(:),ke0(:),r(:),rhs(:)
 open(10,file='models/model_entyu_1/in.dat',status='old')
 read(10,*)n,nelem,unused1,unused2,fact
 allocate(elem(nelem),node(n),material_id(nelem),fixed(n))
 read(10,*)((elem(i)%node(j),j=1,4),i=1,nelem)
 read(10,*)((node(i)%xyz(j),j=1,3),i=1,n)
 close(10)
 do i=1,n
  node(i)%xyz=node(i)%xyz/fact
 end do
 open(10,file='models/model_entyu_1/sin.dat',status='old')
 read(10,*)nb,nr,nm,n2,nz
 do i=1,nb
  read(10,*)unused1,unused2,j,dummy
 end do
 do i=1,nr
  read(10,*)lo,hi,mat
  material_id(lo:hi)=mat
 end do
 do i=1,nm
  read(10,*)vals
  matprop(i)%kx=vals(1);matprop(i)%ky=vals(2);matprop(i)%kz=vals(3)
  matprop(i)%rho=vals(4)/8.854187d-12;matprop(i)%cv=0
 end do
 close(10)
 allocate(p0(n),p24(n),p25(n),ks0(n),ke0(n),r(n),rhs(n))
 open(11,file='output/poweroff_20261029/full/phi0.bin',access='stream',form='unformatted');read(11)p0;close(11)
 open(11,file='output/poweroff_20261029/full/phi24.bin',access='stream',form='unformatted');read(11)p24;close(11)
 open(11,file='output/poweroff_20261029/full/phi25.bin',access='stream',form='unformatted');read(11)p25;close(11)
 call cdecal()
 ks0=0;ke0=0;r=0;rhs=0;fixed=.false.
 do i=1,nelem
  mat=material_id(i)
  if(mat==1.or.mat==3)fixed(elem(i)%node)=.true.
  do j=1,4
   k=elem(i)%node(j)
   do l=1,4
    kval=gtemp0(i,j,l,mat);cval=gtemp2(i,j,l,mat)
    ks0(k)=ks0(k)+kval*p0(elem(i)%node(l))
    ke0(k)=ke0(k)+cval*dt*p0(elem(i)%node(l))
    r(k)=r(k)+(kval+cval)*p25(elem(i)%node(l))-cval*p24(elem(i)%node(l))
    rhs(k)=rhs(k)+cval*p24(elem(i)%node(l))
   end do
  end do
 end do
 write(*,'(a,i0)')'nodes ',n
 write(*,'(a,i0)')'elements ',nelem
 write(*,'(a,es25.16)')'initial_Ksigma_free_max_SI ',maxval(abs(ks0),mask=.not.fixed)
 write(*,'(a,es25.16)')'initial_Kepsilon_free_max_SI ',maxval(abs(ke0),mask=.not.fixed)
 write(*,'(a,es25.16)')'off25_residual_free_max_SI ',maxval(abs(r),mask=.not.fixed)
 write(*,'(a,es25.16)')'off25_rhs_free_max_SI ',maxval(abs(rhs),mask=.not.fixed)
 write(*,'(a,es25.16)')'off25_residual_free_L2_over_rhs_L2 ', &
 sqrt(sum(r*r,mask=.not.fixed)/sum(rhs*rhs,mask=.not.fixed))
end program
