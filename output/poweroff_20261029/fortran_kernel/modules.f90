module elems
 implicit none
 type element
  integer :: node(4)
  real(8) :: cc(4),dd(4),ee(4),de
 end type
 integer :: nelem=1
 type(element) :: elem(1)
end module
module nodes
 implicit none
 type nodetype
  real(8) :: xyz(3)
 end type
 type(nodetype) :: node(4)
end module
module thermos
 implicit none
 type material
  real(8) :: kx,ky,kz,rho,cv
 end type
 real(8) :: dt=1d-9
 type(material) :: matprop(1)
end module
