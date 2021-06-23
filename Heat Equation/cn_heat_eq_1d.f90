!code name: cn_heat_eq_1d.f90
!
!description: This code solves the Heat equation by
!             Crank-Nicolson (CN) scheme and compares
!             the results with the exact solution
!
!     Problem:
!
!            ut=sig*uxx
!
!     Exact solution:
!            u(x,t)=exp(-pi**2*t)*sin(pi*x)
!
! Referece: 1) Prof. Dr. Arnaldo Gammal Notes
!
!           2) An Introduction to the Numerical
!              Solution of Differential Equation
!             (Revised edition)-Douglas Quinney (1987)
!
!date: June 23rd, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!____________________________________________________
!
!define 'global variables'
module variables
implicit none
integer,parameter::Nx=200 !number of grid points
real*8,parameter::Lx=1.0d0  !space
real*8,parameter::dx=Lx/Nx !spatial step dx=0.005
real*8,dimension(1:Nx+1)::x !spatial array
real*8,parameter::dt=0.0005d0 !time step
real*8,parameter::pi=dacos(-1.0d0)
real*8,dimension(1:Nx+1)::u0 !ansatz array
real*8,parameter::sig=1.0d0 !factor of the heat equation
end module variables
!************************************
! main code
program heat_1d_cn
use variables
implicit none
integer,parameter::Nt=40 !number of time steps
real*8,dimension(1:Nx+1)::u !solution
real*8::tf=Nt*dt !tf=0.1
real*8::exact
integer::i,k

!call for initial conditions
!like the x-grid and the ansatz
call initial()
!store the ansatz like a initial guess for
!the solkutin u
u=u0

!time propagation
do k=1,Nt
   call CNx(u,dt)
end do

!plot the solutions at the time tf
open(1,file="CN_solution.txt")
open(2,file="exact_solution.txt")

do i=1,Nx+1
   write(1,*)x(i),u(i)
   write(2,*)x(i),exact(x(i),tf)
end do

end program heat_1d_cn
!************************************
! exact solution of the heat equation
function exact(x0,t0)
use variables, only:pi
implicit none
real*8::exact
real*8::x0,t0

exact=dsin(pi*x0)*dexp(-t0*pi**2)

end function exact
!************************************
! define some varibles like
! x-grid and ansatz
subroutine initial()
use variables, only:pi,Nx,dx,x,u0
implicit none
integer::i

!define the spatial grid x=[0,Lx]
do i=1,Nx+1
   x(i)=(i-1)*dx
end do

!define the ansatz u0(x)
do i=1,Nx+1
   u0(i)=dsin(pi*x(i))
end do

!insert the boundary conditions
!u0(0)=u0(Lx)=0

u0(1)=0.0d0
u0(Nx+1)=0.0d0

end subroutine initial
!************************************
! Crank Nicolson routine
!
subroutine CNx(u,dt0)
use variables, only:Nx,dx,sig
implicit none
real*8,dimension(1:Nx+1)::u !solution
real*8,dimension(1:Nx+1)::fx !auxiliary function
real*8::dt0 ! time-step can be 'dt' or smaller
real*8::rx  ! CN parameter
real*8,dimension(1:Nx+1)::a,b,c ! Matrix A elements
real*8,dimension(1:Nx+1)::l,v,w ! Matrix L and U elements
real*8,dimension(1:Nx+1)::d !right side of the equation (old solution)
real*8,dimension(1:Nx+1)::z ! forward substitution auxiliary variable
integer::i

rx=dt0*sig/(dx**2)

!insert the boundary conditions
!u(0)=u(Lx)=0
u(1)=0.0d0
u(Nx+1)=0.0d0

! store the solution in a auxiliary variable
do i=1,Nx+1
   fx(i)=u(i)
end do

!---------------------------------------------
!define the non-zero elements of the matrix A
!define upper diagonal
do i=2,Nx
   c(i)=-0.5d0*rx
end do
!define the main diagonal
do i=2,Nx
   b(i)=1.0d0+rx
end do
!define the lower diagonal
do i=2,Nx
   a(i)=-0.5d0*rx
end do

!------------------------------------------
! define L and U matrix elements

do i=2,Nx
   v(i)=c(i)
end do

w(2)=b(2)
do i=3,Nx
   l(i)=a(i)/w(i-1)
   w(i)=b(i)-l(i)*v(i-1)
end do

!-------------------------------------

! older (right) side
do i=2,Nx
   d(i)=fx(i)+0.5d0*rx*(fx(i+1)-2*fx(i)+fx(i-1))
end do

!forward substitution
z(2)=d(2)
do i=3,Nx
      z(i)=d(i)-z(i-1)*l(i)
end do

!get the solution
!by backward substitution
fx(Nx)=z(Nx)/w(Nx)
do i=Nx-1,2,-1
   fx(i)=(z(i)-c(i)*fx(i+1))/w(i)
end do

!return and update the original solution
do i=2,Nx
   u(i)=fx(i)
end do

!insert the boundary conditions again
!u(0)=u(Lx)=0.0d0
u(1)=0.0d0
u(Nx+1)=0.0d0

!obs: Some of the variables here can be difined
!     only one time instead of every time the
!     routine is called, because they don't
!     change with 'u(x)', but here we write
!     them everytime to make more easy to
!     undertand.

end subroutine CNx
