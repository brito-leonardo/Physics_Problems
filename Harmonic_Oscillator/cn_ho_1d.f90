!code name: cn_ho_1d.f90
!
!description: This code solves the quantum harmonic
!             oscillator (HO) by Crank-Nicolson (CN)
!             scheme with imaginary-time progation.
!             We find the ground state wave function
!             and its corresponding energy E.
!             We compare the results with the exact
!             solution.
!
!     Problem:
!
!            i*ut=-0.5*uxx+Vtrap(x)*u, Vtrap(x)=0.5*x**2
!
!     Exact solution:

!     wave function:  u(x)=exp(-0.5*x**2)/sqrt(sqrt(pi))
!
!     energy:        E=0.5
!
!
!           2) An Introduction to the Numerical
!              Solution of Differential Equation
!             (Revised edition)-Douglas Quinney (1987)
!
!date: June 27th, 2021
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
real*8,parameter::Lx=16.0d0  !space
real*8,parameter::dx=Lx/Nx !spatial step dx=0.005
real*8,dimension(1:Nx+1)::x !spatial array
real*8,parameter::pi=dacos(-1.0d0)
complex*16,dimension(1:Nx+1)::u0 !ansatz array
real*8,dimension(1:Nx+1)::Vtrap !harmonic oscillator trap potencial
real*8,parameter::sig=-0.5d0 !factor of the heat equation
complex*16,parameter::ci=(0.0d0,1.0d0), cr=(1.0d0,0.0d0)
end module variables
!************************************
! main code
program ho_1d_cn
use variables, only:Nx,x,u0
implicit none
complex*16,dimension(1:Nx+1)::psi
real*8,dimension(1:Nx+1)::psi0
real*8::E0,E,res,Vir
real*8::error
integer::option
integer::i

!call for initial conditions
!like the x-grid and the ansatz
call initial()
!store the ansatz like a initial guess for
!the solution psi
psi=u0

! Find the ground-state
! option=1 means imaginary-time propagation
option=1
call ground_state(psi,option)
!return the energy
call energy(psi,E,res,vir)
!call exact solution and its energy
call exact(psi0,E0)

open(100,file="density_exact.txt")
open(200,file="density_CN.txt")
do i=1,Nx+1
   write(100,*)x(i),psi0(i)**2
   write(200,*)x(i),dreal(dconjg(psi(i))*psi(i))
end do
open(300,file="energy_comparative.txt")
write(300,*)'Energy exact: ',E0
write(300,*)'Energy CN method: ',E
write(300,*)'error: ', dabs(E-E0)

!test you really found a stationary state
!use a long real time for that. The state
!have to return the same energy for all times

option=2
call propag(psi,option)


end program ho_1d_cn
!***********************************************************************
subroutine ground_state(psi,option)
use variables, only:Nx
implicit none
integer::option
complex*16,dimension(1:Nx+1)::psi !wave function
integer,parameter::Nt=10 !number of sequence of 10.000 time steps
integer,parameter::Ntt=10000 !number of each sequence of time steps
real*8,parameter::dt=0.0005d0
real*8::E,res,vir
real*8::Eold,error
real*8::norma
integer::i,k,kk


open(1,file="data_ground_state.txt")

Eold=1.0d0
   !time propagation
do k=1,Nt
   do kk=1,Ntt
      !first nonderivative half-step
      call ND(psi,dt/2,option)
      !derivative Crank-Nicolson step
      call CNx(psi,dt,option)
      !second nonderivative half-step
      call ND(psi,dt/2,option)
      !renormalize teh wave function
      call norm(psi,norma)
      psi=psi/dsqrt(norma)
   end do
  !energy convergence criterium to finish the loop
   call energy(psi,E,res,vir)
   error=dabs(E-Eold)
   Eold=E
   write(1,991)k*Ntt*dt,E,error,res,Vir
   if((error<1.d-4).and.(res<1.d-3).and.(vir<1.d-3))then
      exit
   end if
end do

991 format(5f10.5)

end subroutine ground_state
!***********************************************************************
subroutine propag(psi,option)
use variables, only:Nx
implicit none
integer::option
complex*16,dimension(1:Nx+1)::psi
integer,parameter::Nt=10, Ntt=100000
real*8,parameter::dt=0.001d0
real*8::E,res,vir,norma
integer::k,kk


open(2,file="data_propagation.txt")

   !time propagation
do k=1,Nt
   do kk=1,Ntt
      !first nonderivative half-step
      call ND(psi,dt/2,option)
      !derivative Crank-Nicolson step
      call CNx(psi,dt,option)
      !second nonderivative half-step
      call ND(psi,dt/2,option)
   end do
   !observe the energy and norm, they have to be conserved
   call energy(psi,E,res,vir)
   call norm(psi,norma)
   write(2,992)k*Ntt*dt,norma,E
end do

992 format(3f10.5)


end subroutine propag
!***********************************************************************
! exact solution of the heat equation
subroutine exact(psi0,E0)
use variables, only:Nx,x,pi
implicit none
real*8,dimension(1:Nx)::psi0
real*8::E0
integer::i

!exact wave function solution
do i=1,Nx+1
   psi0(i)=dexp(-0.5d0*x(i)**2)/dsqrt(dsqrt(pi))
end do

!exact energy
E0=0.5d0

end subroutine exact
!************************************
! define some varibles like
! x-grid, ansatz and the trap potential
subroutine initial()
use variables, only:pi,Nx,Vtrap,dx,x,Lx,u0
implicit none
real*8::norma
integer::i

!define the spatial grid x=[-Lx/2,+Lx/2]
x(1)=-Lx/2.0d0
do i=2,Nx+1
   x(i)=x(1)+(i-1)*dx
end do

!define the ansatz u0(x)
do i=1,Nx+1
   u0(i)=dexp(-x(i)**2)
end do
!normalize it
call norm(u0,norma)
u0=u0/dsqrt(norma)

!define the trap potential
do i=1,Nx+1
   Vtrap(i)=0.5d0*x(i)**2
end do



end subroutine initial
!****************************************************************
!  non-derivative problem
!
!    c*ut=V*u
!
!    where
!     c=i (real-time)
!     c=-1 (imaginary time)
!
!
subroutine ND(psi,dt0,option)
use variables, only:Nx,Vtrap,ci,cr
implicit none
complex*16,dimension(1:Nx+1)::psi !wave function
complex*16,dimension(1:Nx+1)::V  !potential
complex*16,parameter::cte_RE=ci,cte_IM=-cr !auxiliary constants
complex*16::cte
real*8::dt0 !time-step parameter
integer::option !states the nature of the time (ral or imaginary)
integer::i


if(option==1)then
   cte=cte_IM
end if
if(option==2)then
   cte=cte_RE
end if

do i=1,Nx+1
   V(i)=Vtrap(i)
end do

do i=1,Nx+1
   psi(i)=cdexp((V(i)*dt0)/cte)*psi(i)
end do

end subroutine ND
!*****************************************************************
! Crank Nicolson routine
!
! Suppose the matrix equation:
!     (left side)      A*x=d  (right side)
!
! where A is a tridiagonal matrix (blanck spaces are zero elements):
!
!
!          (b(2) c(2)                             )
!          (a(2) b(3) c(3)                        )
!          (    a(3)  b(4) c(4)                   )
!          (          a(4) b(5) c(5)              )
!    A =   (                                      )
!          (                                      )
!          (               a(Nx-1) b(Nx-1) c(Nx-1))
!          (                        a(Nx)   b(Nx) )

! The matrx A can be written as A=L*U
! where L and U are matrix too, 'lower' and 'upper' diagonals
! then we can write A*x=d as L*U*x=d , with
!
!          (1                                 )
!          (l(2) 1                            )
!          (    l(3)  1                       )
!          (         l(4) 1                   )
!    L =   (                                  )
!          (                                  )
!          (               l(Nx-1)   1        )
!          (                        l(Nx)   1 )
!
!
!          (w(2) v(2)                             )
!          (     w(3) v(3)                        )
!          (          w(4) v(4)                   )
!          (               w(5) v(5)              )
!    U =   (                                      )
!          (                                      )
!          (                      w(Nx-1)  v(Nx-1))
!          (                                 w(Nx))
!
!    define z=U*x, we get L*z=d, we solve it by  forward substitution
!    and later the U*x=z by backrward substitution, returning x.
!
!    We have to solve the problem:
!
!                         (left side) =(right side)
!                               A *x  =  d
!
!  u(i)-0.5*rx*(u(i+1)-2*u(i)+u(i-1))=ui+0.5*rx*(u(i+1)-2*u(i)+u(i-1))
!           (new solution)                    (old solution)
!
subroutine CNx(psi,dt0,option)
use variables, only:Nx,dx,sig,ci,cr
implicit none
integer::option
complex*16,dimension(1:Nx+1)::psi !solution
complex*16,dimension(1:Nx+1)::fx !auxiliary function
real*8::dt0 ! time-step can be 'dt' or smaller
complex*16::rx  ! CN parameter
complex*16,dimension(1:Nx+1)::a,b,c ! Matrix A elements
complex*16,dimension(1:Nx+1)::l,v,w ! Matrix L and U elements
complex*16,dimension(1:Nx+1)::d !right side of the equation (old solution)
complex*16,dimension(1:Nx+1)::z ! forward substitution auxiliary variable
complex*16,parameter::cte_RE=ci,cte_IM=-cr
complex*16::cte
integer::i


if(option==1)then
   cte=cte_IM
end if
if(option==2)then
   cte=cte_RE
end if

rx=(dt0*sig)/(cte*(dx**2))

!insert the boundary conditions
!u(0)=u(Lx)=0
psi(1)=0.0d0
psi(Nx+1)=0.0d0

! store the solution in a auxiliary variable
do i=1,Nx+1
   fx(i)=psi(i)
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
   psi(i)=fx(i)
end do

!insert the boundary conditions again
!u(0)=u(Lx)=0.0d0
psi(1)=0.0d0
psi(Nx+1)=0.0d0

!obs: Some of the variables here can be defined
!     only one time instead of every time the
!     routine is called, because they don't
!     change with 'u(x)', but here we write
!     them everytime to make more easy to
!     undertand.

end subroutine CNx
!**************************************************
subroutine norm(psi,norma)
use variables, only:Nx,dx
implicit none
complex*16,dimension(1:Nx+1)::psi
real*8,dimension(1:Nx+1)::den
integer,dimension(1:Nx+1)::fx
real*8::norma
integer::i

!simpson factors
fx(1)=1
do i=2,Nx,2
   fx(i)=4
end do
do i=3,Nx-1,2
   fx(i)=2
end do
fx(Nx+1)=1

norma=0.0d0
!simpson integral method
do i=1,Nx+1
      den(i)=dreal(dconjg(psi(i))*psi(i))
      norma=norma+(fx(i)*dx*den(i))/3.0d0
end do

end subroutine norm
!***********************************************
subroutine energy(psi,E,res,vir)
use variables, only:Nx,dx,Vtrap
implicit none
complex*16,dimension(1:Nx+1)::psi,psi_dx
real*8,dimension(1:Nx+1)::denk,denp
complex*16,dimension(1:Nx+1)::res_vec
integer,dimension(1:Nx+1)::fx
real*8::Ekin,Epot,E
real*8::vir
real*8::res
integer::i


!simpson factors
fx(1)=1
do i=2,Nx,2
   fx(i)=4
end do
do i=3,Nx-1,2
   fx(i)=2
end do
fx(Nx+1)=1


!first derivative
!boundary conditions
psi_dx(1)=(0.0d0,0.0d0)
psi_dx(Nx+1)=(0.0d0,0.0d0)
!centered differences
do i=2,Nx
   psi_dx(i)=(psi(i+1)-psi(i-1))/(2*dx)
end do
!kinetic and potential energies
Epot=0.0d0
Ekin=0.0d0
!simpson integral method
do i=1,Nx+1
      !kinetic energy
      denk(i)=dreal(dconjg(psi_dx(i))*psi_dx(i))
      Ekin=Ekin+(fx(i)*dx*0.5d0*denk(i))/3.0d0
      !potential energy
      denp(i)=dreal(dconjg(psi(i))*psi(i))
      Epot=Epot+(fx(i)*dx*Vtrap(i)*denp(i))/3.0d0
end do
!total energy
E=Epot+Ekin
!virial theorem
vir=dabs(-Ekin+Epot)

!residue
do i=2,Nx+1
   res_vec(i)=-0.5d0*(psi(i+1)-2*psi(i)+psi(i-1))/(dx**2) &
             &+Vtrap(i)*psi(i) &
             &-E*psi(i)
end do
res=maxval(zabs(res_vec))


end subroutine energy
