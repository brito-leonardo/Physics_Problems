!code name: LU_scheme_short.f90
!
!description: This code solves the tridiagonal matrix system A*x=d
!             using LU decomposition (short version-without explicit matrices)
!
!     Exercise 1.4.3: Solve the system of equations
!
!                            A*x = d
!
!          ( 5  -1   0   0) (x1)   (4.3)
!          (-1   5  -1   0) (x2) = (3.8)
!          ( 0  -1   5  -1) (x3)   (3.1)
!          ( 0   0  -1   5) (x4)   (4.9)
!
!      We want to find x = (x1 x2  x3  x4)
!      the solution is : x= (1.1  1.2  1.1  1.2)
!
!
! Referece: 1) An Introduction to the Numerical
!              Solution of Differential Equation
!              Section 1.4.3. Tridiagonal Matrices, p.31
!             (Revised edition)-Douglas Quinney (1987)
!
!date: June 26th, 2021
!author: Leonardo Brito
!Physics Phd Student
!Institute of Physics (IFUSP)
!University of SÆo Paulo (USP)
!____________________________________________________
program LU_ex
implicit none
integer,parameter::N=4
real,dimension(1:N,1:N)::A,L,U
real,dimension(1:N)::d,z,x,aa,b,c,ll,uu,v
integer::i,j


write(*,*)"We have the problem 'A*x=d'"
write(*,*)
write(*,*)

!define matrix A
call initial(N,A,d)
!write the matrix A at terminal
write(*,*)"with 'A':"
do i=1,N
   write(*,*)(A(i,j), j=1,N)
end do
write(*,*)!add an empty line

!write the vector v at terminal
write(*,*)"and 'd':"
do i=1,N
   write(*,*)d(i)
end do
write(*,*)!add an empty line

!get the L and U matrices (like example 1.4.2, from Quinney 1987)
call LU_routine(N,A,d,x)

!write the vector x (solution) at the terminal
write(*,*)"The solution is 'x':"
do i=1,N
   write(*,*)x(i)
end do

pause

end program LU_ex
!*****************************************
subroutine initial(N,A,d)
implicit none
integer::N
real,dimension(1:N)::d
real,dimension(1:N,1:N)::A

A(1,1)=5.0; A(1,2)=-1.0;A(1,3)=0.0;A(1,4)=0.0
A(2,1)=-1.0; A(2,2)=5.0;A(2,3)=-1.0;A(2,4)=0.0
A(3,1)=0.0; A(3,2)=-1.0;A(3,3)=5.0;A(3,4)=-1.0
A(4,1)=0.0; A(4,2)=0.0;A(4,3)=-1.0;A(4,4)=5.0

d(1)=4.3;d(2)=3.8;d(3)=3.1;d(4)=4.9

end subroutine initial
!*******************************************
! We want to write A=LU, since
!
!
!     (b1   c1                          )
!  A= (aa2  b2  c2                      )
!     (      .   .   .                  )
!     (          .   .       .          )
!     (               .       .      .  )
!     (             aa_n-1  b_n-1  c_n-1)
!     (                      aa_n   c_n )
!
!  and also
!
!   (1                )         (uu1  v1                 )
!L= (ll2  1           )  and U= (    uu2   v2            )
!   (    .   .        )         (        .    .          )
!   (      .   .      )         (           .    .       )
!   (         .   .   )         (              .    .    )
!   (          lln   1)         (                     uun)
!

subroutine LU_routine(N,A,d,x)
implicit none
integer::N
real,dimension(1:N,1:N)::A !tridiagonal matrix
real,dimension(1:N)::d ! vector from the system A*x*=d
real,dimension(1:N)::x ! solution vector
real,dimension(1:N)::aa,b,c ! lower, main and upper diagonals from A
real,dimension(1:N)::ll !lower diagonal from L
real,dimension(1:N)::uu,v !main and upper diagonals from U
real,dimension(1:N)::z  !auxiliary vector from L*z=d and U*x=z
integer::i

!Store the diagonals from matrix A
!A-upper diagonal
do i=2,N
   c(i-1)=A(i-1,i)
end do
!A-main diagonal
do i=1,N
   b(i)=A(i,i)
end do
!A-lower diagonal
do i=1,N-1
   aa(i+1)=A(i+1,i)
end do


!L and U diagonals !see eq.(1.4.6)
do i=1,N-1
   v(i)=c(i)
end do
uu(1)=b(1)
do i=2,N
   ll(i)=aa(i)/uu(i-1)
   uu(i)=b(i)-ll(i)*v(i-1)
end do


! Solve L*z=d
!forward substitution
z(1)=d(1)
do i=2,N
   z(i)=d(i)-z(i-1)*ll(i)
end do


!Solve U*x=z
!backward substitution
x(N)=z(N)/uu(N)
do i=N-1,1,-1
   x(i)=(z(i)-v(i)*x(i+1))/uu(i)
end do

end subroutine LU_routine
