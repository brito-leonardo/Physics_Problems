#code name: cn_heat_eq_1d.py
#
#description: This code solves the Heat equation by
#             Crank-Nicolson (CN) scheme and compares
#             the results with the exact solution
#
#reference:    1) Prof.Dr.Arnaldo Gammal Notes
#
#              2) An Introduction to the Numerical
#              Solution of Differential Equation
#              Section 1.4.3. Tridiagonal Matrices
#              (Revised edition)-Douglas Quinney (1987)
#
#date: June 29th, 2021
#author: Leonardo Brito, Physics PdD student
#Institute of Physics (IFUSP)
#University of Sâo Paulo (USP), Brazil
#**********************************************************
import numpy as np

#--------------------------------------------------------
#function: LU scheme 

def LU_scheme(N0,A0,d0):
    #we have the problme A*x0=d0
    #get the upper, main and lower diagonals
    c=np.zeros((N0))
    b=np.zeros((N0))
    a=np.zeros((N0))
    for i in range(1,N0):
        c[i-1]=A0[i-1,i]
    for i in range(0,N0):
        b[i]=A0[i,i]
    for i in range(0,N0-1):
        a[i+1]=A0[i+1,i]
    
    #get the L and U diagonals
    # we get A=L*U, L*U*x0=d0
    l=np.zeros((N0))
    u=np.zeros((N0))
    v=np.zeros((N0))
    for i in range(0,N0-1):
        v[i]=c[i]
    u[0]=b[0]    
    for i in range(1,N0):
        l[i]=a[i]/u[i-1]
        u[i]=b[i]-l[i]*v[i-1]
    #solve L*z=d0 by forward substitution    
    z=np.zeros((N0))
    z[0]=d0[0]    
    for i in range(1,N0):
        z[i]=d0[i]-z[i-1]*l[i]
    #get the solution x0 from 'U*x0=z'
    #by backward substitution
    x0=np.zeros((N0))
    x0[N0-1]=z[N0-1]/u[N0-1]
    for i in range(N0-2,-1,-1):
        x0[i]=(z[i]-v[i]*x0[i+1])/u[i]
    return x0

#--------------------------------------------------------
#function:  Crank-Nicolson scheme

def CN_scheme(Nx0,f0,sig,dx0,dt0):
    #set CN parameter
    rx=(dt0*sig)/(dx0**2) 
    #set the boundary conditions
    f0[0]=0.0
    f0[Nx0-1]=0.0
    #define the left side
    #define a reduced matrix A
    A=np.zeros((Nx0-2,Nx0-2))
    d0=np.zeros((Nx0))
    d=np.zeros((Nx0-2))
    #declare the upper diagonal
    for i in range(0,Nx0-3):
        A[i,i+1]=-0.5*rx
    #declare the main diagonal
    for i in range(0,Nx0-2):
        A[i,i]=1.0+rx
    #declare the lower diagonal 
    for i in range(1,Nx0-2):
        A[i,i-1]=-0.5*rx
    #declare right side
    d0[0]=0.0
    d0[Nx0-1]=0.0
    for i in range(1,Nx0-1):
        d0[i]=f0[i]+0.5*rx*(f0[i+1]-2*f0[i]+f0[i-1])
    #inser the result from 1 to Nx-2 into a new vector
    for i in range(0,Nx0-2):
        d[i]=d0[i+1]
    #call LU scheme and solve the problem A*f0r=d
    #f0r is f0 without the first and the last points 0 and Nx0-1   
    f0r=np.zeros(Nx0-2)
    f0r=LU_scheme(Nx0-2,A,d)
    #update f0
    for i in range(1,Nx0-1):
        f0[i]=f0r[i-1]
    #set the boundary conditions again
    f0[0]=0.0
    f0[Nx0-1]=0.0
    return f0

#--------------------------------------------------------
#function: evolution scheme

def evolution(Nx0,f0,sig,Nt0,dx0,dt0):
    f=f0
    for k in range(0,Nt0+1):
        f=CN_scheme(Nx0,f,sig,dx0,dt0)     
    return f    

#__________________________________________________________

#main code

Nx=201 #dimension of the vectors
Lx=1.0
dx=Lx/(Nx-1)
x=np.arange(0.0,Lx+dx,dx)
psi0=np.sin(np.pi*x) # initial guess/anstaz
psi=psi0 #enter the ansatz
#set the boundary conditions
psi[0]=0.0
psi[Nx-1]=0.0
# set the heat equation parameter
sig=1.0
#evolution parameters
Nt=40 #number of time steps
dt=0.0005 #time-step
psi=evolution(Nx,psi,sig,Nt,dx,dt) #time evolution 
tf=Nt*dt
#exact solution
y=np.sin(np.pi*x)*np.exp(-tf*(np.pi)**2)

#----------------------------------------------------------

#write the results in data files

#cn approximation
with open('heat_cn_py.txt','w') as data:
    for i in range(0,Nx-1):
        line=str(x[i])+' '+str(psi[i])
        line=line+'\n'
        data.writelines(line)
#exact result
with open('heat_exact_py.txt','w') as data2:
    for i in range(0,Nx-1):
        line2=str(x[i])+' '+str(y[i])
        line2=line2+'\n'
        data2.writelines(line2)
#__________________________________________________________
