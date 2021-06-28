#code name: LU_scheme_short.py
#
#description: This code solves the tridiagonal matrix
#             system A*x=d using LU decomposition
#            (short version-without explicit matrices)
#
#reference: An Introduction to the Numerical
#           Solution of Differential Equation
#           Section 1.4.3. Tridiagonal Matrices, p.31
#           (Revised edition)-Douglas Quinney (1987)
#
#date: June 28th, 2021
#author: Leonardo Brito, Physics PdD student
#Institute of Physics (IFUSP)
#University of Sâo Paulo (USP), Brazil
#**********************************************************
import numpy as np

#--------------------------------------------------------
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
    z[0]=d[0]    
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

#main code
#problem A*x=d
#define the vector d and matrix A
N=4
d=np.array([4.3, 3.8, 3.1, 4.9])
A=np.array([(5.0,-1.0,0.0,0.0),(-1.0, 5.0, -1.0, 0.0),
(0.0,-1.0,5.0, -1.0),(0.0,0.0,-1.0, 5.0)])
x=np.array(N)
x=LU_scheme(N,A,d)

print('the solution is: ',x)
