"""
Created on Fri Nov 27 2020

@author: Josh Peterson
"""

import numpy as np

t1 = np.array([[0.0+0.0j,0.5+0.0j,0.0+0.0j],[0.5+0.0j,0.0+0.0j,0.0+0.0j],[0.0+0.0j,0.0+0.0j,0.0+0.0j]])
t2 = np.array([[0.0+0.0j,0.0-0.5j,0.0+0.0j],[0.0+0.5j,0.0+0.0j,0.0+0.0j],[0.0+0.0j,0.0+0.0j,0.0+0.0j]])
t3 = np.array([[0.5+0.0j,0.0+0.0j,0.0+0.0j],[0.0+0.0j,-0.5+0.0j,0.0+0.0j],[0.0+0.0j,0.0+0.0j,0.0+0.0j]])
t4 = np.array([[0.0+0.0j,0.0+0.0j,0.5+0.0j],[0.0+0.0j,0.0+0.0j,0.0+0.0j],[0.5+0.0j,0.0+0.0j,0.0+0.0j]])
t5 = np.array([[0.0+0.0j,0.0+0.0j,0.0-0.5j],[0.0+0.0j,0.0+0.0j,0.0+0.0j],[0.0+0.5j,0.0+0.0j,0.0+0.0j]])
t6 = np.array([[0.0+0.0j,0.0+0.0j,0.0+0.0j],[0.0+0.0j,0.0+0.0j,0.5+0.0j],[0.0+0.0j,0.5+0.0j,0.0+0.0j]])
t7 = np.array([[0.0+0.0j,0.0+0.0j,0.0+0.0j],[0.0+0.0j,0.0+0.0j,0.0-0.5j],[0.0+0.0j,0.0+0.5j,0.0+0.0j]])
t8 = np.array([[(1/(2*np.sqrt(3)))+0.0j,0.0+0.0j,0.0+0.0j],[0.0+0.0j,(1/(2*np.sqrt(3)))+0.0j,0.0+0.0j],[0.0+0.0j,0.0+0.0j,-(1/np.sqrt(3))+0.0j]])

basis = [t1, t2, t3, t4, t5, t6, t7, t8]

def commute(ti, tj):
    ij = np.matmul(ti,tj)
    ji = np.matmul(tj,ti)
    return np.add(ij,-1*ji)

def CommuteBasisMatrices():
    for i in range(len(basis)):
        for j in range(len(basis)):
            print("T(",i+1,") into T(",j+1,"):")
            print(commute(basis[i], basis[j]))
         
def Computef(M):
    f8 = M[2][2].imag / t8[2][2].real
    f6 = M[2][1].imag / t6[2][1].real
    f7 = -M[2][1].real / t7[2][1].imag
    f5 = -M[2][0].real / t5[2][0].imag
    f4 = M[2][0].imag / t4[2][0].real
    f2 =  -M[1][0].real / t2[1][0].imag
    f1 = M[1][0].imag / t1[1][0].real
    f3 = (M[1][1].imag - f8*t8[1][1].real) / t3[1][1].real
    return [f1,f2,f3,f4,f5,f6,f7,f8]
    
def ComputeAllf():
    for i in range(len(basis)):
        for j in range(len(basis)):
            f = Computef(commute(basis[i],basis[j]))
            print("f(",i+1,j+1,") =",f)

"""
Part b.
_______________________________________________________________________________
The function below evaluates all f.  The output is formated as f(a,b) = [f(abc)] for
all values of c.  Looking through the list clearly shows that f(abc) is totally
antisymmetric.
"""
print("Part b:")
ComputeAllf()
print("")

def ComputeTraces():
    for i in range(len(basis)):
        for j in range(len(basis)):
            print("Trace of T(",i+1,")T(",j+1,") =", np.trace(np.matmul(basis[i],basis[j])))
            
"""
Part c.
_______________________________________________________________________________
The function below evaluates the traces of t(a)t(b) for all a and b.  Looking
at the output, it is clear that condition 15.78 holds, with C(r) = 0.5.
"""
print("Part c:")
ComputeTraces()
print("")

def ComputeC2():
    c2 = 0
    for t in basis:
        c2 += np.matmul(t,t)
    print("LHS =",c2)

"""
Part d.
_______________________________________________________________________________
The function below computes the LHS of equation 15.92.  From the output it is
appearent that C2(r) = 4/3.
"""      
print("Part d:")
ComputeC2()
print("")