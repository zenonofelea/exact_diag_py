from ED_python_public.spins1D import Hamiltonian1D
from ED_python_public.Basis import Basis1D
from numpy import *

L=10

def A(t):
	return t

Q=[[1.0j,i,(i+1)%L,(i+2)%L] for i in xrange(L)]
Qcc=[[1.0j,i,(i+1)%L,(i+2)%L] for i in xrange(L)]
J=[[1.0,i,(i+1)%L] for i in xrange(L)]
h=[[1.0,i] for i in xrange(L)]

static=[['zz',J],['x',h],['+z-',Q],['-z+',Qcc]] 
dynamic=[['x',h,A]]


basis=Basis1D(L,Nup=L/2,kblock=0)
H_1=Hamiltonian1D(static,dynamic,L,basis=basis)

print H_1.DenseEE()






