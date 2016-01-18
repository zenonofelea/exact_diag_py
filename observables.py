from ED_python_public.spins1D import Hamiltonian1D
from ED_python_public.Basis import Basis1D
from numpy import *
import numpy as np
import scipy as sp

L=8


'''
Q=[[1.0j,i,(i+1)%L,(i+2)%L] for i in xrange(L)]
Qcc=[[-1.0j,i,(i+1)%L,(i+2)%L] for i in xrange(L)]
U=[[1.0,i,(i+1)%L] for i in xrange(L)]
h=[[1.0,i] for i in xrange(L)]
'''

PBC = 0

U = 1.0
J0 = 1
Q=0;
h=0;

zeta = 0.6; #driving amplitude in rot frame
dzeta = 0.2*zeta; #driving amplitude in rot frame
# infinite-frequency effective hoppings
J      = J0*2/pi*(zeta-dzeta)*cos(pi/2*(zeta-dzeta))/(1 - (zeta-dzeta)**2)
Jprime = J0*2/pi*(zeta+dzeta)*cos(pi/2*(zeta+dzeta))/(1 - (zeta+dzeta)**2)

betavec = [1, 1/1000]

#define parameter strings over the lattice

hstr = [[h,i] for i in xrange(L)]
hstaggstr = [[h*(-1)**i,i] for i in xrange(L)]

if PBC==1:

	J0str   = [[J0,i,(i+1)%L] for i in xrange(L)]
	J0ccstr = [[J0.conjugate(),i,(i+1)%L] for i in xrange(L)]

	Ustr    = [[U,i,(i+1)%L] for i in xrange(L)]

	Qstr    = [[Q,i,(i+1)%L,(i+2)%L] for i in xrange(L)]
	Qccstr  = [[Q.conjugate(),i,(i+1)%L,(i+2)%L] for i in xrange(L)]

elif PBC==0:

	Jstr        = [[-J,i,(i+1)%L] for i in xrange(0,L-1,2)]
	Jprimestr   = [[-Jprime,i,(i+1)%L] for i in xrange(1,L-1,2)]

	J0str   = [[J0,i,(i+1)%L] for i in xrange(L-1)]
	J0ccstr = [[J0.conjugate(),i,(i+1)%L] for i in xrange(L-1)]

	Ustr    = [[U,i,(i+1)%L] for i in xrange(L-1)]

	Qstr    = [[Q,i,(i+1)%L,(i+2)%L] for i in xrange(L-1)]
	Qccstr  = [[Q.conjugate(),i,(i+1)%L,(i+2)%L] for i in xrange(L-1)]


#####################################################################

#staticSSH=[['+-',Jstr],['-+',Jstr],['+-',Jprimestr],['-+',Jprimestr], ['zz',Ustr]] 

#static=[['z',hstr],['+-',J0str],['-+',J0ccstr],['zz',Ustr],['+z-',Qstr],['-z+',Qccstr]]
static1=[['+-',J0str],['-+',J0ccstr], ['zz',Ustr]]  
static2=[['+-',J0str],['-+',J0ccstr],['+z-',Qstr],['-z+',Qccstr]] 
dynamic=[]; #[['x',hstr,A]]

#####################################################################

#define symmetry parameters
symm= [0,1,1]


basis=Basis1D(L)#,Nup=L/2,kblock=symm[0],pblock=symm[1],zblock=symm[2])

#basis=Basis1D(L,Nup=L/2,kblock=0,pblock=+1,zblock=-1)
#basis=Basis1D(L,Nup=L/2,kblock=0,pblock=+1,zblock=+1)
#basis=Basis1D(L,Nup=L/2,kblock=0,pblock=+1)
#basis=Basis1D(L,Nup=L/2,kblock=0)

#H_SSH=Hamiltonian1D(staticSSH,dynamic,L,basis=basis)
H1=Hamiltonian1D(static1,dynamic,L,basis=basis)
H2=Hamiltonian1D(static2,dynamic,L,basis=basis)

#print 'hermiticity error is' ,np.linalg.norm( H_1.todense() - H_1.todense().conjugate().transpose() )

#print H1.todense()
#ESSH,VSSH = H_SSH.DenseEV(0)
E1,V1 = H1.DenseEV(0)
E2,V2 = H2.DenseEV(0)
#print 'evalues:', E1
#print '# e''values =', [len(E1), len(E2)]
#print '# evalues:', len(E)

"""
# calculate GS
psiSSHGS = VSSH[:,0]
ESSHGS = ESSH[0]
print ESSHGS/L
"""

Ns = basis.Ns

#print psiGS


def S_ent(L_A,L,psi,symm):
	# L_A: length of subsystem A
	# psi: pure quantum state
	# symm=[kblock,pblock,zblock]: vector with integers to specify symmetries of basis

	Ns_A = 2**L_A
	Ns_B = 2**(L-L_A)

	v = np.reshape(psi, (Ns_A, Ns_B))
	
	alpha = sp.linalg.svd(v, compute_uv=False, overwrite_a=True, check_finite=True)

	return -1./L_A*( ( abs(alpha)**2).dot( 2*log( abs(alpha)  ) ) ).sum()

#print "entanglement entropy:", [0.5*S_ent(L/2,L,psiSSHGS,symm), log(2)]	


def Diag_Ens_Observables(V1,V2,E1,betavec,L,alpha):
	# V1, V2:  matrices with pre and post quench eigenbases
	# E1: vector of energies of pre-quench
	# betavec: vector of inverse temperatures
	# alpha: Renyi entropy parameter

	Ns = len(E1) # Hilbert space dimension

	#define thernal density matrix
	rho = zeros((Ns,len(betavec)),dtype=np.float64)
	for i in xrange(len(betavec)):
		rho[:,i] = exp(-betavec[i]*(E1-E1[0]))/sum( exp(-betavec[i]*(E1-E1[0]) ) )

	#probability amplitudes
	a_n = (V1.conjugate().transpose()).dot(V2);
	V1 = None
	V2 = None
	# transition rates matrix
	T_nm = real( np.multiply(a_n, a_n.conjugate()) )
	# probability rates matrix (H1->H2->H1)
	pn = T_nm.dot(T_nm.transpose() )


	#calculate diagonal energy < H1 > in long time limit
	Ed_GS = (pn.transpose().dot(E1)/L )[0] # GS
	Ed_T  = (pn.transpose().dot(E1)/L ).dot(rho) # finite-temperature
	E_Tinf = E1.sum()/Ns/L # infinite temperature

	#calculate double-quench entropy (H1->H2->H1)
	Sdq_GS = (np.einsum( 'ij,ji->i', -pn.transpose(), log(pn) )/L )[0] # GS
	Sdq_T  = (np.einsum( 'ij,ji->i', -pn.transpose(), log(pn) )/L ).dot(rho) # finite-temperature
	S_Tinf = log(2); # infinite temperature

	# free up memory
	pn = None

	#calculate diagonal (Shannon) entropy (H1->H2)
	Sd_GS = (np.einsum( 'ij,ji->i', -T_nm.transpose(), log(T_nm) )/L )[0] # GS
	Sd_T  = (np.einsum( 'ij,ji->i', -T_nm.transpose(), log(T_nm) )/L ).dot(rho) # finite-temperature

	if alpha != 1:
		#calculate diagonal (Renyi) entropy for parameter alpha (H1->H2)
		Sd_Renyi_GS = (1/(1-alpha)*np.diag( log( np.power( T_nm, alpha )   ) )/L )[0] # GS
		Sd_Renyi_T = 1/(1-alpha)*( (log( np.power( np.diag(T_nm), alpha )  ) )/L  ).dot(rho) # finite-temperature
	else:
		print "Renyi parameter alpha cannot be 1!"
		Sd_Renyi_GS = "nan"
		Sd_Renyi_T = "nan"

	# calculate long-time energy fluctuations

	# <H1^2>
	H1_mn2 = (a_n.conjugate().transpose().dot( np.einsum('i,ij->ij',E1,a_n)) )**2
	a_n = None
	np.fill_diagonal(H1_mn2,0) 

	dE = np.real(np.einsum( 'ij,ji->i', T_nm, H1_mn2.dot(T_nm.transpose()) )/(L**2) )
	# free up memory
	T_nmF = None
	H1_mn2 = None

	deltaE_GS = dE[0] # GS
	deltaE_T  = dE.dot(rho) # finite-temperature
	# free up memory
	dE = None

	return [[Ed_GS, E_Tinf, Ed_T], [Sdq_GS, S_Tinf, Sdq_T], [Sd_GS, Sd_T], [Sd_Renyi_GS, Sd_Renyi_T], [deltaE_GS, deltaE_T]]

#alpha = 2	
#print Diag_Observables(V1,V2,E1,betavec,L,alpha)
	

def Kullback_Leibler_div(p1,p2):
	# p1,p2: probability distributions
	return np.multiply( p1, log( np.divide(p1,p2) ) ).sum()


# calculate energy as a function of time after a quench
def Observable_vs_time(psi,V2,E2,Obs,times):
	# psi: initial state
	# V2, E2: matrix w/ eigenbasis and vector of eogenvalues of post-quench Hamiltonian H2
	# Obs: observable of interest
	# times: vector with time values

	# project initial state onto basis V2
	c_n = V2.conjugate().transpose().dot(psi)


	# define time-evolved state in basis V2
	def psit(a_n,t):
		# a_n: probability amplitudes
		# t: time vector

		return V2.dot( np.multiply(exp(-1j*E2*t), a_n ) )

	
	Lt = len(times)

	# preallocate state
	# psi_time = zeros((Ns,Lt),dtype=np.complex128)

	# preallocate expectation value
	Expt_time = zeros((Lt),dtype=np.float128)

	# loop over time vector
	for m in xrange(Lt):
		# psi_time[:,m] = psit(c_nF,times[m])
		Expt_time[m] = real( psit(c_n,times[m]).conjugate().transpose().dot( Obs.dot(psit(c_n,times[m]),0) )  )

	return Expt_time

#times = [0,1,2,3,4,5.0]
#print Observable_vs_time(V1[:,0],V2,E2,H1,times)

def Mean_Level_Spacing(E):
	# compute consecuive E-differences
	sn = np.diff(E)
	# check for degeneracies
	if len(np.unique(E)) != len(E):
		print "Warning: degeneracies in the spectrum"
	# calculate the ratios of consecutive spacings
	aux = np.zeros((len(E)-1,2),dtype=np.float64)

	aux[:,0] = sn
	aux[:,1] = np.roll(sn,-1)

	return np.mean( np.divide( aux.min(1), aux.max(1) )[0:-1] )



