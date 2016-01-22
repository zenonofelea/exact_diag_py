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
J0 = sqrt(3)
Q=0;
h=0;

zeta = 0.6; #driving amplitude in rot frame
dzeta = 0.2*zeta; #driving amplitude in rot frame
# infinite-frequency effective hoppings
J      = J0*2/pi*(zeta-dzeta)*cos(pi/2*(zeta-dzeta))/(1 - (zeta-dzeta)**2)
Jprime = J0*2/pi*(zeta+dzeta)*cos(pi/2*(zeta+dzeta))/(1 - (zeta+dzeta)**2)

betavec = [1.0, 1.0/100]

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
#static1=[ ['zz',Ustr]]
static1=[['+-',J0str],['-+',J0ccstr], ['zz',Ustr]]    
static2=[['+-',J0str],['-+',J0ccstr],['+z-',Qstr],['-z+',Qccstr]] 
dynamic=[]; #[['x',hstr,A]]

#####################################################################

# use for checking S_ent
basis=Basis1D(L)
# use for checking the rest
basis=Basis1D(L,Nup=L/2,kblock=1,pblock=1,zblock=1)


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


# calculate GS
psi1GS = V1[:,0]
E1GS = E1[0]
print E1GS/L


Ns = basis.Ns

#print psiGS


def Renyi_entropy(L_A,L,psi, init_site=0,alpha=1):
	# L_A: length of subsystem A
	# psi: pure quantum state
	# alpha: Renyi parameter
	# init_site: initial site counted from which subsystem A spans L_A sites to the right
	
	if alpha < 0:
		print "alpha must be a nonnegative real number"

	if init_site > L-1:
		print "init_site cannot exceed lattice site number"	
	

	#calculate H-space dimensions of the subsystem and the system
	Ns_A = 2**L_A

	if init_site == 0:
		Ns_Ac = 2**(L-L_A) # dim of H-space of complement of A
		v = np.real( np.reshape(psi, (Ns_A, Ns_Ac)) ) #cast initial state into a vector of the tensor product space
	else:
		Ns_Ac_l = 2**(init_site) # dim of H-space of complement left of A
		Ns_Ac_r = 2**(L-L_A - init_site)  # dim of H-space of complement right of A
		v = np.real( np.reshape(psi, (Ns_A, Ns_Ac_l*Ns_Ac_r)) ) #cast initial state into a vector of the tensor product space
		v = np.roll(v,Ns_Ac_l, axis=1) # permute columns to put the dimension of subsystem A at the right place

	# apply singular value decomposition
	gamma = sp.linalg.svd(v, compute_uv=False, overwrite_a=True, check_finite=True)

	# calculate Renyi entropy
	if any(gamma == 1.0):
		return 0
	else:
		if alpha == 1.0:
			return -1./L_A*( ( abs(gamma)**2).dot( 2*log( abs(gamma)  ) ) ).sum()
		else:
			return  1./L_A*( 1./(1-alpha)*log( (gamma**alpha).sum() )  )

#S_ent = Renyi_entropy(L/2,L,psi1GS,alpha=3,init_site=2)
#print "entanglement entropy:", [0.5*S_ent, log(2)]	


def Diag_Ens_Observables(L,V1,V2,E1,betavec=[],alpha=1, Obs=False, Ed=False,S_double_quench=False,Sd_Renyi=False,deltaE=False):
	# V1, V2:  matrices with pre and post quench eigenbases
	# E1: vector of energies of pre-quench
	# Obs: any hermitian observable
	# betavec: vector of inverse temperatures
	# alpha: Renyi entropy parameter

	variables_GS = []
	variables_T = []

	if any(Obs):
		variables_GS.append("Obs_GS")
		variables_T.append("Obs_T")
	if Ed:
		variables_GS.append("Ed_GS")
		variables_GS.append("E_Tinf")
		variables_T.append("Ed_T")
	if S_double_quench:
		variables_GS.append("S_double_quench_GS")
		variables_T.append("S_double_quench_T")
	if Sd_Renyi:
		variables_GS.append("Sd_Renyi_GS")
		variables_T.append("Sd_Renyi_T")
	if S_double_quench or Sd_Renyi:
		variables_GS.append("S_Tinf")
	if deltaE:
		variables_GS.append("deltaE_GS")
		variables_T.append("deltaE_T")

	if not variables_GS:
		print "No observables were requested: ..exiting"
		return None

	
	print "WARNING: if no symmetries are used, the probabilities pn have 0 entries, and log(0) gives an error"
	print "I tested the code with symmetries"

	Ns = len(E1) # Hilbert space dimension

	if betavec:
		#define thermal density matrix w.r.t. the basis V1	
		rho = zeros((Ns,len(betavec)),dtype=np.float64)
		for i in xrange(len(betavec)):
			rho[:,i] = exp(-betavec[i]*(E1-E1[0]))/sum( exp(-betavec[i]*(E1-E1[0]) ) )

	# diagonam matrix elements of Obs in the basis V2
	if any(Obs):
		O_mm = np.real( np.einsum( 'ij,ji->i', V2.transpose().conj(), Obs.dot(V2 ) ) )
	#probability amplitudes
	a_n = (V1.conjugate().transpose()).dot(V2);
	V1 = None
	V2 = None
	# transition rates matrix (H1->H2)
	T_nm = real( np.multiply(a_n, a_n.conjugate()) )
	# probability rates matrix (H1->H2->H1)
	if Ed or S_double_quench:
		pn = T_nm.dot(T_nm.transpose() )

	# diagonal ens expectation value of Obs in post-quench basis
	if any(Obs):
		Obs_GS = (np.einsum( 'ij,j->i', T_nm.transpose(), O_mm )/L )[0] # GS
		if betavec:
			Obs_T = (np.einsum( 'ij,j->i', T_nm.transpose(), O_mm )/L ).dot(rho) # finite-temperature


	#calculate diagonal energy <H1> in long time limit
	if Ed:
		Ed_GS = (pn.transpose().dot(E1)/L )[0] # GS
		if betavec:
			Ed_T  = (pn.transpose().dot(E1)/L ).dot(rho) # finite-temperature
		E_Tinf = E1.sum()/Ns/L # infinite temperature

	#calculate double-quench entropy (H1->H2->H1)
	if S_double_quench:
		S_double_quench_GS = (np.einsum( 'ij,ji->i', -pn.transpose(), log(pn) )/L )[0] # GS
		if betavec:
			S_double_quench_T  = (np.einsum( 'ij,ji->i', -pn.transpose(), log(pn) )/L ).dot(rho) # finite-temperature
	

	# free up memory
	pn = None

	#calculate diagonal Renyi entropy for parameter alpha: equals (Shannon) entropy for alpha=1: (H1->H2)
	if Sd_Renyi:
		if alpha != 1.0:
			#calculate diagonal (Renyi) entropy for parameter alpha (H1->H2)
			Sd_Renyi_GS = (1/(1-alpha)*np.diag( log( np.power( T_nm, alpha )   ) )/L )[0] # GS
			if betavec:
				Sd_Renyi_T = 1/(1-alpha)*( (log( np.power( np.diag(T_nm), alpha )  ) )/L  ).dot(rho) # finite-temperature
		else:
			print "Renyi entropy equals diagonal entropy"
			Sd_Renyi_GS = (np.einsum( 'ij,ji->i', -T_nm.transpose(), log(T_nm) )/L )[0] # GS
			if betavec:
				Sd_Renyi_T = (np.einsum( 'ij,ji->i', -T_nm.transpose(), log(T_nm) )/L ).dot(rho) # finite-temperature

	# infinite temperature entropy
	if S_double_quench or Sd_Renyi:
		S_Tinf = log(2); 

	# calculate long-time energy fluctuations
	if deltaE:
		# calculate <H1^2>
		H1_mn2 = (a_n.conjugate().transpose().dot( np.einsum('i,ij->ij',E1,a_n)) )**2
		a_n = None
		np.fill_diagonal(H1_mn2,0) 

		dE = np.real(np.einsum( 'ij,ji->i', T_nm, H1_mn2.dot(T_nm.transpose()) )/(L**2) )
		# free up memory
		T_nmF = None
		H1_mn2 = None

		deltaE_GS = dE[0] # GS
		if betavec:
			deltaE_T  = dE.dot(rho) # finite-temperature
		# free up memory
		dE = None

	return_dict = {}
	for i in range(len(variables_GS)):
		return_dict[variables_GS[i]] = vars()[variables_GS[i]]
	if betavec:
		return_dict_T = {}
		for i in range(len(variables_T)):
			return_dict_T[variables_T[i]] = vars()[variables_T[i]]
    	# merge dictionaries
    	if betavec: # no idea why it wants this if clause, I thought the previous if betavec one should be enough
    		return_dict.update(return_dict_T)
	
	return return_dict
	#return [[Ed_GS, E_Tinf, Ed_T], [S_double_quench_GS, S_Tinf, S_double_quench_T], [Sd_Renyi_GS, Sd_Renyi_T], [deltaE_GS, deltaE_T]]


print Diag_Ens_Observables(L,V1,V2,E1,Obs = V1+V1.transpose().conj(), Ed = True, betavec=betavec	)

	

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
	# compute consecutive E-differences
	sn = np.diff(E)
	# check for degeneracies
	if len(np.unique(E)) != len(E):
		raise Exception("Warning: degeneracies in the spectrum")
	# calculate the ratios of consecutive spacings
	aux = np.zeros((len(E)-1,2),dtype=np.float64)

	aux[:,0] = sn
	aux[:,1] = np.roll(sn,-1)

	return np.mean( np.divide( aux.min(1), aux.max(1) )[0:-1] )



