#from ED_python_public.spins1D import Hamiltonian1D
#from ED_python_public.Basis import Basis1D
from exact_diag_py.hamiltonian import hamiltonian
from exact_diag_py.basis import basis1d
from numpy import *
import numpy as np
import scipy as sp

import sys
import os



L=int(sys.argv[1])


PBC = 0

U = 1.0
J0 = 1.0
Q=sqrt(3);
h=0*13.33;

zeta = 0.6; #driving amplitude in rot frame
dzeta = 0.2*zeta; #driving amplitude in rot frame
# infinite-frequency effective hoppings
J      = J0*2/pi*(zeta-dzeta)*cos(pi/2*(zeta-dzeta))/(1 - (zeta-dzeta)**2)
Jprime = J0*2/pi*(zeta+dzeta)*cos(pi/2*(zeta+dzeta))/(1 - (zeta+dzeta)**2)

betavec = [1.0, 1.0/100]

#define parameter strings over the lattice

hstr = [[h,L-1]] #[[h,i] for i in xrange(L)]
hstaggstr = [[h*(-1)**i,i] for i in xrange(L)]

if PBC==1:

	Jstr        = [[-J,i,(i+1)%L] for i in xrange(0,L,2)]
	Jprimestr   = [[-Jprime,i,(i+1)%L] for i in xrange(1,L,2)]

	J0str   = [[J0,i,(i+1)%L] for i in xrange(L)]
	J0ccstr = [[J0.conjugate(),i,(i+1)%L] for i in xrange(L)]

	Ustr    = [[U,i,(i+1)%L] for i in xrange(L)]

	Qstr    = [[Q,i,(i+1)%L] for i in xrange(L)]

	#basis=basis1d(L,Nup=L/2,kblock=0,pblock=1,zblock=1)
	basis=basis1d(L)

elif PBC==0:

	Jstr        = [[-J,i,(i+1)%L] for i in xrange(0,L-1,2)]
	Jprimestr   = [[-Jprime,i,(i+1)%L] for i in xrange(1,L-1,2)]

	J0str   = [[J0,i,(i+1)%L] for i in xrange(L-1)]
	J0ccstr = [[J0.conjugate(),i,(i+1)%L] for i in xrange(L-1)]

	Ustr    = [[U,i,(i+1)%L] for i in xrange(L-1)]

	Qstr    = [[Q,i,(i+1)%L] for i in xrange(L-1)]

	#basis=basis1d(L,Nup=L/2,pblock=1,zblock=1)
	basis=basis1d(L)

#####################################################################

staticSSH=[['+-',Jstr],['-+',Jstr],['+-',Jprimestr],['-+',Jprimestr], ['zz',Ustr],['x',hstr]] 

#static=[['z',hstr],['+-',J0str],['-+',J0ccstr],['zz',Ustr],['+z-',Qstr],['-z+',Qccstr]]
#static1=[ ['zz',Ustr]]
static1=[['+-',J0str],['-+',J0ccstr],['zz',Ustr]]    
static2=[['+-',J0str],['-+',J0ccstr],['zz',Qstr]] 
dynamic=[]; #[['x',hstr,A]]

Ns = basis.Ns

#####################################################################

# use for checking S_ent
#basis=Basis1D(L)
# use for checking the rest
#basis=Basis1D(L,Nup=L/2,pblock=1)
#basis=Basis1D(L,Nup=L/2,zblock=1)



H_SSH=hamiltonian(staticSSH,dynamic,L,basis=basis,pauli=False)
H1=hamiltonian(static1,dynamic,L,basis=basis)
H2=hamiltonian(static2,dynamic,L,basis=basis)

#print 'hermiticity error is' ,np.linalg.norm( H_1.todense() - H_1.todense().conjugate().transpose() )

#print H1.todense()
ESSH,VSSH = H_SSH.eigh(overwrite_a=True, overwrite_b=True) 
E1,V1 = H1.eigh(overwrite_a=True, overwrite_b=True) 
E2,V2 = H2.eigh(overwrite_a=True, overwrite_b=True) 

# calculate GS
psiGS = np.real( VSSH[:,0] )
E1GS = ESSH[0]

#print psi1GS
print E1GS/L



def Entanglement_entropy(L,psi,subsys=[i for i in xrange( int(L/2) )],basis=None,alpha=1.0, DM=False):
	# psi: pure quantum state
	# subsys: a list of integers modelling the site numbers of the subsystem
	# basis: the basis of the Hamiltonian: needed only when symmetries are used
	# alpha: Renyi parameter
	# DM: if on returns the reduced density matrix corresponding to psi

	variables = ["Sent"]

	if np.any(DM):
		variables.append("DM")


	# turn the following messages into WARNINGS	
	
	if alpha < 0.0:
		print "alpha must be a nonnegative real number"

	if len(subsys) > np.floor(L/2):
		print "subsystem length cannot exceed half the total lattice site number"

	if max(subsys) > L-1:
		print "subsystem definition contains sites exceeding the total lattice site number"

	
	# re-write the state in the initial basis
	if len(psi)<2**L:
		print "need to parse the basis variable"
		psi = np.asarray( basis.get_vec(psi,sparse=True).todense().T )[0,:]
	del basis

	#calculate H-space dimensions of the subsystem and the system
	L_A = len(subsys)
	Ns_A = 2**L_A
	
	# define lattice indices putting the subsystem to the left
	for i in xrange(L):
		if i in subsys:
			continue
		else:
			subsys.append(i)

	'''
	the algorithm for the entanglement entropy of an arbitrary subsystem goes as follows:

	1) the initial state psi has 2^L entries corresponding to the spin-z configs
	2) reshape psi into a 2x2x2x2x...x2 dimensional array (L products in total). Call this array v.
	3) v should satisfy the property that v[0,1,0,0,0,1,...,1,0], total of L entries, should give the entry of psi 
	   corresponding to the spin-z basis vector (0,1,0,0,0,1,...,1,0). This ensures a correspondence of the v-indices
	   (and thus the psi-entries) to the L lattice sites.
	4) fix the lattice sites that define the subsystem L_A, and reshuffle the array v according to this: e.g. if the 
 	   subsystem consistes of sites (k,l) then v should be reshuffled such that v[(k,l), (all other sites)]
 	5) reshape v[(k,l), (all other sites)] into a 2D array of dimension ( L_A x L/L_A ) and proceed with the SVD as below  

	'''

	# performs 2) and 3)
	v = np.reshape(psi, tuple([2 for i in xrange(L)] ) )
	del psi
	# performs 4)
	v = np.transpose(v, axes=subsys) 
	# performs 5)
	v = np.reshape(v, ( Ns_A, Ns/Ns_A) )
	
	
	# apply singular value decomposition
	if DM==False:
		gamma = sp.linalg.svd(v, compute_uv=False, overwrite_a=True, check_finite=True)
	else:
		U, gamma, V = sp.linalg.svd(v, full_matrices=False, overwrite_a=True, check_finite=True)
		# calculate reduced density matrix DM	
		DM =   reduce( np.dot, [U, np.diag(gamma**2) , U.T.conjugate() ] )   #.dot( dot( np.einsum('i,ij->ij',gamma,Vh) ).dot(Vh)  )
		del U, V	
		print "NEED TO TEST this reduced DM against something known!!!"
	del v


	# calculate Renyi entropy
	if any(gamma == 1.0):
		Sent = 0.0
	else:
		if alpha == 1.0:
			Sent = -1./L_A*( ( abs(gamma)**2).dot( 2*log( abs(gamma)  ) ) ).sum()
		else:
			Sent =  1./L_A*( 1./(1-alpha)*log( (gamma**alpha).sum() )  )

	# define dictionary with outputs
	return_dict = {}
	for i in range(len(variables)):
		return_dict[variables[i]] = vars()[variables[i]]

	return return_dict

S_ent0 = Entanglement_entropy(L,psiGS,subsys=[0,1],alpha=1,basis=basis,DM=True)
print "entanglement entropy:", S_ent0['Sent']
print "DM:", S_ent0['DM']

br

def Diag_Ens_Observables(L,V1,V2,E1,betavec=[],alpha=1.0, Obs=False, Ed=False,S_double_quench=False,Sd_Renyi=False,deltaE=False):
	# V1, V2:  matrices with pre and post quench eigenbases
	# E1: vector of energies of pre-quench
	# Obs: any hermitian observable
	# betavec: vector of inverse temperatures
	# alpha: Renyi entropy parameter

	

	print "np.any(Obs) yells at me when Obs is sparse; any ideas?"

	variables_GS = []
	variables_T = []

	if np.any(Obs):
		variables_GS.append("Obs_GS")
		variables_T.append("Obs_T")
	if Ed:
		print "The value of E_Tinf depends on the symmetries used"
		variables_GS.append("Ed_GS")
		variables_GS.append("E_Tinf")
		variables_T.append("Ed_T")
		variables_T.append("E_Tave")
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

	
	Ns = len(E1) # Hilbert space dimension

	if betavec:
		print "All thermal expectation values depend statistically on the symmetry used via the available number of states as part of the system-size dependence"
		#define thermal density matrix w.r.t. the basis V1	
		rho = zeros((Ns,len(betavec)),dtype=np.float64)
		for i in xrange(len(betavec)):
			rho[:,i] = exp(-betavec[i]*(E1-E1[0]))/sum( exp(-betavec[i]*(E1-E1[0]) ) )

	# diagonal matrix elements of Obs in the basis V2
	if np.any(Obs):
		O_mm = np.real( np.einsum( 'ij,ji->i', V2.transpose().conj(), Obs.dot(V2) ) )
	#probability amplitudes
	a_n = V1.conjugate().transpose().dot(V2);
	V1 = None
	V2 = None
	# transition rates matrix (H1->H2)
	T_nm = np.real( np.multiply(a_n, a_n.conjugate()) )
	T_nm[T_nm<=1E-16] = np.finfo(float).eps	
	# probability rates matrix (H1->H2->H1)
	if Ed or S_double_quench:
		pn = T_nm.dot(T_nm.transpose() )


	# diagonal ens expectation value of Obs in post-quench basis
	if np.any(Obs):
		Obs_GS = T_nm[0,:].dot(O_mm)/L # GS
		if betavec:
			Obs_T = (np.einsum( 'ij,j->i', T_nm, O_mm )/L ).dot(rho) # finite-temperature


	#calculate diagonal energy <H1> in long time limit
	if Ed:
		Ed_GS = pn[0,:].dot(E1)/L  # GS
		if betavec:
			Ed_T  = (pn.dot(E1)/L ).dot(rho) # finite-temperature
			E_Tave = E1.dot(rho)/L # average energy density
		E_Tinf = E1.sum()/Ns/L # infinite temperature

	#calculate double-quench entropy (H1->H2->H1)
	if S_double_quench:
		S_double_quench_GS = -pn[0,:].dot(log(pn[0,:]))/L # GS
		if betavec:
			S_double_quench_T  = (np.einsum( 'ij,ji->i', -pn, log(pn) )/L ).dot(rho) # finite-temperature
	

	# free up memory
	pn = None

	#calculate diagonal Renyi entropy for parameter alpha: equals (Shannon) entropy for alpha=1: (H1->H2)
	if Sd_Renyi:
		if alpha != 1.0:
			#calculate diagonal (Renyi) entropy for parameter alpha (H1->H2)
			Sd_Renyi_GS = 1/(1-alpha)*log( np.power( T_nm[0,:], alpha ).sum() )/L  # # GS
			if betavec:
				Sd_Renyi_T = 1/(1-alpha)*( log( np.power( T_nm, alpha ).sum(1)  )/L  ).dot(rho) # finite-temperature
		else:
			print "Renyi entropy equals diagonal entropy"
			Sd_Renyi_GS = -T_nm[0,:].dot( log(T_nm[0,:]) ) /L # GS
			if betavec:
				Sd_Renyi_T = (np.einsum( 'ij,ji->i', -T_nm, log(T_nm.transpose()) )/L ).dot(rho) # finite-temperature

	# infinite temperature entropy
	if S_double_quench or Sd_Renyi:
		S_Tinf = log(2); 

	# calculate long-time energy fluctuations
	if deltaE:
		# calculate <H1^2>
		H1_mn2 = (a_n.conjugate().transpose().dot( np.einsum('i,ij->ij',E1,a_n)) )**2
		a_n = None
		np.fill_diagonal(H1_mn2,0) 

		deltaE_GS = np.real( reduce( np.dot,[T_nm[0,:], H1_mn2, T_nm[0,:] ])  )/L**2  # GS
		if betavec:
			deltaE_T  = np.real(np.einsum( 'ij,ji->i', T_nm, H1_mn2.dot(T_nm.transpose()) )/(L**2) ).dot(rho) # finite-temperature
		# free up memory
		T_nmF = None
		H1_mn2 = None

	return_dict = {}
	for i in range(len(variables_GS)):
		return_dict[variables_GS[i]] = vars()[variables_GS[i]]
	if betavec:
		for i in range(len(variables_T)):
			return_dict[variables_T[i]] = vars()[variables_T[i]]
			

	
	return return_dict
	#return [[Ed_GS, E_Tinf, Ed_T], [S_double_quench_GS, S_Tinf, S_double_quench_T], [Sd_Renyi_GS, Sd_Renyi_T], [deltaE_GS, deltaE_T]]


#print Diag_Ens_Observables(L,V1,V2,E1,Obs = V1+V1.transpose().conj(), Ed = True, betavec=betavec	)
print Diag_Ens_Observables(L,V1,V2,E1, Obs=H1.todense(),Ed=True, betavec=betavec	)

	

def Kullback_Leibler_div(p1,p2):
	# p1,p2: probability distributions
	return np.multiply( p1, log( np.divide(p1,p2) ) ).sum()


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
		Expt_time[m] = np.real( reduce( np.dot, [psit(c_n,times[m]).conjugate().T, Obs, psit(c_n,times[m]) ]  )  )

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



