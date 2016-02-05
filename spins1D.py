import sys # needed for sys.stop('error message')
#local modules:
from Basis import Basis1D
#from py_lapack import eigh # used to diagonalize hermitian and symmetric matricies

#python 2.7 modules
from scipy.linalg import norm, eigh
from scipy.sparse import coo_matrix	# needed as the initial format that the Hamiltonian matrices are stored as
from scipy.sparse import csr_matrix	# the final version the sparse matrices are stored as, good format for dot produces with vectors.
from scipy.sparse.linalg  import eigsh	# needed for the sparse linear algebra packages
from scipy.integrate import complex_ode,ode	# ode solver used in evolve wave function.
from numpy import isscalar,dot,asarray, array, int32, int64, float32, float64, complex64, complex128
from copy import deepcopy


#global names:
supported_dtypes=(float32, float64, complex64, complex128)



def StaticH(B,static,dtype):
	"""
	args:
		static=[[Sopstr_1,Popstr_1,indx_1],...,[Sopstr_n,Popstr_n,indx_n]], list of Sopstr,Popstr,indx to add up for static piece of Hamiltonian.
		dtype = the low level C-type which the matrix should store its values with.
	returns:
		H: a csr_matrix representation of the list static

	description:
		this function takes the list static and creates a list of matrix elements is coordinate format. it does
		this by calling the basis method Op which takes a state in the basis, acts with opstr and returns a matrix 
		element and the state which it is connected to. This function is called for ever opstr in static and for every 
		state in the basis until the entire hamiltonian is mapped out. It takes those matrix elements (which need not be 
		sorted or even unique) and creates a coo_matrix from the scipy.sparse library. It then converts this coo_matrix
		to a csr_matrix class which has optimal sparse matrix vector multiplication.
	"""
	ME_list=[] # this is a list which stores the matrix elements as lists [[row,col,ME],...] for the whole hamiltonian. 
	for i in xrange(len(static)): 
		List=static[i]
		Sopstr=List[0]
		Popstr=List[1]
		bonds=List[2]
		for bond in bonds:
			J=bond[0]
			indx=bond[1:]
			ME_list.extend(B.Op(J,Sopstr,Popstr,indx))

	if static: # if static is not an empty list []:
		# there is no way to tranpose a list so we must convert to array, this process will convert all parts of the list to the most compatible type.
		ME_list=asarray(ME_list).T.tolist() # transpose list so that it is now [[row,...],[col,...],[ME,...]] which is how coo_matrix is constructed.
		ME_list[1]=map( lambda a:int(abs(a)), ME_list[1]) # convert the indices back to integers 
		ME_list[2]=map( lambda a:int(abs(a)), ME_list[2])	# convert the indices back to integers
		#print "ME_list[0]", ME_list[0]
		#print "ME_list[1]", ME_list[1]
		#print "ME_list[2]", ME_list[2]
		H=coo_matrix((ME_list[0],(ME_list[1],ME_list[2])),shape=(B.Ns,B.Ns),dtype=dtype) # construct coo_matrix
		H=H.tocsr() # convert to csr_matrix
		H.sum_duplicates() # sum duplicate matrix elements
		H.eliminate_zeros() # remove all zero matrix elements
		return H 
	else: # else return None which indicates there is no static part of Hamiltonian.
		return None








def DynamicHs(B,dynamic,dtype):
	"""
	args:
	dynamic=[[Sopstr_1, Popstr_1,indx_1,func_1],...,[Sopstr_n,Popstr_n,indx_n,func_1]], list of Sopstr,Popstr,indx and functions to drive with
	dtype = the low level C-type which the matrix should store its values with.

	returns:
	tuple((func_1,H_1),...,(func_n,H_n))

	H_i: a csr_matrix representation of Sopstr_i,Popstr_i,indx_i
	func_i: callable function of time which is the drive term in front of H_i

	description:
		This function works the same as static, but instead of adding all of the elements 
		of the dynamic list together, it returns a tuple which contains each individual csr_matrix 
		representation of all the different driven parts. This way one can construct the time dependent 
		Hamiltonian simply by looping over the tuple returned by this function. 
	"""
	Dynamic_Hs=[]
	for i in xrange(len(dynamic)):
		ME_list=[]
		List=dynamic[i]
		Sopstr=List[0]
		Popstr=List[1]
		bonds=List[2]
		for bond in bonds:
			J=bond[0]
			indx=bond[1:]
			ME_list.extend(B.Op(J,Sopstr,Popstr,indx))
	
		ME_list=asarray(ME_list).T.tolist()
		ME_list[1]=map( lambda a:int(abs(a)), ME_list[1])
		ME_list[2]=map( lambda a:int(abs(a)), ME_list[2])
		H=coo_matrix((ME_list[0],(ME_list[1],ME_list[2])),shape=(B.Ns,B.Ns),dtype=dtype)
		H=H.tocsr()
		H.sum_duplicates()
		H.eliminate_zeros()
		Dynamic_Hs.append((List[3],H))

	return tuple(Dynamic_Hs)










class Hamiltonian1D:
	def __init__(self,static,dynamic,L,Nph,**init_params):
		"""
		This function intializes the Hamtilonian. You can either initialize with symmetries, or an instance of Basis1D.
		Note that if you initialize with a basis it will ignore all symmetry inputs.
		"""
		Nph=init_params.get("Nph")

		Nup=init_params.get("Nup")
		kblock=init_params.get("kblock")
		zblock=init_params.get("zblock")
		pblock=init_params.get("pblock")
		pzblock=init_params.get("pzblock")
		dtype=init_params.get("dtype")
		a=init_params.get("a")
		basis=init_params.get("basis")
		if a == None:
			a=1
		if dtype == None:
			dtype=complex128
		if basis == None:  
			basis=Basis1D(L,Nph=Nph,Nup=Nup,a=a,kblock=kblock,zblock=zblock,pblock=pblock,pzblock=pzblock)
		if not isinstance(basis,Basis1D):
			raise TypeError("basis is not instance of Basis1D")
		if dtype not in supported_dtypes:
			raise TypeError("Hamiltonian1D doesn't support type: "+str(dtype))


		self.static=static
		self.dynamic=dynamic
		self.L=L
		self.Nph=Nph
		self.Ns=basis.Ns
		self.dtype=dtype
		if self.Ns > 0:
			self.Static_H=StaticH(basis,static,dtype)
			self.Dynamic_Hs=DynamicHs(basis,dynamic,dtype)
			self.shape=(self.Ns,self.Ns)





	def tocsr(self,time=0):
		"""
		args:
			time=0, the time to evalute drive at.

		description:
			this function simply returns a copy of the Hamiltonian as a csr_matrix evaluated at the desired time.
		"""
		if self.Ns <= 0:
			return csr_matrix(asarray([[]]))
		if not isscalar(time):
			raise Exception("time must be a scaler")

		if self.Static_H != None: # if there is a static Hamiltonian...
			H=self.Static_H	
			for ele in self.Dynamic_Hs:
				H = H + ele[1]*ele[0](time)
		else: # if there isn't...
			for ele in self.Dynamic_Hs:
				H = H + ele[1]*ele[0](time)

		return H





	def todense(self,time=0):
		"""
		args:
			time=0, the time to evalute drive at.

		description:
			this function simply returns a copy of the Hamiltonian as a dense matrix evaluated at the desired time.
			This function can overflow memory if not careful.
		"""
		if self.Ns <= 0:
			return matrix([])
		if not isscalar(time):
			raise Exception("time must be a scaler")

		return self.tocsr(time=time).todense()





	def dot(self,V,time=0):
		"""
		args:
			V, the vector to multiple with
			time=0, the time to evalute drive at.

		description:
			This function does the spare matrix vector multiplication of V with the Hamiltonian evaluated at 
			the specified time. It is faster in this case to multiple each individual parts of the Hamiltonian 
			first, then add all those vectors together.
		"""

		if self.Ns <= 0:
			return array([])
		if not isscalar(time):
			raise Exception("time must be a scaler")

		V=asarray(V)
		if self.Static_H != None: # if there is a static Hamiltonian...
			V_dot = self.Static_H.dot(V)	
			for ele in self.Dynamic_Hs:
				J=ele[0](time)
				V_dot += J*(ele[1].dot(V))
		else: # if there isn't...
			for ele in self.Dynamic_Hs:
				J=ele[0](time)
				V_dot += J*(ele[1].dot(V))

		return V_dot





	def MatrixElement(self,Vl,Vr,time=0):
		"""
		args:
			Vl, the vector to multiple with on left side
			Vr, the vector to multiple with on the right side
			time=0, the time to evalute drive at.

		description:
			This function takes the matrix element of the Hamiltonian at the specified time
			between Vl and Vr.
		"""
		if self.Ns <=0:
			return array([])

		Vl=asarray(Vl)
		Vr=asarray(Vr)
		HVr=self.dot(Vr,time=time)
		ME=dot(Vl.T.conj(),HVr)
		return ME





	def SparseEV(self,time=0,k=6,sigma=None,which='SA',maxiter=10000):
		"""
		args:
			time=0, the time to evalute drive at.
			other arguements see documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.eigsh.html
			
		description:
			function which diagonalizes hamiltonian using sparse methods
			solves for eigen values and eigen vectors, but can only solve for a few of them accurately.
			uses the scipy.sparse.linalg.eigsh function which is a wrapper for ARPACK
		"""
		if self.Ns <= 0:
			return array([]), array([[]])

		return eigsh(self.tocsr(time=time),k=k,sigma=sigma,which=which,maxiter=maxiter)





	def DenseEE(self,time=0):
		"""
		args:
			time=0, time to evaluate drive at.

		description:
			function which diagonalizes hamiltonian using dense methods solves for eigen values. 
			uses wrapped lapack functions which are contained in module py_lapack
		"""
		
		if self.Ns <= 0:
			return array([])

		return eigh(self.todense(time=time),JOBZ='N')





	def DenseEV(self,time=0):
		"""
		args:
			time=0, time to evaluate drive at.

		description:
			function which diagonalizes hamiltonian using dense methods solves for eigen values 
			and eigen vectors. uses wrapped lapack functions which are contained in module py_lapack
		"""
		if self.Ns <= 0:
			return array([]), array([[]])

		return eigh(self.todense(time=time))





	def evolve(self,v0,t0,time,real_time=True,verbose=False,**integrator_params):
		"""
		args:
			v0, intial wavefunction to evolve.
			t0, intial time 
			time, iterable or scalar, or time to evolve v0 to
			real_time, evolve real or imaginary time
			verbose, print times out as you evolve
			**integrator_params, the parameters used for the dop853 explicit rung-kutta solver.
			see documentation http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.integrate.ode.html

		description:
			This function uses complex_ode to evolve an input wavefunction.

		"""
		if self.Ns <= 0:
			return array([])

		v0=asarray(v0)

		if real_time:
			solver=complex_ode(lambda t,y:-1j*self.dot(y,time=t))
		else:
			solver=complex_ode(lambda t,y:-self.dot(y,time=t))

		solver.set_integrator("dop853", **integrator_params)
		solver.set_initial_value(v0,t=t0)
		
		if isscalar(time):
			if time==t0: return v0
			solver.integrate(time)
			if solver.successful():
				return solver.y
			else:
				raise Exception('failed to integrate')		
		else:
			sol=[]
			for t in time:
				if verbose: print t
				if t==t0: 
					sol.append(v0)
					continue
				solver.integrate(t)
				if solver.successful():
					sol.append(solver.y)
				else:
					raise Exception('failed to integrate')
			return sol




	# possiply impliment this in fortran using naive csr matrix vector dot product, might speed things up,
	# but maybe not since the exponential taylor converges pretty quickly. 
	def Exponential(self,V,z,time=0,n=1,atol=10**(-8)):
		"""
		args:
			V, vector to apply the matrix exponential on.
			a, the parameter in the exponential exp(aH)V
			time, time to evaluate drive at.
			n, the number of steps to break the expoential into exp(aH/n)^n V
			error, if the norm the vector of the taylor series is less than this number
			then the taylor series is truncated.

		description:
			this function computes exp(zH)V as a taylor series in zH. not useful for long time evolution.

		"""
		if self.Ns <= 0:
			return array([])
		if not isscaler(time):
			raise Exception("time must be a scaler")

		if n <= 0: raise Exception("n must be >= 0")

		V=asarray(V)
		for j in xrange(n):
			V1=array(V)
			e=1.0; i=1		
			while e > atol:
				V1=(z/(n*i))*self.dot(V1,time=time)
				V+=V1
				if i%2 == 0:
					e=norm(V1)
				i+=1
		return V


	# special methods to be added later
	"""
	def __add__(self,other):
		if isinstance(other,Hamiltonian1D):
			if self.Ns != other.Ns: raise Exception("cannot add Hamiltonians of different dimensions")
			new=deepcopy(other)
			new.Static_H+=self.Static+H
			new.Dynamic_Hs+=self.Dynamic+Hs
			return new
		else:
			raise Exception("Not Implimented")


	def __sub__(self,other):
		if isinstance(other,Hamiltonian1D):
			if self.Ns != other.Ns: raise Exception("cannot add Hamiltonians of different dimensions")
			new=deepcopy(other)
			new.StaticH-=self.Static_H
			for ele in self.Dynamic_Hs:
				new.DynamicHs.append((ele[0],-ele[1]))
			return new
		else:
			raise Exception("Not Implimented")
	"""

	



