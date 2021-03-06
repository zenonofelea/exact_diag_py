from .base import basis

import numpy as _np
from scipy import sparse as _sp
from scipy.sparse import linalg as _sla
from scipy import linalg as _la


# gives the basis for the kronecker/Tensor product of two basis: b1 (x) b2
class tensor_basis(basis):

	def __init__(self,b1,b2):
		if not isinstance(b1,basis):
			raise ValueError("b1 must be instance of basis class")
		if not isinstance(b2,basis):
			raise ValueError("b2 must be instance of basis class")
		if isinstance(b1,tensor_basis): 
			raise TypeError("Can only create tensor basis with non-tensor type basis")
		if isinstance(b2,tensor_basis): 
			raise TypeError("Can only create tensor basis with non-tensor type basis")
		self._b1=b1
		self._b2=b2

		self._Ns = b1.Ns*b2.Ns
		self._dtype = _np.min_scalar_type(-self._Ns)

		self._operators = self._b1._operators +"\n"+ self._b2._operators



	def _get__str__(self):
		n_digits = int(_np.ceil(_np.log10(self._Ns)))
		str_list_1 = self._b1._get__str__()
		str_list_2 = self._b2._get__str__()
		Ns2 = self._b2.Ns
		temp = "\t{0:"+str(n_digits)+"d}  "
		str_list=[]
		for b1 in str_list_1:
			b1 = b1.split()
			s1 = b1[1]
			i1 = int(b1[0])
			for b2 in str_list_2:
				b2 = b2.split()
				s2 = b2[1]
				i2 = int(b2[0])
				str_list.append((temp.format(i2+Ns2*i1))+"\t"+s1+s2)

		if self._Ns > MAXPRINT:
			half = MAXPRINT//2
			str_list_1 = str_list[:half]
			str_list_2 = str_list[-half:]

			str_list = str_list_1
			str_list.extend(str_list_2)	

		return str_list		





	def get_vec(self,v0,sparse=True):
		if self._Ns <= 0:
			return _np.array([])

		if v0.ndim == 1:
			if v0.shape[0] != self._Ns:
				raise ValueError("v0 has incompatible dimensions with basis")
			v0 = v0.reshape((-1,1))
			return _combine_get_vec(self,v0,sparse)
		elif v0.ndim == 2:
			if v0.shape[0] != self._Ns:
				raise ValueError("v0 has incompatible dimensions with basis")
			return _combine_get_vecs(self,v0,sparse)
		else:
			raise ValueError("excpecting v0 to have ndim at most 2")




	def get_proj(self,dtype):
		proj1 = self._b1.get_proj(dtype)
		proj2 = self._b2.get_proj(dtype)

		return _sp.kron(proj1,proj2)


	def Op(self,opstr,indx,J,dtype,pauli):
		n=opstr.count("|")
		if n > 1: 
			raise ValueError("only one '|' charactor allowed")
		i = opstr.index("|")
		indx1 = indx[:i]
		indx2 = indx[i:]

		opstr1,opstr2=opstr.split("|")
		if not _np.can_cast(J,_np.dtype(dtype)):
			raise TypeError("can't cast J to proper dtype")

		if self._b1._Ns < self._b2._Ns:
			ME1,row1,col1 = self._b1.Op(opstr1,indx1,J,dtype,pauli)
			ME2,row2,col2 = self._b2.Op(opstr2,indx2,1.0,dtype,pauli)
		else:
			ME1,row1,col1 = self._b1.Op(opstr1,indx1,1.0,dtype,pauli)
			ME2,row2,col2 = self._b2.Op(opstr2,indx2,J,dtype,pauli)
			

		n1 = row1.shape[0]
		n2 = row2.shape[0]

		row1 = row1.astype(self._dtype)
		row1 *= self._b2.Ns
		row = _np.kron(row1,_np.ones_like(row2))
		row += _np.kron(_np.ones_like(row1),row2)

		del row1,row2

		col1 = col1.astype(self._dtype)
		col1 *= self._b2.Ns
		col = _np.kron(col1,_np.ones_like(col2))
		col += _np.kron(_np.ones_like(col1),col2)

		del col1,col2

		ME = _np.kron(ME1,ME2)

		del ME1,ME2

		return ME,row,col


		






def _combine_get_vec(basis,v0,sparse):
	Ns1=basis.b1.Ns
	Ns2=basis.b2.Ns

	Ns = min(Ns1,Ns2)

	# reshape vector to matrix to rewrite vector as an outer product.
	v0=_np.reshape(v0,(Ns1,Ns2))
	# take singular value decomposition to get which decomposes the matrix into separate parts.
	# the outer/tensor product of the cols of V1 and V2 are the product states which make up the original vector 
	if k<=0:
		V1,S,V2=_la.svd(v0)
	else:
		V1,S,V2=_sla.svds(v0,k=Ns-1,which='SM',maxiter=10**10)
		V12,[S2],V22=_sla.svds(v0,k=1,which='LM',maxiter=10**10)

		S.resize((Ns,))
		S[-1] = S2
		V1.resize((Ns1,Ns))
		V1[:,-1] = V12[:,0]
		V2.resize((Ns,Ns2))
		V2[-1,:] = V22[0,:]
		
		#V1 = _np.hstack(V1,V12)
		#V2 = _np.hstack(V2,V22)
	# svd returns V2.H so take the hc to reverse that
	V2=V2.T.conj()
	eps = _np.finfo(S.dtype).eps
	# for any values of s which are 0, remove those vectors because they do not contribute.
	mask=(S >= 10*eps)
	V1=V1[:,mask]
	V2=V2[:,mask]
	S=S[mask]


	# Next thing to do is take those vectors and convert them to their full hilbert space
	V1=basis.b1.get_vec(V1,sparse)
	V2=basis.b2.get_vec(V2,sparse)


	# calculate the dimension total hilbert space with no symmetries
	Ns=V2.shape[0]*V1.shape[0]		


	if sparse:
		v0=_sp.csr_matrix((Ns,1),dtype=V2.dtype)
		# combining all the vectors together with the tensor product as opposed to the outer product
		for i,s in enumerate(S):
			v1=V1.getcol(i)
			v2=V2.getcol(i)
			v=_sp.kron(v1,v2)
			v0 = v0 + s*v
		n=_np.sqrt(v0.multiply(v0.conj()).sum())
#		v0=v0/n
		v0=v0.astype(V1.dtype)
		
		
	else:
		v0=_np.zeros((Ns,),dtype=V2.dtype)
		for i,s in enumerate(S):
			v1=V1[:,i]
			v2=V2[:,i]
			v=_np.kron(v1,v2)
			v0 += s*v
		v0 /= _la.norm(v0)


	return v0




def _combine_get_vecs(basis,V0,sparse):
	v0_list=[]
	V0=V0.T
	for v0 in V0:
		v0_list.append(_combine_get_vec(basis,v0,sparse))

	if sparse:
		V0=_sp.hstack(v0_list)
	else:
		V0=_np.hstack(v0_list)

	return V0


