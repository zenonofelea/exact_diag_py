from basis1d import basis
import numpy as _np


# gives the basis for the kronecker/Tensor product of two basis: b1 (x) b2
class tensor(basis):
	def __init__(self,b1,b2):
		if not isinstance(b1,basis):
			raise ValueError("b1 must be instance of basis class")
		if not isinstance(b2,basis):
			raise ValueError("b2 must be instance of basis class")

		self.b1=b1
		self.b2=b2

		self.Ns = b1.Ns*b2.Ns
		self.dtype = _np.min_scalar_type(-self.Ns)



	def Op(self,dtype,J,opstr,indx,pauli):
		i = opstr.index("|")
		indx1 = indx[:i]
		indx2 = indx[i:]

		opstr1,opstr2=opstr.split("|")

		ME1,row1,col1 = self.b1.Op(dtype,1.0,opstr1,indx1,pauli)
		ME2,row2,col2 = self.b2.Op(dtype,1.0,opstr2,indx2,pauli)

#		print ME1,row1,col1
#		print ME2,row2,col2

		n1 = row1.shape[0]
		n2 = row2.shape[0]

		row1 = _np.array(_np.broadcast_to(row1,(n2,n1)).T,dtype=self.dtype).ravel()
		row2 = _np.array(_np.broadcast_to(row2,(n1,n2)),dtype=self.dtype).ravel()
		row1 *= self.b2.Ns
		row = row1+row2

		del row1,row2

		col1 = _np.array(_np.broadcast_to(col1,(n2,n1)).T,dtype=self.dtype).ravel()
		col2 = _np.array(_np.broadcast_to(col2,(n1,n2)),dtype=self.dtype).ravel()
		col1 *= self.b2.Ns
		col = col1+col2

		del col1,col2

		ME1 = _np.array(_np.broadcast_to(ME1,(n2,n1)).T).ravel()
		ME2 = _np.array(_np.broadcast_to(ME2,(n1,n2))).ravel()
		ME = J*ME1*ME2

		del ME1,ME2

		return ME,row,col


		





