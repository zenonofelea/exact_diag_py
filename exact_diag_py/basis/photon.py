import numpy as _np
from basis1d import basis



class photons(basis):
	def __init__(self,Np):
		if (type(Np) is not int):
			raise ValueError("expecting integer for Np")

		self.Np = Np
		self.Ns = Np+1
		self.dtype = _np.min_scalar_type(-self.Ns)
		self.basis = _np.fromiter(xrange(self.Ns),dtype=self.dtype,count=self.Ns)

	def Op(self,dtype,J,opstr,*args):

		row = _np.array(self.basis)
		col = _np.array(self.basis)
		ME = _np.ones((self.Np+1,),dtype=dtype)
		for o in opstr[::-1]:
			if o == "I":
				continue
			elif o == "+":
				col += 1
				ME *= _np.sqrt(dtype(_np.abs(col)))
			elif o == "-":
				ME *= _np.sqrt(dtype(_np.abs(col)))
				col -= 1
			else:
				raise Exception("operator symbol {0} not recognized".format(o))

		mask = ( col >= 0)
		mask *= (col < (self.Ns))
		row = row[mask]
		col = col[mask]
		ME = ME[mask]
		ME *= J

		return ME,row,col		

			
