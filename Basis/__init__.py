from Z_Basis import BasisError, Basis
from PBC_Basis import PeriodicBasis1D
from OBC_Basis import OpenBasis1D
import BitOps

class Basis1D:
	def __init__(self,Length,Nup=None,kblock=None,a=1,zblock=None,zAblock=None,zBblock=None,pblock=None,pzblock=None):

		# testing blocks for basis
		if (type(kblock) is int):
			if (type(zblock) is int) or (type(zAblock) is int) or (type(zBblock) is int) or (type(pblock) is int) or (type(pzblock) is int):
				#raise BasisError("Translation, spin inversion, and parity symmetries are not implimented at this time.")
				self.B=PeriodicBasis1D(Length,Nup=Nup,kblock=kblock,pblock=pblock,zblock=zblock,zAblock=zAblock,zBblock=zBblock,pzblock=pzblock,a=a)
			else:
				self.B=PeriodicBasis1D(Length,Nup=Nup,kblock=kblock,a=a)
		elif (type(zblock) is int) or (type(zAblock) is int) or (type(zBblock) is int) or (type(pblock) is int) or (type(pzblock) is int):
			self.B=OpenBasis1D(Length,Nup=Nup,zblock=zblock,zAblock=zAblock,zBblock=zBblock,pblock=pblock,pzblock=pzblock)
		else:
			self.B=Basis(Length,Nup=Nup)
		
		self.Nup=Nup
		self.kblock=kblock
		self.a=a
		self.zblock=zblock
		self.zAblock=zAblock
		self.zBblock=zBblock
		self.pblock=pblock
		self.pzblock=pzblock
		self.Ns=self.B.Ns

	def RefState(self,s):
		return self.B.RefState(s)

	def Op(self,J,opstr,indx):
		return self.B.Op(J,opstr,indx)





