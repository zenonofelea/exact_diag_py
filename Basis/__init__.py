from Z_Basis import BasisError, Basis
from PBC_Basis import PeriodicBasis1D
from OBC_Basis import OpenBasis1D
import BitOps

class Basis1D:
	def __init__(self,Length,Nph,Nup=None,kblock=None,a=1,zblock=None,pblock=None,pzblock=None):

		# testing blocks for basis
		if (type(kblock) is int):
			if (type(zblock) is int) or (type(pblock) is int) or (type(pzblock) is int):
				self.B=PeriodicBasis1D(Length,Nup=Nup,kblock=kblock,pblock=pblock,zblock=zblock,pzblock=pzblock,a=a)
			else:
				self.B=PeriodicBasis1D(Length,Nup=Nup,kblock=kblock,a=a)
		elif (type(zblock) is int) or (type(pblock) is int) or (type(pzblock) is int):
			self.B=OpenBasis1D(Length,Nup=Nup,zblock=zblock,pblock=pblock,pzblock=pzblock)
		else:
			self.B=Basis(Length,Nph=Nph,Nup=Nup)
		
		self.Nph=Nph

		self.Nup=Nup
		self.kblock=kblock
		self.a=a
		self.zblock=zblock
		self.pblock=pblock
		self.pzblock=pzblock
		self.Ns=self.B.Ns
		self.Ns_tot=self.B.Ns_tot

	def RefState(self,s):
		return self.B.RefState(s)

	def Op(self,J,Sopstr,Popstr,indx,Nph):
		return self.B.Op(J,Sopstr,Popstr,indx)





