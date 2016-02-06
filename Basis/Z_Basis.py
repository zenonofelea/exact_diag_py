# python 2.7 modules
import operator as op # needed to calculate n choose r in function ncr(n,r).
# array is a very nice data structure which stores values in a c or fortran like array, saving memory, but has all features of a list.
# it is not good for array type operations like multiplication, for those use numpy arrays.
from array import array as vec

# local modules
from SpinOps import SpinPhotonOp # needed to act with opstr
from BitOps import * # loading modules for bit operations.

import numpy as np

# References:
# [1]: A. W. Sandvik, AIP Conf. Proc. 1297, 135 (2010)





class BasisError(Exception):
	# this class defines an exception which can be raised whenever there is some sort of error which we can
	# see will obviously break the code. 
	def __init__(self,message):
		self.message=message
	def __str__(self):
		return self.message



def ncr(n, r):
# this function calculates n choose r used to find the total number of basis states when the magnetization is conserved.
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

# Parent Class Basis: This class is the basic template for all other basis classes. It only has Magnetization symmetry as an option.
# all 'child' classes will inherit its functionality, but that functionality can be overwritten in the child class.
# Basis classes must have the functionality of finding the matrix elements built in. This way, the function for constructing 
# the hamiltonian is universal and the basis object takes care of constructing the correct matrix elements based on its internal symmetry. 
class Basis:
	def __init__(self,L,Nph=0,Nup=None):
		# This is the constructor of the class Basis:
		#		L: length of chain
		#		Nup: number of up spins if restricting magnetization sector. 
		self.L=L
		self.Nph=Nph
		if (type(Nup) is int) and (type(Nph) is int):
			if Nup < 0 or Nup > L: raise BasisError("0 <= Nup <= "+str(L))
			if Nph < 0: raise BasisError("0 <= Nph ")
			self.Nup=Nup
			self.Mcon=True 
			self.symm=True # Symmetry exsists so one must use the search functionality when calculating matrix elements
			self.Ns=ncr(L,Nup)
			self.Ns_tot=(Nph+1)*ncr(L,Nup) 
			zbasis=vec('L')
			sp_zbasis = vec('L')
			sp=sum([2**i for i in xrange(0,Nup)])
			sp_zbasis.append(sp)
			for ph in xrange(self.Nph+1):
				zbasis.append( ElegantPair(sp,ph) )
			for i in xrange(self.Ns-1):
				t = (sp | (sp - 1)) + 1
				sp = t | ((((t & -t) / (sp & -sp)) >> 1) - 1)
				sp_zbasis.append(sp)
				for ph in xrange(self.Nph+1):
					s = ElegantPair(sp,ph)
					print [sp,ph], s
					zbasis.append(s)
		else:
			self.Ns=2**L
			self.Ns_tot=(Nph+1)*2**L
			self.Mcon=False
			self.symm=False # No symmetries here. at all so each integer corresponds to the number in the hilbert space.
			sp_zbasis=xrange(self.Ns)
			zbasis=[]
			for sp in xrange(self.Ns):
				for ph in xrange(self.Nph+1):
					print [sp,ph]
					s = ElegantPair(sp,ph)
					zbasis.append(s)
			#zbasis=xrange(self.Ns)

		self.basis=zbasis # total spin + photon basis
		self.sp_basis=sp_zbasis #spin basis only
		print zbasis
		print sp_zbasis


	def FindZstate(self,s):
		if self.symm:
			bmin=0;bmax=self.Ns-1
			while True:
				b=(bmin+bmax)/2
				if s < self.sp_basis[b]:
					bmax=b-1
				elif s > self.sp_basis[b]:
					bmin=b+1
				else:
					return b
				if bmin > bmax:
					return -1
		else: return s


	def Op(self,J,Sopstr,Popstr,indx):
		ME_list=[]
		for st in xrange(self.Ns_tot):
			s1=self.basis[st]
			ME,sp2,ph2=SpinPhotonOp(s1,Sopstr,Popstr,indx,self.Nph)
			#print st, s1
			#print [sp2, ph2], ElegantPair( sp2, ph2 )
			
			stt = self.FindZstate(sp2)

			if stt>=0:
				stt = self.basis.index( ElegantPair( sp2, ph2 ) )
			ME_list.append([J*ME,st,stt])

		return ME_list



