# python 2.7 modules
from array import array as vec
from numpy import sqrt
# local modules
from BitOps import * # loading modules for bit operations.
from SpinOps import SpinPhotonOp
from Z_Basis import Basis, BasisError

# References:
# [1]: A. W. Sandvik, AIP Conf. Proc. 1297, 135 (2010)


def CheckStatePZ(pz,s,L,rpz=2):
	t=s
	t=fliplr(t,L)
	t=flip_all(t,L)
	if t==s:
		if pz != -1:
			rpz*=2
		else:
			rpz=-1*abs(rpz)
	elif t > s:
		rpz*=1
	else:
		rpz=-1*abs(rpz)

	return rpz
		

def CheckStateP(p,s,L,rp=2):
	t=s
	t=fliplr(t,L)
	if t == s:
		if p != -1:
			rp*=2
		else:
			rp=-1*abs(rp)
	elif t > s: 
		rp*=1
	else:
		rp=-1*abs(rp)

	return rp;


def CheckStateZ(z,s,L,rz=2):
	t=s
	t=flip_all(t,L)
	if t > s:
		rz*=1;
	else:
		rz=-1*abs(rz)
	return rz;



class OpenBasis1D(Basis):
	def __init__(self,L,Nph=0,Nup=None,Ntot=None,pblock=None,zblock=None,pzblock=None):
		# This function in the constructor of the class:
		#		L: length of the chain
		#		Nup: number of up spins if restricting magnetization sector. 
		#		Nph: number of photons
		# 		Ntot: total number of spin ups plus the number of photons
		#		pblock: the number associated with parity quantum number of the block
		#		zblock: the number associated with spin inversion quantum number of the block
		#		pzblock: the number associated with parity + spin inversion quantum number of the block

		#	Note: the PZ block assumes the Hamiltonian is invariant under the total transformation PZ, 
		#				but not each transformation separately.

		Basis.__init__(self,L,Nph,Nup,Ntot) # this calls the initialization of the basis class which initializes the basis list given Nup and Mcon/symm
		
		sp_basis=self.sp_basis # take initialized basis from Basis class and store in separate array to access, then overwrite basis.
		basis=self.basis

		self.basis=vec('L')
		self.sp_basis=vec('L')

		self.Kcon=False # is false by definition
		self.kblock=None
		self.a=1

		# if symmetry is needed, the reference states must be found.
		# This is done through the CheckState function. Depending on
		# the symmetry, a different function must be used. Also if multiple
		# symmetries are used, the Checkstate functions be called
		# sequentially in order to check the state for all symmetries used.
		if (type(pblock) is int) and (type(zblock) is int):
			if abs(pblock) != 1:
				raise BasisError("pblock must be either +/- 1")
			if abs(zblock) != 1:
				raise BasisError("zblock must be either +/- 1")
			self.Pcon = True
			self.Zcon = True
			self.PZcon = True
			self.symm = True
			self.pblock = pblock
			self.zblock = zblock
			self.pzblock = pblock*zblock
			if (type(pzblock) is int) and (self.pzblock != self.pblock*self.zblock):
				print "OpenBasis1D wanring: contradiction between pzblock and pblock*zblock, assuming the block denoted by pblock and zblock" 
			self.Npz = []
			for s in sp_basis:
				rpz = CheckStateZ(zblock,s,self.L)
				rpz = CheckStateP(pblock,s,self.L,rp=rpz)
				rpz = CheckStatePZ(pblock*zblock,s,self.L,rpz=rpz)
				if rpz > 0:
					self.sp_basis.append(s)
					self.Npz.append(rpz)
			self.Ns=len(self.sp_basis)
		elif type(pblock) is int:
			if abs(pblock) != 1:
				raise BasisError("pblock must be either +/- 1")
			self.Pcon = True
			self.Zcon = False
			self.PZcon = False
			self.symm = True
			self.pblock = pblock
			self.z = zblock
			self.Np = []
			for s in sp_basis:
				rp=CheckStateP(pblock,s,self.L)
				if rp > 0:
					self.sp_basis.append(s)
					self.Np.append(rp)
			self.Ns=len(self.sp_basis)
		elif type(zblock) is int:
			if abs(zblock) != 1:
				raise BasisError("zblock must be either +/- 1")
			self.Pcon = False
			self.Zcon = True
			self.PZcon = False
			self.symm = True
			self.z = zblock
			for s in sp_basis:
				rz=CheckStateZ(zblock,s,self.L)
				if rz > 0:
					self.sp_basis.append(s)
			self.Ns=len(self.sp_basis)
		elif type(pzblock) is int:
			if abs(pzblock) != 1:
				raise BasisError("pzblock must be either +/- 1")
			self.PZcon = True
			self.Zcon = False
			self.Pcon = False
			self.symm = True
			self.pzblock = pzblock
			self.Npz = []
			for s in sp_basis:
				rpz = CheckStatePZ(pzblock,s,self.L)
				if rpz > 0:
					self.sp_basis.append(s)
					self.Npz.append(rpz)
			self.Ns=len(self.sp_basis)	
		else: 
			self.Pcon=False
			self.Zcon=False
			self.PZcon=False

		# add photon counterpart to the basis
		for sp in self.sp_basis:
			for ph in xrange(self.Nph+1):
				s = ElegantPair(sp,ph)
				if s in basis: #includes the symmetries from Z_Basis
					self.basis.append(s)
		self.Ns_tot=len(self.basis)	
	


	def RefState(self,s):
		# this function takes an integer s which represents a spin configuration in the Sz basis, then tries to find its 
		# reference state depending on the symmetries specified by the user. it does this by applying the various symmetry 
		# operations on the state and seeing whether a smaller integer is produced. This smaller integer by definition is the
		# reference state.
		# it returns r which is the reference state. g,q, and qg are the number of times the P,Z and PZ operators had to act.
		# This information is needed to calculate the matrix element s between states in this basis [1].
		t=s; r=s; g=0; q=0; qg=0;
		if self.Pcon and self.Zcon:
			t = flip_all(t,self.L)
			if t < r:
				r=t; g=1;q=0;
			t=s
			t = fliplr(t,self.L)
			if t < r:
				r=t; q=1; g=0;
			t=flip_all(t,self.L)
			if t < r:
				r=t; q=1; g=1;
		elif self.Pcon:
			t = fliplr(t,self.L)
			if t < s:
				r=t; q=1;
		elif self.Zcon:
			t = flip_all(t,self.L)
			if t < s:
				r=t; g=1;
		elif self.PZcon:
			t = fliplr(t,self.L)
			t = flip_all(t,self.L)
			if t < s:
				r=t; qg=1;		

		return r,q,g,qg


	def Op(self,J,Sopstr,Popstr,indx):
		# This function find the matrix elemement and state which opstr creates
		# after acting on an inputed state index.
		#		J: coupling in front of opstr
		#		st: index of a local state in the basis for which the opstor will act on
		#		opstr: string which contains a list of operators which  
		#		indx: a list of ordered indices which tell which operator in opstr live on the lattice.
		if self.Pcon or self.Zcon or self.PZcon: # if the user wants to use any symmetries, special care must be taken [1]
			ME_list=[]
			for st_tot in xrange(self.Ns_tot):

				s1=self.basis[st_tot]
				sp1, ph1 = ElegantUnpair(s1)
				# need the index in the spin basis, but some states can appear twice when symmetries are on, in which 
				# case ElegantUnpair returns the first found index. 
				if self.sp_basis.count(sp1)==2:
					if aux==True:
						st = self.sp_basis.index(sp1) + 1
					else:
						st = self.sp_basis.index(sp1)
						aux=True
				else:
					st = self.sp_basis.index(sp1)
					aux=False

				ME,sp2,ph2=SpinPhotonOp(s1,Sopstr,Popstr,indx,self.Nph)

				sp2,q,g,qg=self.RefState(sp2)
				stt=self.FindZstate(sp2)
				
				if stt >= 0:
					z = ElegantPair(sp2,ph2)
					if z in self.basis:
						stt_tot = self.basis.index(z) 
						if self.Pcon and self.Zcon:
							ME *= sqrt( float(self.Npz[stt])/self.Npz[st])*J*self.pblock**(q)*self.zblock**(g)
						elif self.Pcon:
							ME *= sqrt( float(self.Np[stt])/(self.Np[st]))*J*self.pblock**(q)
						elif self.Zcon:
							ME *=  J*self.z**(g)
						elif self.PZcon:
							ME *= sqrt( float(self.Npz[stt])/self.Npz[st] )*J*self.pzblock**(qg)
					else:
						ME = 0.0
						stt_tot = st_tot		
				else:
					ME = 0.0
					stt_tot = st_tot
				ME_list.append([ME,st_tot,stt_tot])	
			return ME_list
		else: # else just use method from parent class.
			return Basis.Op(self,J,Sopstr,Popstr,indx)






