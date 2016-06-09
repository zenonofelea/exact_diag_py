# python 2.7 modules
from array import array as vec
from numpy import sqrt
# local modules
from BitOps import * # loading modules for bit operations.
from SpinOps import SpinOp
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


def CheckStateZ(z,s,L,rz=2,sublat=None):
	#if sublat: rz=4;
	t=s
	t=flip_all(t,L,sublat=sublat)
	if t > s:
		rz*=1;
	else:
		rz=-1*abs(rz)
	return rz;



class OpenBasis1D(Basis):
	def __init__(self,L,Nup=None,pblock=None,zblock=None,zAblock=None,zBblock=None,pzblock=None):
		# This function in the constructor of the class:
		#		L: length of the chain
		#		Nup: number of up spins if restricting magnetization sector. 
		#		pblock: the number associated with parity quantum number of the block
		#		zblock: the number associated with spin inversion quantum number of the block
		#		pzblock: the number associated with parity + spin inversion quantum number of the block

		#	Note: the PZ block assumes the Hamiltonian is invariant under the total transformation PZ, 
		#				but not each transformation separately.

		Basis.__init__(self,L,Nup) # this calls the initialization of the basis class which initializes the basis list given Nup and Mcon/symm
		zbasis=self.basis # take initialized basis from Basis class and store in separate array to access, then overwrite basis.
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
			self.ZAcon = False
			self.ZBcon = False
			self.PZcon = True
			self.symm = True
			self.pblock = pblock
			self.zblock = zblock
			self.pzblock = pblock*zblock
			if (type(pzblock) is int) and (self.pzblock != self.pblock*self.zblock):
				print "OpenBasis1D wanring: contradiction between pzblock and pblock*zblock, assuming the block denoted by pblock and zblock" 
			self.Npz = []
			self.basis = []
			for s in zbasis:
				rpz = CheckStateZ(zblock,s,self.L)
				#if s==3: print "prelim rpz 1", s, rpz
				rpz = CheckStateP(pblock,s,self.L,rp=rpz)
				#if s==3: print "prelim rpz 2", s, rpz
				rpz = CheckStatePZ(pblock*zblock,s,self.L,rpz=rpz)
				#if s==3: print "prelim rpz 3", s, rpz
				#print rpz, int2bin(s,self.L)
				if rpz > 0:
					self.basis.append(s)
					self.Npz.append(rpz)
			self.Ns=len(self.basis)
		elif type(pblock) is int:
			if abs(pblock) != 1:
				raise BasisError("pblock must be either +/- 1")
			self.Pcon = True
			self.Zcon = False
			self.ZAcon = False
			self.ZBcon = False
			self.PZcon = False
			self.symm = True
			self.pblock = pblock
			self.z = zblock
			self.Np = []
			self.basis = []
			for s in zbasis:
				rp=CheckStateP(pblock,s,self.L)
#				print rp, int2bin(s,self.L)
				if rp > 0:
					self.basis.append(s)
					self.Np.append(rp)
			self.Ns=len(self.basis)
		elif type(zblock) is int:
			if abs(zblock) != 1:
				raise BasisError("zblock must be either +/- 1")
			self.Pcon = False
			self.Zcon = True
			self.ZAcon = False
			self.ZBcon = False
			self.PZcon = False
			self.symm = True
			self.z = zblock
			self.basis = []
			for s in zbasis:
				rz=CheckStateZ(zblock,s,self.L)
#				print rz, int2bin(s,self.L)
				if rz > 0:
					self.basis.append(s)
			self.Ns=len(self.basis)
		elif type(zAblock) is int and type(zBblock) is int:
			if abs(zAblock) != 1 or abs(zBblock) != 1:
				raise BasisError("zAblock and zBblock must be either +/- 1")
			self.Pcon = False
			self.Zcon = False
			self.ZAcon = True
			self.ZBcon = True
			self.PZcon = False
			self.symm = True
			self.zA = zAblock
			self.zB = zBblock
			self.basis = []
			for s in zbasis:
				rzA = CheckStateZ(zAblock,s,self.L,sublat='A')
				rzB = CheckStateZ(zBblock,s,self.L,sublat='B')
				rz  = CheckStateZ(zAblock*zBblock,s,self.L)
#				print rz, int2bin(s,self.L)
				if rzA > 0 and rzB > 0 and rz > 0:
					self.basis.append(s)
			self.Ns=len(self.basis)
		elif type(zAblock) is int:
			if abs(zAblock) != 1:
				raise BasisError("zAblock must be either +/- 1")
			self.Pcon = False
			self.Zcon = False
			self.ZAcon = True
			self.ZBcon = False
			self.PZcon = False
			self.symm = True
			self.zA = zAblock
			self.basis = []
			for s in zbasis:
				rz=CheckStateZ(zAblock,s,self.L,sublat='A')
#				print rz, int2bin(s,self.L)
				if rz > 0:
					self.basis.append(s)
			self.Ns=len(self.basis)
		elif type(zBblock) is int:
			if abs(zBblock) != 1:
				raise BasisError("zBblock must be either +/- 1")
			self.Pcon = False
			self.Zcon = False
			self.ZAcon = False
			self.ZBcon = True
			self.PZcon = False
			self.symm = True
			self.zB = zBblock
			self.basis = []
			for s in zbasis:
				rz=CheckStateZ(zBblock,s,self.L,sublat='B')
#				print rz, int2bin(s,self.L)
				if rz > 0:
					self.basis.append(s)
			self.Ns=len(self.basis)
		elif type(pzblock) is int:
			if abs(pzblock) != 1:
				raise BasisError("pzblock must be either +/- 1")
			self.PZcon = True
			self.Zcon = False
			self.Pcon = False
			self.ZAcon = False
			self.ZBcon = False
			self.symm = True
			self.pzblock = pzblock
			self.Npz = []
			self.basis = []
			for s in zbasis:
				rpz = CheckStatePZ(pzblock,s,self.L)
#				print rpz, int2bin(s,self.L)
				if rpz > 0:
					self.basis.append(s)
					self.Npz.append(rpz)
			self.Ns=len(self.basis)	
		else: 
			self.Pcon=False
			self.Zcon=False
			self.ZAcon=False
			self.ZBcon=False
			self.PZcon=False
		'''
		print "basis:"
		for i in self.basis:
			print i, int2bin(i,self.L)
		print "______"
		'''

	def RefState(self,s):
		# this function takes an integer s which represents a spin configuration in the Sz basis, then tries to find its 
		# reference state depending on the symmetries specified by the user. it does this by applying the various symmetry 
		# operations on the state and seeing whether a smaller integer is produced. This smaller integer by definition is the
		# reference state.
		# it returns r which is the reference state. g,q, and qg are the number of times the P,Z and PZ operators had to act.
		# This information is needed to calculate the matrix element s between states in this basis [1].
		t=s; r=s; g=0; gA=0; gB=0; q=0; qg=0;
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
		elif self.ZAcon and self.ZBcon:
			
			t = flip_all(t,self.L,sublat='A')
			if t < r:
				r=t; gA=1; gB=0;			
			t=s
			t = flip_all(t,self.L,sublat='B')
			if t < r:
				r=t; gA=0; gB=1;
			t=s
			t = flip_all(t,self.L)
			if t < r:
				r=t; gA=1; gB=1;
		elif self.ZAcon:
			t = flip_all(t,self.L,sublat='A')
			if t < s:
				r=t; gA=1;
		elif self.ZBcon:
			t = flip_all(t,self.L,sublat='B')
			if t < s:
				r=t; gB=1;
		elif self.PZcon:
			t = fliplr(t,self.L)
			t = flip_all(t,self.L)
			if t < s:
				r=t; qg=1;		

		return r,q,g,gA,gB,qg


	def Op(self,J,opstr,indx):
		# This function find the matrix elemement and state which opstr creates
		# after acting on an inputed state index.
		#		J: coupling in front of opstr
		#		st: index of a local state in the basis for which the opstor will act on
		#		opstr: string which contains a list of operators which  
		#		indx: a list of ordered indices which tell which operator in opstr live on the lattice.
		if self.Pcon or self.Zcon or self.ZAcon or self.ZBcon or self.PZcon: # if the user wants to use any symmetries, special care must be taken [1]
			ME_list=[]
			for st in xrange(self.Ns):
				s1=self.basis[st]
				ME,s2=SpinOp(s1,opstr,indx)
				s2,q,g,gA,gB,qg=self.RefState(s2)
				stt=self.FindZstate(s2)
				if stt >= 0:
					#print "(st, stt)=", (st, stt), "<--------" 
					if self.Pcon and self.Zcon:
						ME *= sqrt( float(self.Npz[stt])/self.Npz[st])*J*self.pblock**(q)*self.zblock**(g)
					elif self.Pcon:
						ME *= sqrt( float(self.Np[stt])/(self.Np[st]))*J*self.pblock**(q)
					elif self.Zcon:
						ME *= J*self.z**(g)
					elif self.ZAcon and self.ZBcon:
						#print "ME:", self.zA**(gA)*self.zB**(gB)
						ME *= J*self.zA**(gA)*self.zB**(gB)
					elif self.ZAcon:
						ME *= J*self.zA**(gA)
					elif self.ZBcon:
						ME *= J*self.zB**(gB)
					elif self.PZcon:
						ME *= sqrt( float(self.Npz[stt])/self.Npz[st] )*J*self.pzblock**(qg)		
				else:
					ME = 0.0
					stt = st
				ME_list.append([ME,st,stt])	
			return ME_list
		else: # else just use method from parent class.
			return Basis.Op(self,J,opstr,indx)






