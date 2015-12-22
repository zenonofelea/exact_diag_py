from BitOps import * # loading modules for bit operations.
from SpinOps import SpinOp
from Z_Basis import Basis, BasisError
from array import array as vec
from numpy import pi,exp,sin,cos,sqrt,sign


# References:
# [1]: A. W. Sandvik, AIP Conf. Proc. 1297, 135 (2010)

'''
def CheckStateT(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	for i in xrange(1,L+1,T):
		t = shift(t,-T,L)
		if t < s:
			return -1
		elif t==s:
			if kblock % (L/i) != 0: return -1 # need to check the shift condition 
			return i
'''

def CheckStateT(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	R=-1
	for i in xrange(1,L+1,T):
		t = shift(t,-T,L)
		if t < s:
			return R
		elif t==s:
			if kblock % (L/i) != 0: return R # need to check the shift condition 
			R = i
			return R			



def CheckStateTP(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	R=-1
	m = -1 
	for i in xrange(1,L+1,T):
		t = shift(t,-T,L)
		if t < s:
			return R,m
		elif t==s:
			if kblock % (L/i) != 0: # need to check the shift condition 
				return R,m
			R = i
			#print R,s
			#return R,m
			t = s
			t = fliplr(t,L)
			for j in xrange(0,R):
				if t < s:
					R = -1
					return R,m
				elif t == s:
					m = j
					return R,m
				t = shift(t,-T,L) 
			return R,m




# child class of Basis, this is the momentum conserving basis:
# because it is a child class of Basis, it inherits its methods like FindZstate which searches basis for states
class PeriodicBasis1D(Basis):
	def __init__(self,L,Nup=None,kblock=None,pblock=None,zblock=None,pzblock=None,a=1):
		# This function in the constructor of the class:
		#		L: length of the chain
		#		Nup: number of up spins if restricting magnetization sector. 
		#		kblock: the number associated with the momentum block which basis is restricted to (i.e. k=2*pi*kblock/L)
 		#		a: number of lattice spaces between unit cells.

		Basis.__init__(self,L,Nup) # this calls the initialization of the basis class which initializes the basis list given Nup and Mcon/symm
		zbasis=self.basis # take initialized basis from Basis class and store in separate array to access, then overwrite basis.
		self.Pcon=False
		self.Zcon=False
		self.PZcon=False
		self.pblock=None
		self.zblock=None
		self.pzblock=None


		# if symmetry is needed, the reference states must be found.
		# This is done through the CheckState function. Depending on
		# the symmetry, a different function must be used. Also if multiple
		# symmetries are used, the Checkstate functions be called
		# sequentially in order to check the state for all symmetries used.

		if type(kblock) is int and type(pblock) is int:
			if kblock < 0 or kblock >= L: raise BasisError("0<= kblock < "+str(L))
			if abs(pblock) != 1: raise BasisError("pblock must have integer values +/-1")
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.pblock=pblock
			self.Kcon=True
			self.Pcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.R=vec('I') 
			self.m=[]#vec('I')
			self.basis=vec('L')

			if abs( sin(self.k) ) <= 1E-14:
				sigma_r = [+1]
				self.gk = 2
			else:
				sigma_r = [-1,+1]
				self.gk = 1
			for s in zbasis:
				r,m=CheckStateTP(kblock,L,s,T=a)
				#print [r,m]
				#some redundancy: if r=-1 I guess we need not do the following
				for sigma in sigma_r:
					if m != -1:
						if 1 + sigma*pblock*cos(kblock*m*2*pi/L) == 0:
							r = -1
						if (sigma == -1) and (1 - sigma*pblock*cos(kblock*m*2*pi/L) != 0):
							r = -1
					if r>0:
						#print [r,m], sigma*r
						self.R.append(sigma*r)
						self.m.append(m)	
						self.basis.append(s)
			self.Ns=len(self.basis)
		elif type(kblock) is int:
			if kblock < 0 or kblock >= L: raise BasisError("0<= kblock < "+str(L))
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.Kcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.R=vec('I') 
			self.basis=vec('L')
			for s in zbasis:
				r=CheckStateT(kblock,L,s,T=a)
				#print r
				if r > 0:
					self.R.append(r)
					self.basis.append(s)
			self.Ns=len(self.basis)
		else: 
			self.Kcon=False # do not change symm to False since there may be Magnetization conservation.
		
		print 'R =', self.R
		print 'list basis vectors:'
		for i in xrange(self.Ns):
			print [self.basis[i],int2bin( self.basis[i] ,L)]


	def RefState(self,s):
		# this function takes an integer s which represents a spin configuration in the Sz basis, then tries to find its 
		# reference state depending on the symmetries specified by the user. it does this by applying the various symmetry 
		# operations on the state and seeing whether a smaller integer is produced. This smaller integer by definition is the
		# reference state.
		# it returns r which is the reference state. l is the number of times the translation operator had to act.
		# This information is needed to calculate the matrix element s between states in this basis [1].
		t=s; r=s; l=0; q=0;
		if self.Kcon and self.Pcon:
			for i in xrange(1,self.L+1,self.a):
				t=shift(t,-self.a,self.L)
				if t < r:
					r=t; l=i;		
			t = s;
			t = fliplr(t,self.L)
			for i in xrange(1,self.L+1,self.a):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; q=1;		
		elif self.Kcon:
			for i in xrange(1,self.L+1,self.a):
				t=shift(t,-self.a,self.L)
				if t < r:
					r=t; l=i;
		return r,l,q


	def helement(self, st, stt, s1, s2, l):
		sigma = sign(self.R[st])
		if sign(self.R[st]) == sign(self.R[stt]): #sigma-diagonal elements
			if fliplr(s2,self.L) != shift(s2,-self.m[stt],self.L):
				return cos(self.k*l)
			else:
				return (cos(self.k*l) + sigma*self.pblock*cos(self.k*(l-self.m[stt]))  )/(1 + sigma*self.pblock*cos(self.k*self.m[stt]) )
		else:
			if fliplr(s2,self.L) != shift(s2,-self.m[stt],self.L):
				return -sigma*sin(self.k*l)
			else:
				return (-sigma*sin(self.k*l) + self.pblock*sin(self.k*(l-self.m[stt]))  )/(1 - sigma*self.pblock*cos(self.k*self.m[stt]) )

	
	def Nratios(self, st, stt): #Eq. (150) in [1], or sqrt( N_b^{-sigma}/ N_b^sigma )
		if abs( sin(self.k) ) >= 1E-14: #if semimomentum states exist
			sigma = sign(self.R[st])
			if fliplr(stt,self.L) != shift(stt,-self.m[stt],self.L):
				return 1.0
			else:
				return abs( sin(self.k*self.m[stt])) / ( 1.0 + sigma*self.pblock*cos(self.k*self.m[stt]) )
		else:
			return sqrt(float(self.R[st])/self.R[stt])

	def Nasigma(self,st): #Eq. (144) from [1]
		s1 = self.basis[st]
		sigma = sign(self.R[st])
		r = float( abs(self.R[st]) )
		gk = float(self.gk)
		t = fliplr(s1,self.L) #apply P
		if self.m[st] >= 0:
			return self.gk/r*(1 + sigma*self.pblock*cos(self.k*self.m[st]))
		else:
			return self.gk/r



	def Op(self,J,st,opstr,indx):	
		# This function finds the matrix elemement and states which opstr creates
		# after acting on an inputed state index.
		#		J: coupling in front of opstr
		#		st: index of a local state in the basis for which the opstor will act on
		#		opstr: string which contains a list of operators which  
		#		indx: a list of ordered indices which tell which operator in opstr live on the lattice.

		if self.Kcon or self.Pcon or self.Zcon or self.PZcon: # if the user wants to use symmetries, special care must be taken [1]
			s1=self.basis[st]
			ME,s2=SpinOp(s1,opstr,indx)
			s2,l,q=self.RefState(s2)
			#print [s1,s2]
			stt=self.FindZstate(s2) # if reference state not found in basis, this is not a valid matrix element.
			#print [st,stt]
			if stt >= 0:
				if self.Kcon and self.Pcon:
					"""
					#semimomentum states k \neq 0,\pi #
					#diagonal matrix elements: Eq. {13} in [1]
					if st>1 and (self.basis[st]==self.basis[st-1]):
						0
					elif st<self.Ns and (self.basis[st]==self.basis[st+1]):
						n=2
					else:
						n=1

					#offdiagonal matrix elements Eq. {16} in [1]
					if stt>1 and (self.basis[stt]==self.basis[stt-1]):
						m=2; stt = stt-1;
					elif stt<self.Ns and (self.basis[stt]==self.basis[stt+1]):
						m=2
					else:
						m=1
					for stt_j in xrange(stt, stt+m-1 + 1, 1):
						for st_i in xrange(st, st+n-1 + 1, 1):
							ME *= J*(sign(st)*self.pblock)**q*Nratios(st_i,stt_j) * helement(st_i,stt_j,l)
					"""
					sigma = sign(self.R[st])
					#print [s1,s2]
					#print [st,stt], [self.R[st],self.R[stt]]
					#print self.Nasigma(st), self.Nasigma(stt)
					#print self.helement(st,stt,l)
					ME *= J*(sigma*self.pblock)**q*sqrt( float( self.Nasigma(stt)/self.Nasigma(st)  ) ) * self.helement(st,stt,s1,s2,l)
					#ME *= J*(sigma*self.pblock)**q*self.Nratios(st, stt) * self.helement(st,stt,l)
				elif self.Kcon:
					#check sign of 1j in teh exponential
					#print [st,stt], [self.R[st],self.R[stt]]
					ME *= sqrt(float(self.R[st])/self.R[stt])*J*exp(-1j*self.k*l)
			else:
				ME=0.0;	stt=st
			return [ME,st,stt]
		else: # else, no special care is needed, just use the equivilant method from Basis class 
			return Basis.Op(self,J,st,opstr,indx)
		
		


