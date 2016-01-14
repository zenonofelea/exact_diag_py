from BitOps import * # loading modules for bit operations.
from SpinOps import SpinOp
from Z_Basis import Basis, BasisError
from array import array as vec
from numpy import pi,exp,sin,cos,sqrt,sign


# References:
# [1]: A. W. Sandvik, AIP Conf. Proc. 1297, 135 (2010)


def CheckStateT(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	R=-1
	for i in xrange(1,L/T+1):
		t = shift(t,-T,L)
		if t < s:
			return R
		elif t==s:
			if kblock % (L/(T*i)) != 0: return R # need to check the shift condition 
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
	for i in xrange(1,L/T+1):
		t = shift(t,-T,L)
		if t < s:
			return R,m
		elif t==s:
			if kblock % (L/(T*i)) != 0: # need to check the shift condition 
				return R,m
			R = i
			break
	t = s
	t = fliplr(t,L)
	for j in xrange(0,R):
		if t < s:
			R = -1
			return R,m
		elif t == s:
			m = j
			break
		t = shift(t,-T,L) 
	return R,m


def CheckStateTZ(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	R=-1
	m = -1 
	for i in xrange(1,L/T+1):
		t = shift(t,-T,L)
		if t < s:
			return R,m
		elif t==s:
			if kblock % (L/(T*i)) != 0: # need to check the shift condition 
				return R,m
			R = i
			break
	t = s
	t = flip_all(t,L)
	for j in xrange(0,R):
		if t < s:
			R = -1
			return R,m
		elif t == s:
			m = j
			break
		t = shift(t,-T,L) 
	return R,m


def CheckStateT_PZ(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	R=-1
	m = -1 
	for i in xrange(1,L/T+1):
		t = shift(t,-T,L)
		if t < s:
			return R,m
		elif t==s:
			if kblock % (L/(T*i)) != 0: # need to check the shift condition 
				return R,m
			R = i
			break
	t = s
	t = flip_all(t,L)
	t = fliplr(t,L)
	for j in xrange(0,R):
		if t < s:
			R = -1
			return R,m
		elif t == s:
			m = j
			break
		t = shift(t,-T,L) 
	return R,m


def CheckStateTPZ(kblock,L,s,T=1):
	# this is a function defined in [1]
	# It is used to check if the integer inputed is a reference state for a state with momentum k.
	#		kblock: the number associated with the momentum (i.e. k=2*pi*kblock/L)
	#		L: length of the system
	#		s: integer which represents a spin config in Sz basis
	#		T: number of sites to translate by, not 1 if the unit cell on the lattice has 2 sites in it.
	t=s
	R=-1
	mpz = -1; mp = -1; mz = -1; 
	for i in xrange(1,L/T+1):
		t = shift(t,-T,L)
		if t < s:
			return R,mp,mz,mpz
		elif t==s:
			if kblock % (L/(T*i)) != 0: # need to check the shift condition 
				return R,mp,mz,mpz
			R = i
			break	
	t = s
	t = fliplr(t,L)
	for j in xrange(0,R):
		if t < s:
			R = -1
			return R,mp,mz,mpz
		elif t == s:
			mp = j
			break
		t = shift(t,-T,L) 

	t = s
	t = flip_all(t,L)
	for j in xrange(0,R):
		if t < s:
			R = -1
			return R,mp,mz,mpz
		elif t == s:
			mz = j
			break
		t = shift(t,-T,L)

	t = s
	t = flip_all(t,L)
	t = fliplr(t,L)
	for j in xrange(0,R):
		if t < s:
			R = -1
			return R,mp,mz,mpz
		elif t == s:
			mpz = j
			break
		t = shift(t,-T,L)	

	return R,mp,mz,mpz


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
		if (type(kblock) is int) and (L%a != 0):
			raise BasisError("lattice spacing a must be divisor of L.")

		# if symmetry is needed, the reference states must be found.
		# This is done through the CheckState function. Depending on
		# the symmetry, a different function must be used. Also if multiple
		# symmetries are used, the Checkstate functions be called
		# sequentially in order to check the state for all symmetries used.
		if type(kblock) is int and type(pblock) is int and type(zblock) is int:
			if abs(pblock) != 1: raise BasisError("pblock must have integer values +/-1")
			if abs(zblock) != 1: raise BasisError("zblock must have integer values +/-1")
			if self.L % 2 != 0: raise BasisError("chain length must be even!")
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.pblock=pblock
			self.zblock=zblock
			self.Kcon=True
			self.Pcon=True
			self.Zcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.Nasigma=[]#
			self.c=vec('I') 
			self.m=[]#vec('I')
			self.basis=vec('L')

			if abs( sin(self.k) ) <= 1E-14: #picks up k = 0, pi modes
				sigma_r = [+1]
			else:
				sigma_r = [-1,+1]
			for s in zbasis:
				for sigma in sigma_r:
					r,mp,mz,mpz=CheckStateTPZ(kblock,L,s,T=a)
					#print "before: ", sigma, [s,r], [mp,mz,mpz]
					if mp==-1 and mz==-1 and mpz==-1:
						if r>0:
							self.c.append(1)
							self.m.append(-1)
							Na = 2/float(r)
							self.Nasigma.append(sigma*Na)				
							self.basis.append(s)
					if mp != -1 and mz == -1 and mpz == -1:
						if 1 + sigma*pblock*cos(mp*self.k) == 0:
							r = -1
						if (sigma == -1) and (1 - sigma*pblock*cos(mp*self.k) != 0):
							r = -1
						if r>0:
							self.c.append(2)
							self.m.append(mp)
							Na = 2/float(r)*(1 + sigma*pblock*cos(mp*self.k))
							self.Nasigma.append(sigma*Na)				
							self.basis.append(s)
					if mp == -1 and mz != -1 and mpz == -1:
						if 1 + zblock*cos(self.k*mz) == 0:
							r=-1
						if r>0:		
							self.c.append(3)
							self.m.append(mz)
							Na = 2/float(r)*(1 + zblock*cos(mz*self.k))
							self.Nasigma.append(sigma*Na)				
							self.basis.append(s)
					if mp == -1 and mz == -1 and mpz != -1:
						if 1 + sigma*pblock*zblock*cos(mpz*self.k) == 0:
							r = -1
						if (sigma == -1) and (1 - sigma*pblock*zblock*cos(mpz*self.k) != 0):
							r = -1
						if r>0:
							self.c.append(4)
							self.m.append(mpz)
							Na = 2/float(r)*(1 + sigma*pblock*zblock*cos(mpz*self.k))
							self.Nasigma.append(sigma*Na)				
							self.basis.append(s)
					if (mp != -1) and (mz != -1):
						if (1 + sigma*pblock*cos(mp*self.k) == 0) or (1 + zblock*cos(mz*self.k) == 0):
							r = -1
						if (sigma == -1) and ( (1 - sigma*pblock*cos(mp*self.k) != 0) or (1 - zblock*cos(mz*self.k) != 0) ):
							r = -1
						if r>0:
							if (mpz - (mp+mz) % self.L) != 0:
								print "condition mpz = mp + mz failed: [mp,mz,mpz] =", [mp,mz,mpz]
							self.c.append(5)
							self.m.append(mp)
							Na = 2/float(r)*(1 + sigma*pblock*cos(mp*self.k))*(1 + zblock*cos(mz*self.k))
							self.Nasigma.append(sigma*Na)				
							self.basis.append(s)	
			self.Ns=len(self.basis)
		elif type(kblock) is int and type(pzblock) is int:
			if abs(pzblock) != 1: raise BasisError("pzblock must have integer values +/-1")
			if self.L % 2 != 0: raise BasisError("chain length must be even!")
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.pzblock=pzblock
			self.Kcon=True
			self.PZcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.Nasigma=[]#vec('I') 
			self.m=[]#vec('I')
			self.basis=vec('L')

			if abs( sin(self.k) ) <= 1E-14: #picks up k = 0, pi modes
				sigma_r = [+1]
			else:
				sigma_r = [-1,+1]
			for s in zbasis:
				for sigma in sigma_r:
					r,m=CheckStateT_PZ(kblock,L,s,T=a)
					#print "before: ", sigma, [s,r], [m]
					if m != -1:
						if 1 + sigma*pzblock*cos(m*self.k) == 0:
							r = -1
						if (sigma == -1) and (1 - sigma*pzblock*cos(m*self.k) != 0):
							r = -1
					#print "after: ", sigma, [s,r], [m]		
					if r>0:
						if m>=0:
							Na = 2/float(r)*(1 + sigma*self.pzblock*cos(self.k*m))
						else:
							Na = 2/float(r)
						self.Nasigma.append(sigma*Na)
						self.m.append(m)	
						self.basis.append(s)
			self.Ns=len(self.basis)	
		elif type(kblock) is int and type(pblock) is int:
			if abs(pblock) != 1: raise BasisError("pblock must have integer values +/-1")
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.pblock=pblock
			self.Kcon=True
			self.Pcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.Nasigma=[]#vec('I') 
			self.m=[]#vec('I')
			self.basis=vec('L')

			if abs( sin(self.k) ) <= 1E-14: #picks up k = 0, pi modes
				sigma_r = [+1]
			else:
				sigma_r = [-1,+1]
			for s in zbasis:
				for sigma in sigma_r:
					r,m=CheckStateTP(kblock,L,s,T=a)
					#print "before: ", sigma, [s,r], [m]
					if m != -1:
						if 1 + sigma*pblock*cos(m*self.k) == 0:
							r = -1
						if (sigma == -1) and (1 - sigma*pblock*cos(m*self.k) != 0):
							r = -1
					#print "after: ", sigma, [s,r], [m]		
					if r>0:
						if m>=0:
							Na = 1/float(r)*(1 + sigma*self.pblock*cos(self.k*m))
							#self.Nasigma.append(sigma*Na)
						else:
							Na = 1/float(r)
						self.Nasigma.append(sigma*Na)
						self.m.append(m)	
						self.basis.append(s)
			self.Ns=len(self.basis)
		elif type(kblock) is int and type(zblock) is int:
			if abs(zblock) != 1: raise BasisError("zblock must have integer values +/-1")
			if self.L % 2 != 0: raise BasisError("chain length must be even!")
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.zblock=zblock
			self.Kcon=True
			self.Zcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.NaTZ=[]#vec('I') 
			self.m=[]#vec('I')
			self.basis=vec('L')
			for s in zbasis:
				r,m=CheckStateTZ(kblock,L,s,T=a)
				#print "before:", [s,r,m]
				if m != -1:
					if 1 + zblock*cos(self.k*m) == 0:
						r=-1
				#print "after:", [s,r,m]		
				if r>0:
					if m >= 0:
						Na = 2/float(r)*(1 + self.zblock*cos(self.k*m))	
					else:
						Na = 2/float(r)
					self.NaTZ.append(Na)
					self.m.append(m)	
					self.basis.append(s)
			self.Ns=len(self.basis)	
		elif type(kblock) is int:
			self.a=a
			self.kblock=kblock
			self.k=2*pi*a*kblock/L
			self.Kcon=True
			self.symm=True # even if Mcon=False there is a symmetry therefore we must search through basis list.
			self.R=vec('I') 
			self.basis=vec('L')
			for s in zbasis:
				r=CheckStateT(kblock,L,s,T=a)
				#print [s,r]
				if r > 0:
					self.R.append(r)
					self.basis.append(s)
			self.Ns=len(self.basis)
		else: 
			self.Kcon=False # do not change symm to False since there may be Magnetization conservation.
			self.Pcon=False
			self.Zcon=False

		"""	
		if self.Kcon==True and self.Pcon==False and self.Zcon==False and self.PZcon==False:
			print 'R =', self.R
			print 'list basis vectors:'
			for i in xrange(self.Ns):
				print [self.basis[i],int2bin( self.basis[i] ,L)]
		elif self.Kcon and self.Zcon:
			print 'NaTZ =', self.NaTZ
			print 'm =', self.m
			print 'list basis vectors:'
			for i in xrange(self.Ns):
				print [self.basis[i],int2bin( self.basis[i] ,L)]
		else:
			print 'Nasigma =', self.Nasigma
			if self.Zcon or self.Pcon:
				print 'm =', self.m
			print 'list basis vectors:'
			for i in xrange(self.Ns):
				print [self.basis[i],int2bin( self.basis[i] ,L)]

		"""

	def RefState(self,s):
		# this function takes an integer s which represents a spin configuration in the Sz basis, then tries to find its 
		# reference state depending on the symmetries specified by the user. it does this by applying the various symmetry 
		# operations on the state and seeing whether a smaller integer is produced. This smaller integer by definition is the
		# reference state.
		# it returns r which is the reference state. l is the number of times the translation operator had to act.
		# This information is needed to calculate the matrix element s between states in this basis [1].
		t=s; r=s; l=0; q=0; g=0; qg=0;
		if self.Kcon and self.Pcon and self.Zcon:
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; g=0; q=0;	
			t = s;
			t = flip_all(t,self.L)
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; g=1; q=0;
			t = s;
			t = fliplr(t,self.L)
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; g=0; q=1;
			t = s;
			t = fliplr(t,self.L)
			t = flip_all(t,self.L)
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; g=1; q=1;		
		elif self.Kcon and self.Pcon:
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t < r:
					r=t; l=i;		
			t = s;
			t = fliplr(t,self.L)
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; q=1;	
		elif self.Kcon and self.Zcon:
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t < r:
					r=t; l=i;		
			t = s;
			t = flip_all(t,self.L)
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; g=1;
		elif self.Kcon and self.PZcon:
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t < r:
					r=t; l=i;		
			t = s;
			t = flip_all(t,self.L)
			t = fliplr(t,self.L)
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t<r:
					r=t; l=i; qg=1;								
		elif self.Kcon:
			for i in xrange(1,self.L/self.a+1):
				t=shift(t,-self.a,self.L)
				if t < r:
					r=t; l=i;

		return r,l,q,g,qg



	def helement(self, st, stt, l, q):
		sigma = sign(self.Nasigma[st])
		Nratios = sqrt( abs( float( self.Nasigma[stt]/self.Nasigma[st]  )   )  )

		if sign(self.Nasigma[st]) == sign(self.Nasigma[stt]): #sigma-diagonal elements
			if self.m[stt] < 0:
				cb = cos(self.k*l) #curly bracket
			else:
				cb = (cos(self.k*l) + sigma*self.pblock*cos(self.k*(l-self.m[stt]))  )/(1 + sigma*self.pblock*cos(self.k*self.m[stt]) )
		else: #sigma-offdiagonal elements
			if self.m[stt] < 0: #fliplr(s2,self.L) != shift(s2,-self.m[stt],self.L):
				cb = -sigma*sin(self.k*l)
			else:
				cb = (-sigma*sin(self.k*l) + self.pblock*sin(self.k*(l-self.m[stt]))  )/(1 - sigma*self.pblock*cos(self.k*self.m[stt]) )

		return (sigma*self.pblock)**q*Nratios*cb


	def helementT_PZ(self, st, stt, l, qg):
		sigma = sign(self.Nasigma[st])
		Nratios = sqrt( abs( float( self.Nasigma[stt]/self.Nasigma[st]  )   )  )

		if sign(self.Nasigma[st]) == sign(self.Nasigma[stt]): #sigma-diagonal elements
			if self.m[stt] < 0:
				cb = cos(self.k*l) #curly bracket
			else:
				cb = (cos(self.k*l) + sigma*self.pzblock*cos(self.k*(l-self.m[stt]))  )/(1 + sigma*self.pzblock*cos(self.k*self.m[stt]) )
		else: #sigma-offdiagonal elements
			if self.m[stt] < 0: #fliplr(s2,self.L) != shift(s2,-self.m[stt],self.L):
				cb = -sigma*sin(self.k*l)
			else:
				cb = (-sigma*sin(self.k*l) + self.pzblock*sin(self.k*(l-self.m[stt]))  )/(1 - sigma*self.pzblock*cos(self.k*self.m[stt]) )

		return (sigma*self.pzblock)**qg*Nratios*cb	
		

	def helementTPZ(self, st, stt, l, q, g):
		sigma = sign(self.Nasigma[st])
		Nratios = sqrt( abs( float( self.Nasigma[stt]/self.Nasigma[st]  )   )  )

		if sign(self.Nasigma[st]) == sign(self.Nasigma[stt]): #sigma-diagonal elements
			if self.c[stt]==1 or self.c[stt]==3:
				cb = cos(self.k*l) #curly bracket
			elif self.c[stt]==2 or self.c[stt]==5:
				cb = (cos(self.k*l) + sigma*self.pblock*cos(self.k*(l-self.m[stt]))  )/(1 + sigma*self.pblock*cos(self.k*self.m[stt]) )
			elif self.c[stt]==4:	
				cb = (cos(self.k*l) + sigma*self.pblock*self.zblock*cos(self.k*(l-self.m[stt]))  )/(1 + sigma*self.pblock*self.zblock*cos(self.k*self.m[stt]) )
		else: #sigma-offdiagonal elements
			if self.c[stt]==1 or self.c[stt]==3:
				cb = -sigma*sin(self.k*l)
			elif self.c[stt]==2 or self.c[stt]==5:
				cb = (-sigma*sin(self.k*l) + self.pblock*sin(self.k*(l-self.m[stt]))  )/(1 - sigma*self.pblock*cos(self.k*self.m[stt]) )
			elif self.c[stt]==4:
				cb = (-sigma*sin(self.k*l) + self.pblock*self.zblock*sin(self.k*(l-self.m[stt]))  )/(1 - sigma*self.pblock*self.zblock*cos(self.k*self.m[stt]) )
				
		return (sigma*self.pblock)**q*self.zblock**g*Nratios*cb

	def Op(self,J,opstr,indx):	
		# This function finds the matrix elemement and states which opstr creates
		# after acting on an inputed state index.
		#		J: coupling in front of opstr
		#		st: index of a local state in the basis for which the opstor will act on
		#		opstr: string which contains a list of operators which  
		#		indx: a list of ordered indices which tell which operator in opstr live on the lattice.

		if self.Kcon or self.Pcon or self.Zcon or self.PZcon: # if the user wants to use symmetries, special care must be taken [1]
			ME_list=[]
			for st in xrange(self.Ns):
				
				s1=self.basis[st]
				#offdiagonal matrix elements
				ME,s2=SpinOp(s1,opstr,indx); 
				#find reference state
				s2,l,q,g,qg=self.RefState(s2)
				stt=self.FindZstate(s2) # if reference state not found in basis, this is not a valid matrix element.

				#print [s1,s2], [st,stt]
				#print l,q,g,self.m[stt]
				if stt >= 0:
					if (self.Kcon and self.Pcon) or (self.Kcon and self.PZcon):

						if abs( sin(self.k) ) <= 1E-14: # k = 0, pi
							if self.Kcon and self.Pcon and self.Zcon:
								ME *= J*self.helementTPZ(st, stt, l, q, g)
							elif self.Kcon and self.PZcon:
								ME *= J*self.helementT_PZ(st, stt, l, qg)
							elif self.Kcon and self.Pcon:
								ME *= J*self.helement(st, stt, l, q)
							ME_list.append([ME,st,stt])						
						else:		
							# Eq. {13} in [1]
							if (st>0) and (self.basis[st]==self.basis[st-1]):
								continue # continue with next loop iteration
							elif (st<self.Ns-1) and (self.basis[st]==self.basis[st+1]):
								n=2
							else:									
								n=1

							me = 0;				
							if st == stt and l==0 and q==0 and g==0 and qg==0: #diagonal matrix elements:

								for st_i in xrange(st, st+n-1 +1, 1):
									if self.Kcon and self.Pcon and self.Zcon:
										me = (J*self.helementTPZ(st_i, st_i, l, q, g) )*ME
									elif self.Kcon and self.PZcon:
										me = (J*self.helementT_PZ(st_i, st_i, l, qg) )*ME
									elif self.Kcon and self.Pcon:
										me = (J*self.helement(st_i, st_i, l, q) )*ME 

									#print l,ME;
									ME_list.append([me,st_i,st_i])
							# Eq. {16} in [1]	
							else: #offdiagonal matrix elements 

								if (stt>0) and (self.basis[stt]==self.basis[stt-1]):
									m=2; stt = stt-1;
								elif (stt<self.Ns-1) and (self.basis[stt]==self.basis[stt+1]):
									m=2
								else:
									m=1
								for stt_j in xrange(stt, stt+m-1 +1, 1):
									for st_i in xrange(st, st+n-1+1 , 1):
										if self.Kcon and self.Pcon and self.Zcon:
											me = (J*self.helementTPZ(st_i, stt_j, l, q, g) )*ME
										elif self.Kcon and self.PZcon:
											me = (J*self.helementT_PZ(st_i, stt_j, l, qg) )*ME
										elif self.Kcon and self.Pcon:
											me = (J*self.helement(st_i, stt_j, l, q) )*ME
										ME_list.append([me,st_i,stt_j])				
					elif self.Kcon and self.Zcon:
						ME *= sqrt(float(self.NaTZ[stt]/self.NaTZ[st] ) )*J*self.zblock**g*exp(-1j*self.k*l)
						ME_list.append([ME,st,stt])		
					elif self.Kcon:
						ME *= sqrt(float(self.R[st])/self.R[stt])*J*exp(-1j*self.k*l)
						ME_list.append([ME,st,stt])
						#if st==stt: print l,ME;
				
				else:
					ME=0.0;	stt=st
					ME_list.append([ME,st,stt])

			return ME_list
		else: # else, no special care is needed, just use the equivilant method from Basis class 
			return Basis.Op(self,J,opstr,indx)