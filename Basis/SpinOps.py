from BitOps import *
import numpy as np


class SpinPhotonOpError(Exception):
	def __init__(self,message):
		self.message=message
	def __str__(self):
		return self.message



def SpinPhotonOp(s,Sopstr,Popstr,indx,Nph=0):
	if len(indx) != len(Sopstr):
		raise SpinPhotonOpError("Dimension mismatch of Sopstr and indx")

	#print "TYPE:", type(s)
	# decompose integer s into the spin sp and photon ph integers
	sp, ph = ElegantUnpair(s)

	#print "VALUES:", sp, ph
	#print "TYPE:", type(sp), type(ph)

	ME=1.0; r=sp; n = ph;
	NSops=len(Sopstr)
	NPops=len(Popstr)

	for i in xrange(NSops-1,-1,-1): #string is written left to right, but operators act in reverse order.
		if Sopstr[i] == "I": #constant
			continue
		elif Sopstr[i] == "z":
			a=testBit(r,indx[i])
			ME*=(-1)**(a+1)/2.0
		elif Sopstr[i] == "x":
			r=flipBit(r,indx[i])
			ME*=1/2.0
		elif Sopstr[i] == "y":
			a=testBit(r,indx[i])
			r=flipBit(r,indx[i])
			ME*=(-1)**(a)/2.0j
		elif Sopstr[i] == "+":
			a=testBit(r,indx[i])
			if a == 1: r=sp; ME=0.0; break # + operator kills the state --> end loop
			else: r=flipBit(r,indx[i])
		elif Sopstr[i] == "-":
			a=testBit(r,indx[i])
			if a == 0: r=sp; ME=0.0; break # - operator kills the state --> end loop
			else: r=flipBit(r,indx[i])
		else:
			raise SpinPhotonOpError("operator symbol "+Sopstr[i]+" not recognized")


	for i in xrange(NPops-1,-1,-1): # string is written left to right, but operators act in reverse order.
		if Popstr[i] == "I": # constant
			continue
		elif Popstr[i] == "d": # photon destruction operator
			if n == 0: n=ph; ME=0.0; break # d operator kills the state --> end loop
			else: ME*=np.sqrt(n); n = n-1; 
		elif Popstr[i] == "c": # photon creation operator
			if n == Nph: n=ph; ME=0.0; break # c operator kills the state --> end loop // finite-size constraint on HO basis
			else: ME*=np.sqrt(n+1); n = n+1;
		else:
			raise SpinPhotonOpError("operator symbol "+Popstr[i]+" not recognized")


#	print Sopstr, indx, ME
	if ME.imag == 0.0:
		ME=ME.real

	return ME, r, n

