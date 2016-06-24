from exact_diag_py.basis import basis1d



# 1d spin symmetries
blocks=[]

for L in xrange(1,33,1):
	for a in xrange(1,L):
		if L%a == 0:
			for kblock in xrange(L/a):
				blocks.append((L,{'kblock':kblock,'a':a}))


	for Nup in xrange(L+1):
		for a in xrange(1,L):
			if L%a == 0:
				for kblock in xrange(L/a):
					blocks.append((L,{'kblock':kblock,'Nup':Nup,'a':a}))	




#(L,a,Nup,kblock)

with open("exact_diag_py/basis/basis1d/basis1d_Ns.py","w") as IO:
	line = "spin_basis_1d_Ns={ \n"
	IO.write(line)
	line = "\t {0}: {1}, \n"

	for L,block in blocks:
		try:
			b = basis1d(L,**block)
			kblock=b._blocks.get('kblock')
			Nup=b._blocks.get('Nup')
			a=b._blocks.get('a')
			tup = (L,a,Nup,kblock)
			if b.Ns != 0:
				print b.Ns,tup
				IO.write(line.format(tup,b.Ns))
		except ValueError:	continue

	IO.write("}\n")
	




