from exact_diag_py.basis import basis1d



# 1d spin symmetries
blocks=[]

for L in xrange(1,33,1):
	for b1 in [-1,1]:
		for b2 in [-1,1]:
			blocks.append((L,{'zblock':b1,'pblock':b2}))
			blocks.append((L,{'zAblock':b1,'zBblock':b2}))

	for b1 in [-1,1]:
		blocks.append((L,{'pblock':b1}))
		blocks.append((L,{'pzblock':b1}))
		blocks.append((L,{'zblock':b1}))
		blocks.append((L,{'zAblock':b1}))
		blocks.append((L,{'zBblock':b1}))

	for Nup in xrange(L+1):
		for b1 in [-1,1]:
			for b2 in [-1,1]:
				blocks.append((L,{'zblock':b1,'pblock':b2,'Nup':Nup}))
				blocks.append((L,{'zAblock':b1,'zBblock':b2,'Nup':Nup}))

		for b1 in [-1,1]:
			blocks.append((L,{'pblock':b1,'Nup':Nup}))
			blocks.append((L,{'pzblock':b1,'Nup':Nup}))
			blocks.append((L,{'zblock':b1,'Nup':Nup}))
			blocks.append((L,{'zAblock':b1,'Nup':Nup}))
			blocks.append((L,{'zBblock':b1,'Nup':Nup}))



	for a in xrange(1,L):
		for kblock in xrange(L):
			blocks.append((L,{'kblock':kblock,'a':a}))
			for b1 in [-1,1]:
				for b2 in [-1,1]:
					blocks.append((L,{'zblock':b1,'pblock':b2,'kblock':kblock,'a':a}))
					blocks.append((L,{'zAblock':b1,'zBblock':b2,'kblock':kblock,'a':a}))

			for b1 in [-1,1]:
				blocks.append((L,{'pblock':b1,'kblock':kblock,'a':a}))
				blocks.append((L,{'zblock':b1,'kblock':kblock,'a':a}))
				blocks.append((L,{'zAblock':b1,'kblock':kblock,'a':a}))
				blocks.append((L,{'zBblock':b1,'kblock':kblock,'a':a}))	


	for Nup in xrange(L+1):
		for a in xrange(1,L):
			for kblock in xrange(L):
				blocks.append((L,{'kblock':kblock,'Nup':Nup,'a':a}))
				for b1 in [-1,1]:
					for b2 in [-1,1]:
						blocks.append((L,{'zblock':b1,'pblock':b2,'kblock':kblock,'Nup':Nup,'a':a}))
						blocks.append((L,{'zAblock':b1,'zBblock':b2,'kblock':kblock,'Nup':Nup,'a':a}))

				for b1 in [-1,1]:
					blocks.append((L,{'pblock':b1,'kblock':kblock,'Nup':Nup,'a':a}))
					blocks.append((L,{'zblock':b1,'kblock':kblock,'Nup':Nup,'a':a}))
					blocks.append((L,{'zAblock':b1,'kblock':kblock,'Nup':Nup,'a':a}))
					blocks.append((L,{'zBblock':b1,'kblock':kblock,'Nup':Nup,'a':a}))		




#(L,a,Nup,kblock,pblock,zblock,pzblock,zAblock,zBblock)

with open("exact_diag_py/basis/basis1d/basis1d_Ns.py","w") as IO:
	line = "spin_basis_1d_Ns={ \n"
	IO.write(line)
	line = "\t {0}: {1}, \n"

	for L,block in blocks:
		try:
			b = basis1d(L,**block)
			pblock=b._blocks.get('pblock')
			zblock=b._blocks.get('zblock')
			pzblock=b._blocks.get('pzblock')
			zAblock=b._blocks.get('zAblock')
			zBblock=b._blocks.get('zBblock')
			kblock=b._blocks.get('kblock')
			Nup=b._blocks.get('Nup')
			a=b._blocks.get('a')
			tup = (L,a,Nup,kblock,pblock,zblock,pzblock,zAblock,zBblock)
			if b.Ns != 0:
				print b.Ns,tup
				IO.write(line.format(tup,b.Ns))
		except ValueError:	continue

	IO.write("}\n")
	




