import numpy as np
import sympy as sp

def taylorcoeffs(f,varis,order=None,point=None):
	"""Computes the Taylor coefficients of f

	f=f(varis) up to order about the point, s.th.
	eg. f''[h,h] = f'' h kron h
	"""

	n = len(varis)
	
	#default values for order and expansion point
	if order is None:
		order = 2
	if point is None:
		point = [0]*n

	pointDict = dict(zip(varis,point))

	#f evaluated at point
	f0 = f.subs(pointDict)
	
	ft = [[f0]]

	fo = [f]
	fac = 1.
	for o in range(order):
		#faculty factor in the expansion
		fac = fac*1/(o+1) 
		fco = []
		#Derivation 
		for fc in fo:
			for var in varis:
				fco.append(sp.diff(fc,var))

			fo = fco

		#Evaluation of current der. at point
		f0a = []
		for f0 in fco:
			f0 = f0.subs(pointDict)
			f0a.append(fac*f0)
		ft.append(f0a)
		#append a list of the current order derivations 

	return ft

def taylorcoeffs2matrix(taylC,order=None,parVals=None):
	"""list of taylC is evaluated to np.matrix

	1st list dim = state dim. 
	"""
	if order is None: 
		order = 2

	matT = []
	for o in range(order+1):
		swpl = []
		for tidim in taylC:
			swpl.append(tidim[o])
		A = sp.Matrix(swpl);
		if parVals is not None:
			A = A.subs(parVals)
		matT.append(np.array(np.array(A), np.float))

	return matT 

def fgcomponentij(Fj,i):
	"""compute components of the bilinearized sys matrix

	both for F and G 
	"""
	if i is 1:
		return Fj

	EyeM = np.eye(Fj.shape[0])
	SumFj = np.zeros([Fj.shape[0]**i,Fj.shape[1]*Fj.shape[0]**(i-1)])
	print([SumFj.shape])
	for summand in range(i):
		SFj = Fj
		for rFac in range(0,summand):
			SFj = np.kron(EyeM,SFj)
		for lFac in range(summand+1,i):
			SFj = np.kron(SFj,EyeM)

		SumFj = SumFj + SFj

	return SumFj

def setupbilinsystem(FJ,GKJ):
	"""setup the coefficients A,N,B from the taylor coeffs

	"""
	Order = len(FJ)
	sDim = FJ[0].shape[0]

#setup of the A matrix 	
	A = np.zeros((sDim,0))
	Az = np.zeros((sDim,0))
	for Fj in FJ:
		A = np.concatenate((A,Fj),axis=1)
	sfDim = A.shape[1]

	for i in range(1,Order):
		Az = np.zeros((sDim**(i+1),Az.shape[1]+sDim**i))
		Ak = np.copy(Az)
		for Fj in FJ[:i]:

			Ak = np.concatenate((Ak,fgcomponentij(Fj,i+1)),axis=1)
			
		A = np.concatenate((A,Ak),axis=0)

#setup of B and the list of Nk matrices 	
	B = np.zeros((sDim,0))
	for GkJ in GKJ:
		B = np.concatenate((B,GkJ[1]),axis=1)
		NK = []
		for Gkj in GkJ


	B = np.concatenate((B,np.zeros(sfDim-sDim,sDim))



	return A B

def evarhsbilinsystem(A,NK,B,u,x):

	return x
