import numpy as np
import sympy as sp
import warnings

def taylorcoeffs(f,Vars,Order=None,point=None):
	"""Computes the Taylor coefficients of f

	f=f(Vars) up to order about the point, s.th.
	eg. f''[h,h] = f'' h kron h
	"""

	n = len(Vars)
	
	#default values for order and expansion point
	if Order is None:
		Order = 2
	if point is None:
		point = [0]*n

	pointDict = dict(zip(Vars,point))

	#f evaluated at point
	f0 = f.subs(pointDict)
	
	ft = [[f0]]

	fo = [f]
	fac = 1.
	for o in range(Order):
		#faculty factor in the expansion
		fac = fac*1/(o+1) 
		fco = []
		#Derivation 
		for fc in fo:
			for var in Vars:
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

def taylorcoeffs2matrix(taylC,Order=None,parPlusVals=None):
	"""list of taylC is evaluated to np.matrix

	1st list dim = state dim. 
	"""
	if Order is None: 
		Order = 2

	matT = []
	for o in range(Order+1):
		swpl = []
		for tidim in taylC:
			swpl.append(tidim[o])
		A = sp.Matrix(swpl);
		if parPlusVals is not None:
			A = A.subs(parPlusVals)
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

#TODO exception handling
	if not FJ[0].shape[0] == FJ[0].shape[1]:
		warnings.warn("Looks like f(0) is still in the list, I'll remove it for you",UserWarning)
		FJ = FJ[1:]
	if FJ[0].shape == (1,1):
		warnings.warn("If you encounter an error in this function try removing f(0) from the list",UserWarning)
	
	Order = len(FJ)
	sDim = FJ[0].shape[0]

#
#setup of the A matrix 	
	A = np.zeros((sDim,0))
#1st block row
	Az = np.zeros((sDim,0))
	for Fj in FJ:
		A = np.concatenate((A,Fj),axis=1)
	sfDim = A.shape[1]

#2nd+ rows
	for i in range(1,Order):
		Az = np.zeros((sDim**(i+1),Az.shape[1]+sDim**i))
		Ak = np.copy(Az)
		for Fj in FJ[:-i]:
			Ak = np.concatenate((Ak,fgcomponentij(Fj,i+1)),axis=1)
			
		A = np.concatenate((A,Ak),axis=0)
#
#setup of B and the list of Nk matrices 	
	B = np.zeros((sDim,0))
	NK = []
	for GkJ in GKJ:
		B = np.concatenate((B,GkJ[0]),axis=1)

#1st block row 
		Nk1 = GkJ[1]
		for Gkj in GkJ[2:-1]:
			Nk1 = np.concatenate((Nk1,Gkj),axis=1)

#2nd is still hard coded
		Nk2 = fgcomponentij(GkJ[0],2)
		for Gkj in GkJ[1:-2]:
			Nk2 = np.concatenate((Nk2,fgcomponentij(Gkj,2)),axis=1)

		Nk = np.concatenate((Nk1,Nk2),axis=0)

#3rd+ block rows
		if Order > 2:
			Nz = np.zeros((0,0))
			for i in range(3,Order+1):
				Nz = np.zeros((sDim**i,Nz.shape[1]+sDim**(i-2)))
				Nki = np.copy(Nz)

				for Gkj in GkJ[:-i]:
					Nki = np.concatenate((Nki,fgcomponentij(Gkj,i)),axis=1)

				Nk = np.concatenate((Nk,Nki),axis=0)

		Nk = np.concatenate((Nk,np.zeros((sfDim,sfDim-Nk.shape[1]))),axis=1)
		NK.append(Nk)

	B = np.concatenate((B,np.zeros((sfDim-sDim,len(GKJ)))),axis=0)

	return A, NK, B

def var0tobilinx0(Var0,Order):
	"""computes the initial value for the bilin approx
	"""
	
	X0bilin = np.copy(Var0)
	Aux0 = np.copy(Var0)
	for k in range(Order-1):
		Aux0 = np.kron(Aux0,Var0)
		X0bilin = np.concatenate((X0bilin,Aux0),0)

	return X0bilin



def evarhsbilinsystem(x,t,A,NK,B,u):

	Xnku = np.zeros((x.shape[0],1))
	for k, Nk in enumerate(NK):
		Xnku = Xnku + u[k]*np.dot(Nk,x)

	Xdot = np.dot(A,x)+Xnku+np.dot(B,u)

	return Xdot

def evarhsmodel(varVals,t,F,GK,parPlusVals,u,Vars):
	"""Evaluate the nonlinear xdot in the actual model
	"""
	varPlusVals = {}
	for i in range(len(F)):
		varPlusVals.update({Vars[i]:varVals[(i,0)]})

	Xdot = np.zeros((len(F),1))

	for i, fi in enumerate(F):
		fip = fi.subs(parPlusVals)
		Xdot[i] = sp.N(fip.subs(varPlusVals))
	
	for k, Gk in enumerate(GK):
		for i, gki in enumerate(Gk):
			gkip = gki.subs(parPlusVals)
			Xdot[i] = Xdot[i]+u[k]*sp.N(gkip.subs(varPlusVals))

	return Xdot

