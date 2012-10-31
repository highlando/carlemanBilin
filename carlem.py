import numpy as np
import sympy as sp

def taylorCoeffs(f,varis,order=None,point=None):
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

def taylorCoeffs2Matrix(taylC,order=None,parsVals=None):
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
		if parsVals is not None:
			A = A.subs(parsVals)
		matT.append(np.array(np.array(A), np.float))

	return matT 

