
## this is main.py but for nonzero f(0) and x(0)
# the bilin system is formulated for xd to give 
# the actual solution as x = xd + x0 , where x0 
# solves the nonlinear system with u=0

# naming conventions 
# i - state dimension
# j - order of expansion
# k - dimension of input
# Variables containing IJK are lists over all indeces


import sympy as sp
import numpy as np
import pylab as pl

from scipy import integrate

## importing the current model
# TODO use a class to define the models

import modeltest as model
#import modelmff as model
#import modelvdposc as model

## reload for use in ipython shell
reload(model)
import carlem as carlem
reload(carlem)

# set up the order of approximation: Order = 2, we will solve for [x,x*x]
Order = 2;

## Evaluation 
#  TODO handling and definition of control functions
#TODO parameter handling

parPlusVals = dict(zip(model.Pars,[1]*len(model.Pars)))

InpU = np.ones((1,1))

# the initial value of xd is 0! 
StateDim = len(model.Vars)
StateBilin = Order*StateDim if StateDim == 1 \
		else (1-StateDim**(Order+1))/(1-StateDim)
#cf. geometrical series
xd0 = np.zeros((StateBilin,1))

xdBil0 = np.concatenate((xd0,model.Var0),0)

# definition of the integration interval
t = np.arange(0,1.1,0.5)

# eva of the actual nonlinear model (1) and the shifted bilin 
xMod = integrate.odeint(nonhomocarlem, xdBil0, t, args = (model.F,model.G,parPlusVals,model.Vars,InpU))

## Plots 
# 
pl.figure(1)
pl.clf()
#pl.plot(xMod[:,0], xMod[:,1])
pl.plot(t,xMod[:,0])
pl.xlabel('$x$')
pl.ylabel('$\\dot{x}$')

#pl.figure(2)
#pl.clf()
#pl.plot(xBil[:,0], xBil[:,1])
#pl.plot(t,xBil[:,0])
#pl.xlabel('$x$')
pl.ylabel('$\\dot{x}$')

pl.show(block=False)


def nonhomocarlem(xModBil,t,F,G,Order,parPlusVals,Vars,InpU):
	StateDim = len(Vars)
	xMod = xModBil[:StateDim]
	xD = xModBil[StateDim:]
	## Taylor Coefficients 
	tfIJ = [] 
	for fi in F:
		tfIJ.append(carlem.taylorcoeffs(fi,Vars,point=xModBil,Order=Order))
	
	tgKIJ = []
	for gkI in G:
		tgkIJ = []
		for gki in gkI:
			tgkIJ.append(carlem.taylorcoeffs(gki,Vars,point=xModBil,Order=Order))
		tgKIJ.append(tgkIJ)
	
	## compute the np matrices of the coefficients
	#  thereto possible parameters are substituted
	
	FJ = carlem.taylorcoeffs2matrix(tfIJ,Order=Order,parPlusVals=parPlusVals)
	GKJ = []
	for gk in tgKIJ:
		GKJ.append(carlem.taylorcoeffs2matrix(gk,Order=Order,parPlusVals=parPlusVals))
	
	
	# setup of the coefficients of the approximation (2)
	A, NK, B = carlem.setupbilinsystem(FJ,GKJ)
	
	xdDot = carlem.evarhsbilinsystem( xD.flatten(), t, A,NK,B,InpU )

	# eva of the actual nonlinear model (1)
	xModDot = carlem.evarhsmodel( xMod.flatten(),t,F,G,parPlusVals,Vars,InpU)

	xDot = np.concatenate((xModDot,xdDot),0)

	return xDot.flatten()
