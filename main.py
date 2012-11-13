# This script sets up and evaluates a bilinear approximation
# of a nonlinear, input affine system of type
#
# 				x' = f(x) + g_k(x)*u_k , x(0) = 0 					(1)
#
# for nonlinear functions f,g_k , k = 1,...,K, a K-dimensional
# input u = [u_k] and f(0) = 0 . The approximation reads
#
# 				X' = A*X+u_k*N_k*X + B*u , X(0) = [x(0),x(0)*x(0),...]
# 																										(2)
#
# We use the Carleman bilinearization as described e.g. in 
# [Rugh1981] : 'Nonlinear System Theory' in Ch. 3.3

# The model (1) is provided in terms of symbolic variables and we
# use sympy to compute the coefficients of the Taylor expansion 
# symbolically

# Also (1) may still contain parameters

# For the evaluation the coefficients are transformed into
# numpy arrays

# the functions used are defined in carlem.py

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

## Taylor Coefficients 
tfIJ = [] 
for fi in model.F:
	tfIJ.append(carlem.taylorcoeffs(fi,model.Vars,Order=Order))

tgKIJ = []
for gkI in model.G:
	tgkIJ = []
	for gki in gkI:
		tgkIJ.append(carlem.taylorcoeffs(gki,model.Vars,Order=Order))
	tgKIJ.append(tgkIJ)

## compute the np matrices of the coefficients
#  thereto possible parameters are substituted
#TODO parameter handling
parPlusVals = dict(zip(model.Pars,[1]*len(model.Pars)))

FJ = carlem.taylorcoeffs2matrix(tfIJ,Order=Order,parPlusVals=parPlusVals)
GKJ = []
for gk in tgKIJ:
	GKJ.append(carlem.taylorcoeffs2matrix(gk,Order=Order,parPlusVals=parPlusVals))

# the initial value is to be extended as well
X0bilin = carlem.var0tobilinx0(model.Var0,Order)

# setup of the coefficients of the approximation (2)
A, NK, B = carlem.setupbilinsystem(FJ,GKJ)

## Evaluation 
#  TODO handling and definition of control functions

InpU = np.ones((1,1))

# definition of the integration interval
t = np.arange(0,1.1,0.5)

# eva of the bilin appr. (2)
xBil = integrate.odeint(carlem.evarhsbilinsystem, X0bilin.flatten(), t, args = (A,NK,B,InpU))

# eva of the actual nonlinear model (1)
xMod = integrate.odeint(carlem.evarhsmodel, model.Var0.flatten(), t, args = (model.F,model.G,parPlusVals,model.Vars,InpU))

## Plots 
# 
pl.figure(1)
pl.clf()
#pl.plot(xMod[:,0], xMod[:,1])
pl.plot(t,xMod[:,0])
pl.xlabel('$x$')
pl.ylabel('$\\dot{x}$')

pl.figure(2)
pl.clf()
#pl.plot(xBil[:,0], xBil[:,1])
pl.plot(t,xBil[:,0])
pl.xlabel('$x$')
pl.ylabel('$\\dot{x}$')

pl.show(block=False)
