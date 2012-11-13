# naming conventions 
# i - state dimension
# j - order of expansion
# k - dimension of input
# Variables containing IJK are lists over all indeces

import sympy as sp
import numpy as np
import pylab as pl

from scipy import integrate

#import modeltest as model
import modelmff as model
reload(model)
import carlem as carlem
reload(carlem)

Order = 2;

tfIJ = [] #taylor coefficients
for fi in model.F:
	tfIJ.append(carlem.taylorcoeffs(fi,model.Vars,Order=Order))

tgKIJ = []
for gkI in model.G:
	tgkIJ = []
	for gki in gkI:
		tgkIJ.append(carlem.taylorcoeffs(gki,model.Vars,Order=Order))
	tgKIJ.append(tgkIJ)

#compute the np matrices of the coefficients

#TODO parameter handling
parPlusVals = dict(zip(model.Pars,[1]*len(model.Pars)))

FJ = carlem.taylorcoeffs2matrix(tfIJ,Order=Order,parPlusVals=parPlusVals)
GKJ = []
for gk in tgKIJ:
	GKJ.append(carlem.taylorcoeffs2matrix(gk,Order=Order,parPlusVals=parPlusVals))

A, NK, B = carlem.setupbilinsystem(FJ,GKJ)

u = np.ones((1,1))
X0bilin = carlem.var0tobilinx0(model.Var0,Order)

#Xdotb = carlem.evarhsbilinsystem(X0bilin,t,A,NK,B,u,)
#Xdotm = carlem.evarhsmodel(model.Var0,t,model.F,model.G,parPlusVals,u,model.Vars)

t = np.arange(0,1,0.1)
x = integrate.odeint(carlem.evarhsbilinsystem, X0bilin, t, args = (A,NK,B,u))

t = arange(0,1,0.1)


