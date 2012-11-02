# naming conventions 
# i - state dimension
# j - order of expansion
# k - dimension of input
# Variables containing IJK are lists over all indeces

import sympy as sp
#import modeltest as model
import modelmff as model
reload(model)
import carlem as carlem
reload(carlem)


order = 2;

tfIJ = [] #taylor coefficients
for fi in model.F:
	tfIJ.append(carlem.taylorcoeffs(fi,model.varis))

tgKIJ = []
for gkI in model.G:
	tgkIJ = []
	for gki in gkI:
		tgkIJ.append(carlem.taylorcoeffs(gki,model.varis))
	tgKIJ.append(tgkIJ)

#compute the np matrices of the coefficients

parVals = dict(zip(model.params,[1]*len(model.params)))

FJ = carlem.taylorcoeffs2matrix(tfIJ,parVals=parVals)

GKJ = []
for gk in tgKIJ:
	GKJ.append(carlem.taylorcoeffs2matrix(gk,parVals=parVals))

