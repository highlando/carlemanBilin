import sympy as sp
import modelMFF as model
reload(model)
import carlem
reload(carlem)

order = 2;

tF = [] #taylor coefficients
for fi in model.F:
	tF.append(carlem.taylorCoeffs(fi,model.varis))

tG = []
for gk in model.G:
	tGk = []
	for gki in gk:
		tGk.append(carlem.taylorCoeffs(gki,model.varis))
	tG.append(tGk)




	






