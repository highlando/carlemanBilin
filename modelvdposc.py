import sympy as sp
import numpy as np
#the mean field flow model

#definition of the variables
a1, a2 = sp.symbols('a1, a2')
Vars = [a1, a2]
#definition of the parameters
mu = sp.symbols('mu')
Pars = [mu]

#definition of the F-function components
f1 = a1.subs(a1,1)
f2 = mu*(1-a1**2)*a2-a1

F = [f1, f2]

#definition of the k G-functions and its n components

#all coefficients must be symbols, not ints
g1, g2 = sp.symbols('g1, g2')
g1 = g1.subs(g1,1)
g2 = g2.subs(g2,0)

#G = [[g1, g2, g3],[g1, g2, g3]];
G = [[g1, g2]];

Var0 = np.array([[1],[0]])
