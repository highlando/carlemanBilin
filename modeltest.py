import sympy as sp
import numpy as np
#test model x'=alph*x*x+u

#definition of the variables
a1 = sp.Symbol('a1')
Vars = [a1]
#definition of the parameters
alph = sp.Symbol('alph')
Pars = [alph]

#definition of the F-function components
f1 = alph*a1*a1 

F = [f1]

#definition of the k G-functions and its n components

#all coefficients must be symbols, not ints
g1 = sp.Symbol('g1')
g1 = g1.subs(g1,1)

G = [[g1]];

Var0 = np.zeros((1,1))
