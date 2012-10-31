import sympy as sp
#the mean field flow model

#definition of the variables
a1, a2, a3 = sp.symbols('a1, a2, a3')
varis = [a1, a2, a3]
#definition of the parameters
sigs, oms, bet, gam = sp.symbols('sigs, oms, bet, gam')
params = [sigs, oms, bet, gam]

#definition of the F-function components
f1 = (sigs-bet*a3)*a1 + (oms + gam*a3)*a2
f2 = (sigs-bet*a3)*a2 - (oms + gam*a3)*a1
f3 = (sigs-bet*a3)*a3 + bet*(a1*a1 + a2*a2)

F = [f1, f2, f3]

#definition of the k G-functions and its n components

#all coefficients must be symbols, not ints
g1, g2, g3 = sp.symbols('g1, g2, g3')
g1 = g1.subs(g1,1)
g2 = g2.subs(g2,0)
g3 = g3.subs(g3,0)

G = [[g1, g2, g3],[g1, g2, g3]];


