import sympy as sp
import numpy as np
import pylab as pl

from scipy import integrate

def dydt(x,t):
	A = np.ones((2,2))

	xdot = np.dot(A,x)#+np.dot(B,t*u)

	return xdot

t = np.arange(0,1,0.1)

A = np.ones((2,2))
B = np.ones((2,1))
x0 = np.ones((2,1))
u = np.ones((1,1))

x = integrate.odeint(dydt, x0, t)
