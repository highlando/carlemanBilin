from pylab import *
from scipy import *

from scipy import integrate

def myfun(xvect, t):
    dxdt = zeros((2,))
    dxdt[0] = xvect[1]
    dxdt[1] = -10*sin(xvect[0]) - xvect[1]
    return dxdt

x0 = [-2,4]

t = arange(0,10,0.05)

x = integrate.odeint(myfun, x0, t)

figure(1)
clf()
plot(x[:,0], x[:,1])
xlabel('$x$')
ylabel('$\\dot{x}$')

show()