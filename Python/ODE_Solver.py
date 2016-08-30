#ODE SIR Model ODE solver

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

beta = .1
gamma = .9

def f(y,t):
    dy0 = -beta*y[0]*y[1]
    dy1 = beta*y[0]*y[1]-gamma*y[1]
    dy2 = gamma*y[1]
    
    return [dy0,dy1,dy2]
    
yinit = [80,10,10]
tspan = np.linspace(0,7,100)
y = odeint(f,yinit,tspan)

plt.plot(tspan,y[:,0],tspan,y[:,1],tspan,y[:,2])
plt.xlabel('time')
plt.ylabel('Population')
plt.title('SIR Model')


    