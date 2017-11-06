#!/usr/bin/env python

"""
Author: David Ruffner, NYU
Date: 9/21/09

From template by: Jonathan Zrake, NYU CCPP


This program integrates a simple ODE using the following methods
forward Euler
Midpoint Method
RK4 Method
And it plots the results for a given initial condition.
"""
import pylab


k=1
nu=1
a=1
w=1

# ODE of form
#  x'(t) = -k*x(t)/nu + A*cos(w*t)/nu
#  let x'(t)=f(x,t)

# The function equal to the derivative
#--------------------------------------------------------------------
def f(x,t):
    return -x + pylab.cos(t)

a       = 0 # Lower and upper bounds of time the independent variable , t
b       =  80
n_samp  =  80  # Number of sample points to use in the integration

t, dt   =  pylab.linspace(a, b, n_samp, retstep=True)

x_euler = pylab.zeros(n_samp)
x_mid = pylab.zeros(n_samp)
x_RK4 = pylab.zeros(n_samp)

#Initial Conditions
x_initial= 3

#ODE integrator

#Forward Euler integration
#-----------------------------------------------------------------------
x_euler[0]= x_initial
for i in range(n_samp-1):
    x_euler[i+1] = x_euler[i] + f(x_euler[i],t[i])*dt

#-----------------------------------------------------------------------

#Midpoint method integration
#-----------------------------------------------------------------------
x_mid[0]= x_initial
for i in range(n_samp-1):
    x1 = x_mid[i]+ dt*.5*f(x_mid[i],t[i])
    t1 = t[i]+dt*.5
    x_mid[i+1] = x_mid[i] + f(x1,t1)*dt
#-----------------------------------------------------------------------

#RK4 integration
#-----------------------------------------------------------------------
x_RK4[0]=x_initial
for i in range(n_samp-1):
    k1 = f(x_RK4[i],t[i])
    k2 = f(x_RK4[i] + dt*.5*k1, t[i] + dt*.5)
    k3 = f(x_RK4[i] + dt*.5*k2, t[i] + dt*.5)
    k4 = f(x_RK4[i] + dt*k3, t[i] + dt)
    x_RK4[i+1] = x_RK4[i] + (k1+2*k2+2*k3+k4)*dt/6
#----------------------------------------------------------------------

#Plotting the Results
#-----------------------------------------------------------------------
pylab.plot(t,x_euler,'-o', label=r"$Estimate of Solution, Eulers$")
pylab.plot(t,x_mid,'-o', label=r"$Estimate of Solution, Midpoint Method$")
pylab.plot(t,x_RK4,'-o', label=r"$Estimat of Solution, RK4$")
pylab.legend(loc='best')
pylab.xlabel(r"$t$")
pylab.ylabel(r"$x$")
pylab.show()
