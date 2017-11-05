#!/usr/bin/env python

"""
Author: David Ruffner, NYU
Date: 9/25/09

From template by: Jonathan Zrake, NYU CCPP


This program integrates a the Lane Emden using the following methods
forward Euler
Midpoint Method
RK4 Method
And it plots the results for a given initial condition.
"""
import pylab
import math



# Lane Emden Equation
#  w'' + 2*w'/z + w = 0
#  which can be rewritten in terms of v = w'
#  (d/dz)(w,v) = (v,-2*v/z-w**n)
#  Let function f(w,z,z) = (d/dz)(w,z)
#--------------------------------------------------------------------
n = 1
def f(x,z):
    '''this function receives input x=[w,v] and z, and returns the deriva
tives of w and v with respect to z in an array'''
    return pylab.array([x[1],-2*x[1]/z-((x[0])**n)])
hvalues = list()
L1_values_euler = list()
L1_values_mid = list()
L1_values_RK4 = list()
for i in range(5):
    

    a       = 0.0000001 # Lower and upper bounds of time the independent variable , z
    b       =  3.14
    n_samp  =  10**i + 1# Number of sample points to use in the integration
    z, dz   =  pylab.linspace(a, b, n_samp, retstep=True)
    x_euler = pylab.zeros([n_samp,2])
    x_mid = pylab.zeros([n_samp,2])
    x_RK4 = pylab.zeros([n_samp,2])
    hvalues.append(dz)

#Initial Conditions
    x_initial = pylab.array([1,0])


    x_euler = pylab.zeros([n_samp,2])
    x_mid = pylab.zeros([n_samp,2])
    x_RK4 = pylab.zeros([n_samp,2])

#Initial Conditions
    x_initial = pylab.array([1,0])

    x_euler = pylab.zeros([n_samp,2])
    x_mid = pylab.zeros([n_samp,2])
    x_RK4 = pylab.zeros([n_samp,2])

#Initial Conditions
    x_initial = pylab.array([1,0])


    x_euler = pylab.zeros([n_samp,2])
    x_mid = pylab.zeros([n_samp,2])
    x_RK4 = pylab.zeros([n_samp,2])

#Initial Conditions
    x_initial = pylab.array([1,0])


#ODE integrator

#Forward Euler integration
#-----------------------------------------------------------------------
    x_euler[0]= x_initial
    for i in range(n_samp-1):
        x_euler[i+1] = x_euler[i] + f(x_euler[i],z[i])*dz

#-----------------------------------------------------------------------

#Midpoint method integration
#-----------------------------------------------------------------------
    x_mid[0]= x_initial
    for i in range(n_samp-1):
        x1 = x_mid[i]+ dz*.5*f(x_mid[i],z[i])
        z1 = z[i]+dz*.5
        x_mid[i+1] = x_mid[i] + f(x1,z1)*dz
#-----------------------------------------------------------------------

#RK4 integration
#-----------------------------------------------------------------------
    x_RK4[0]=x_initial
    for i in range(n_samp-1):
        k1 = f(x_RK4[i],z[i])
        k2 = f(x_RK4[i] + dz*.5*k1, z[i] + dz*.5)
        k3 = f(x_RK4[i] + dz*.5*k2, z[i] + dz*.5)
        k4 = f(x_RK4[i] + dz*k3, z[i] + dz)
        x_RK4[i+1] = x_RK4[i] + (k1+2*k2+2*k3+k4)*dz/6
#----------------------------------------------------------------------

#Results
#-----------------------------------------------------------------------

    w_euler = [b[0] for b in x_euler]
    w_mid = [b[0] for b in x_mid]
    w_RK4 = [b[0] for b in x_RK4]

#w_exact = [1-(b**2)/6 for b in z]#n=0
    w_exact = [pylab.sin(b)/b for b in z]#n=1
#w_exact = [(1+(b**2)/3)**(-.5) for b in z]#n=5 watch out for using 1/2, that means 1 not .5!!

#Error estimates
#-----------------------------------------------------------------------
    absDw_euler = pylab.zeros(len(w_euler))
    for i in range(len(w_euler)):
        absDw_euler[i] = abs(w_euler[i]-w_exact[i])
    L1_euler = dz*sum(absDw_euler)

    absDw_mid = pylab.zeros(len(w_mid))
    for i in range(len(w_mid)):
        absDw_mid[i] = abs(w_mid[i]-w_exact[i])
    L1_mid = dz*sum(absDw_mid)

    absDw_RK4 = pylab.zeros(len(w_RK4))
    for i in range(len(w_RK4)):
        absDw_RK4[i] = abs(w_RK4[i]-w_exact[i])
    L1_RK4 = dz*sum(absDw_RK4)

    L1_values_euler.append(L1_euler)
    L1_values_mid.append(L1_mid)
    L1_values_RK4.append(L1_RK4)
    
print hvalues

#Plotting Results
#-------------------------------------------------------------------------------
pylab.subplot('211')
pylab.title("Integration of LaneEmden Eq. and analysis of convergence")
pylab.plot(z,w_euler,'-o', label=r"$Estimate of Solution, Eulers$")
pylab.plot(z,w_mid,'-o', label=r"$Estimate of Solution, Mid$")
pylab.plot(z,w_RK4,'-o', label=r"$Estimate of Solution, RK4$")
pylab.plot(z,w_exact,'-o', label=r"$exact solution$")
pylab.legend(loc='best')
pylab.xlabel(r"$z$")
pylab.ylabel(r"$w$")

pylab.subplot('212')
pylab.loglog(hvalues,L1_values_euler, '-o', label=r"$L1norm, Eulers$")
pylab.loglog(hvalues,L1_values_mid, '-o', label=r"$L1norm, Mid$")
pylab.loglog(hvalues,L1_values_RK4, '-o', label=r"$L1norm, RK4$")
pylab.legend(loc='best')
pylab.xlabel(r"$h$")
pylab.ylabel(r"$L1norm$")

pylab.show()
