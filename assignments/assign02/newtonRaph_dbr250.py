#!/usr/bin/env python

'''
Author: David Ruffner, NYU
Date: 9/28/09

This program numerically finds the root of non-linear functions using the newton-Rapheson method

'''

import pylab

#the function to integrate
#------------------------------------------------------------------------
def func(x):
    return x - x**(1/3.) - 2

def funcPrime(x):
    return 1 - (1/3.)*x**(-2/3.)

def newtonStep(x):
    return x - func(x)/funcPrime(x)


#plot function to get a guess
#------------------------------------------------------------------------
#xvalues = pylab.linspace(.1,6,20)
#yvalues = [func(x) for x in xvalues]
#print xvalues
#print yvalues
#pylab.plot(xvalues,yvalues)
#pylab.show()


#Guess of root
#---------------------------------------------------------------------
print "Finding the root of f(x) = x - x^(1/3) - 2 = 0 by Newton's method"

x0 = 3

print "First guess: %f" % x0
print "iterations that converge on the root: "

tol = 10**(-20)#find root to this tolerance

#Iterate root finding method
#--------------------------------------------------------------
x = newtonStep(x0)


diff = abs(x0-x)
count = 1

while diff>=tol:
    print x
    x_old = x
    x = newtonStep(x)
    diff = abs(x-x_old)
    count = count + 1
    if count >= 100:
        break

