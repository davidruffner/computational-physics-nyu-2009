David Ruffner

Euler Equations Solver

This program uses the F_HLL finite volume method to solve the Euler equations 
for specific initial conditions

In this program there are the strong shock initial conditions.


INTRUCTIONS:
	
1) make			(This runs the Makefile which compiles the code)

2) python plot.py
			(This runs plot.py which calls eulers.c
			 which runs the numerical calculation. It also plots
			 the results)
3) If you want to run different shock conditions, these can be adjusted in 
   param.cfg


Other Necessary files
config.c	contains functions for getting parameters to and from files

euler.c		contains functions needed in eulerSolve.c

config.h	the header files for the above utility files
euler.h 

