#!/usr/bin/env python

#The printed value given by this program should approach pi, since the ratio of # darts that hit the circle to the total number of darts, should be the same as the ratio of the area of the circle to the total area, which is PI. This is true# because there is an equal probability of a dart hitting any point in the area.Therefore in the limit of large numbers of darts, area is proportional to the number of darts.

def throw_darts(N):

    X = [ ]
    Y = [ ]

    import random

    for i in range(0,N):
        xValue = random.random()
        yValue = random.random()
        X.append(xValue)
        Y.append(yValue)
   
    return X, Y


def find_hits(X,Y):

    hits = 0

    for i in range(0,len(X)):
        if ((X[i]-.5)**2 + (Y[i]-.5)**2) <= .5**2:
            hits = hits+1

    return hits



N_darts = 10000
X_coords, Y_coords = throw_darts(N_darts)
print 4.0 * find_hits(X_coords, Y_coords) / N_darts
