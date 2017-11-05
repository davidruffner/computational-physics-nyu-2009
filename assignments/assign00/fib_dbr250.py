#!/usr/bin/env python

def fibonacci(N):

    terms = [0,1]
    
    for i in range(0,N-2):
       next = terms[-1]+terms[-2]
       terms.append(next)
    return terms

print fibonacci(10)
