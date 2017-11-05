#!/usr/bin/python  
'''the -O turns off the assertions'''

def myFunc(a,b,c,d):
    ''' this is my function with
    the formula
    (2*a+b)/(c-d)
    '''
    assert c!=d, "c should not equal d"
    
    return (2*a+b)/(c-d)

myFunc(1,2,3,3)
