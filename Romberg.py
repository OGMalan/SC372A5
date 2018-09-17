import math
from math import e
from scipy.integrate import romberg as rb

def Trap(f,a,b,n):
    '''Returns the integral of a function for the range a and b estimated by means of the trapezoidal method using n steps.'''
    h = (b - a) / float(2**n)
    s = 0.5 * (f(a) + f(b))
    for i in range(1,n,1):
        s = s + f(a + i*h)
    print 'Integral of step '+str(n)+' is: '+ str(h*s)
    return h*s

def Richardson(f,a,b,tol):
    '''Returns the integral of a function for the range a,b '''
    i = 1
    Rj = Trap(f,a,b,i)
    running = True
    while running:
        j = i
        i += 1
        Rk = Trap(f,a,b,i)
        #Rj = (4 ** (j - 1) * Rk - Rj) / (4 ** (i - 1) - 1)
        Rj = (4 **(i - j) * Rk - Rj) / (4 ** (i - j) - 1)
        #E = abs((Rk - Rj)/((Rj / Rk)**2 - 1))
        E = abs((Rk - Rj)/ 2**(i))
        if E < tol:
            running = False
        Rj = Rk
        print 'Error is '+str(E)
    return Rk,j

def romberg(f, a, b, tol=1.0e-6):
    '''I, n_panels = romberg(f, a, b, tol=1.0e-6)
    Return the integral from a to b of f(x), by Romberg integration,
    as well as the number of panels used'''
    return Richardson(f,a,b,tol)

def f(x):
    return math.sin(x)


print 'Romberg value is: ' + str(romberg(f,1,10,tol=1e-6))
print 'Scipy value is: ' + str(rb(f,1,10,tol=1e-6))
