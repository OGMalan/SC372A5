import math
from math import e
from scipy.integrate import romberg as rb

def Trap(f,a,b,n):
    '''Returns the integral of a function for the range a and b estimated by means of the trapezoidal method using n steps.'''
    h = (b - a) / float(n)
    s = 0.5 * (f(a) + f(b))
    for i in range(1,n):
        s = s + f(a + i*h)
    print 'Integral of step '+str(n)+' is: '+ str(h*s)
    return h*s

def Richardson(f,a,b,tol):
    '''Returns the integral of a function for the range a,b '''
    n = 1
    Rj = Trap(f,a,b,n)
    running = True
    while running:
        j = n
        n += 1
        Rk = Trap(f,a,b,n)
        E = abs((Rk - Rj) / ((Rk / Rj)**2 - 1))
        if E < tol:
            running = False
        Rj = Rk
        print 'Error is '+str(E)
    return Rk,n

def romberg(f, a, b, tol=1.0e-6):
    '''I, n_panels = romberg(f, a, b, tol=1.0e-6)
    Return the integral from a to b of f(x), by Romberg integration,
    as well as the number of panels used'''
    return Richardson(f,a,b,tol)

def f(x):
    return math.sin(x)


print 'Romberg value is: ' + str(romberg(f,1,10,tol=1e-6))
print 'Scipy value is: ' + str(rb(f,1,10,tol=1e-6))
