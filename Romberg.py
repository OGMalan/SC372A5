import math, time
from math import e
from scipy.integrate import romberg as rb

def Trap(f,a,b,n):
    '''Returns the integral of a function for the range a and b estimated by means of the trapezoidal method using n steps.'''
    h = (b - a) / float(n)
    s = 0.5 * (f(a) + f(b))
    for i in range(1,n):
        s = s + f(a + i*h)
    #print 'Integral of step '+str(n)+' is: '+ str(h*s)
    return h*s

def Richardson(f,a,b,tol):
    '''Iteratively applies Richardson extrapolation to increase the accuracy of the estimate of the integral.'''
    j = 1
    k = 1
    Rj = Trap(f,a,b,j)
    running = True
    while running:
        j = k
        k += 1
        Rk = Trap(f,a,b,k)
        Rj = (4 **(k - j) * Rk - Rj) / (4 ** (k - j) - 1)
        E = abs((Rk - Rj))
        if E < tol:
            running = False
        Rj = Rk
        #print 'Error is '+str(E)
    return Rk,j

def romberg(f, a, b, tol):
    '''I, n_panels = romberg(f, a, b, tol=1.0e-6)
    Return the integral from a to b of f(x), by Romberg integration,
    as well as the number of panels used'''
    return Richardson(f,a,b,tol)

def f(x):
    return e**x

start = time.clock()
my_rom = romberg(f,1,2,1.0e-6)[0]
print 'My romberg integral is: ' + str(my_rom)
print 'Took ' + str(time.clock() - start) + ' seconds'
start = time.clock()
scipy_rom = rb(f,1,2,tol=1.0e-6)
print 'Scipy romberg integral is: ' + str(scipy_rom)
print 'Took ' + str(time.clock() - start) + ' seconds'
diff = abs(my_rom - scipy_rom)
tol = 1*10**-6
print 'The difference is ' + str(diff)
