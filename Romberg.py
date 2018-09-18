import math, time
from math import e
from scipy.integrate import romberg as rb

def Trap(f,a,b,k,Iold):
    '''Returns the integral of a function for the range a and b estimated by means of the trapezoidal method using k steps.'''
    H = (b - a) / (2.0)
    if k == 1:
        Inew = (f(a) + f(b)) / H / 2
        return Inew
    else:
        n = 2 ** (k - 1)
        h = (b - a) / n
        x = a + h / 2.0
        s = 0
        for i in range(n):
            s = s + f(x)
            x = x + h
        Inew = (Iold + h * s) / 2.0
        #print 'Integral of step '+str(n)+' is: '+ str(h*s)
    return s

def Richardson(R,k):
    for j in range(k-1,0,-1):
        c = 4.0**(k - j)
        R[j] = (c * R[j+1] - R [j]) / (c-1.0)
    return R

def romberg(f, a, b, tol=1.0e-6):
    '''I, n_panels = romberg(f, a, b, tol=1.0e-6)
    Return the integral from a to b of f(x), by Romberg integration,
    as well as the number of panels used'''
    R = []
    R[1] = Trap(f,a,b,0,1)
    oldR = R[1]
    running = True
    k = 1
    while running:
        k += 1
        R[k] = trapezoid(f,a,b,R[k-1],k)
        R = richardson(R,k)
        if abs(R[1] - oldR) < tol:
            return R[1],2**(k-1)
        Rold = R[1]

def f(x):
	return x**3

start = time.clock()
my_rom = romberg(f,0,1)
print 'My romberg integral is: ' + str(my_rom[0])
print 'Took ' + str(time.clock() - start) + ' seconds and '+str(my_rom[1])+' iterations'
start = time.clock()
scipy_rom = rb(f,0,1,tol=1.0e-6)
print 'Scipy romberg integral is: ' + str(scipy_rom)
print 'Took ' + str(time.clock() - start) + ' seconds'
diff = abs(my_rom[0] - scipy_rom)
tol = 1*10**-6
print 'The difference is ' + str(diff)
