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

def Richard(k,R,Rk):
	outlist = []
	for j in range(k-1,0,-1):
		if j + 1 == k:
			output = (4 ** (k - j) * Rk - R[j-1]) / (4 ** (k - j) - 1)
		else:
			output = (4 ** (k - j) * R[j] - R[j-1]) / (4 ** (k - j) - 1)
		outlist.insert(0,output)
	return outlist

def romberg(f, a, b, tol=1.0e-6):
    '''I, n_panels = romberg(f, a, b, tol=1.0e-6)
    Return the integral from a to b of f(x), by Romberg integration,
    as well as the number of panels used'''
    k = 1
    R = [Trap(f, a, b, 1), Trap(f, a, b, 2)]
    running = True
    while running:
    	k += 1
    	Rk = R[0]
    	R = Richard(k,R,Rk)
    	E = abs(Rk - R[0])
    	print 'The error is '+str(E)
    	if E < tol:
    		running = False
    return R[0],k

def f(x):
	return 3*x**2 - 4*x**1 + 5

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
