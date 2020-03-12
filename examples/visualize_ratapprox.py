#!/usr/bin/python3
import numpy as np
from sys import argv,exit
from matplotlib import pyplot as plt

try:
    f = open(argv[1]).read()
    exponent = float(argv[2])
    minmass = float(argv[3])
except:
    print("Usage: {} remez_file exponent min_mass".format(argv[0]))
    exit(1)

lines = [ l for l in f.split('\n') if len(l) is not 0 ] 
direct,inverse  = lines[:len(lines)//2] , lines[len(lines)//2:]

def get_ratapprox(lines):
    a0 = float(lines[0])
    coeffs = np.array([[float(x) for x in l.split(',') ] for l in lines[1:]])
    def ratapprox(a0,coeffs,x):
        return a0 + sum([ a/(x+b) for a,b in coeffs])
    return lambda x : ratapprox(a0,coeffs,x)

x = np.exp(np.log(10)*np.arange(-10,2,0.2))


# direct
direct_ratapprox = get_ratapprox(direct)(x)
direct_true = np.power(x,exponent)

diff = np.abs(direct_ratapprox - direct_true)
squaresum = np.hypot(direct_ratapprox,direct_true)

plt.title("exponent {}".format(exponent))
plt.xscale('log')
plt.yscale('log')
plt.ylabel('diff/squaresum')
plt.plot(x, diff/squaresum)
plt.plot([minmass**2,minmass**2],[-1,1])

plt.show()
# inverse
inverse_ratapprox = get_ratapprox(inverse)(x)
inverse_true = np.power(x,-exponent)

diff = np.abs(inverse_ratapprox - inverse_true)
squaresum = np.hypot(inverse_ratapprox,inverse_true)

plt.title("exponent {}".format(-exponent))
plt.xscale('log')
plt.yscale('log')
plt.plot(x, diff/squaresum)
plt.ylabel('diff/squaresum')
plt.plot([minmass**2,minmass**2],[-1,1])
plt.show()





