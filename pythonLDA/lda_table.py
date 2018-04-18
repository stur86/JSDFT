#!/usr/bin/python

import numpy as np
from scipy.special import hyp2f1
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

import json

inf = np.inf

# Define the gamma parameters

# gLDA parameters
g_c = {
0: [0.025979, 0.025979, 0.033891, 0.642367, -0.35379],
1: [33.0265,  0.896251, 24.2518,  16.1820,  -12.5392],
inf: [0.163723, 0.163723, 0.301135, 0.661217, 0.152167]
}

"""
# revLDA parameters
g_c = {
0: [0.025873, 0.025873, 0.032541, 0.741760, -0.498560],
1: [18.3407, -0.154372, 13.2193,  8.807757, -6.681718],
inf: [0.164037, 0.164037, 0.261152, 0.519097, 0.055756]
}
"""

# Define the gamma function

def gamma(eta, gtype):
	global g_c
	sqrtfac = np.sqrt(1.0-eta)
	return (g_c[gtype][0]-g_c[gtype][1]*sqrtfac-g_c[gtype][2]*eta)/(g_c[gtype][3]+sqrtfac+g_c[gtype][4]*eta)

# Define the GLDA energy

def eLDA (rs, eta, gamma0, gamma1, gammainf):

	return gamma0*hyp2f1(1.0, 1.5, gamma1, 2.0*gamma0*(1.0-gamma1)/gammainf*rs)

# Tests

print "Test 1 - gamma0(0): {0} == {1}".format(gamma(0.0, 0), 0)
print "Test 2 - gammainf(0): {0} == {1}".format(gamma(0.0, inf), 0)
print "Test 3 - gamma0(1): {0} == {1}".format(gamma(1.0, 0), -np.pi**2/360.0)
print "Test 4 - gammainf(1): {0} == {1}".format(gamma(1.0, inf), np.log(np.sqrt(2.0*np.pi))-0.75)

# Fitting function

def fitf(x, A, B1, t1, B2, t2, B3, t3):
	return A + B1*np.exp(-x/t1) + B2*np.exp(-x/t2) +  B3*np.exp(-x/t3)

# Calculate everything over a range

# Box width and number of electrons

plt.ion()

allpars = {}

for n in range(2, 20):
	eta = 1.0 - 1.0/n**2
	gamma0 = gamma(eta, 0)
	gamma1 = gamma(eta, 1)
	gammainf = gamma(eta, inf)
	outf = open("LDA_{0}.dat".format(n), 'w')
	xdata = np.array(10**(np.linspace(-3,3, 100)))*2*n
	ydata = []
	for L in xdata:
		ec = eLDA(L/(2*n), eta, gamma0, gamma1, gammainf)
		ydata.append(ec*1000)
		outf.write("{0}\t{1}\n".format(L/np.pi, ec*1000))
	ydata = np.array(ydata)
	fitsol = curve_fit(fitf, xdata, ydata, p0=[ydata[-1], (ydata[0]-ydata[-1])*0.2, 1000.0, (ydata[0]-ydata[-1])*0.2, 100.0, (ydata[0]-ydata[-1])*0.6, 20.0])

	allpars[n] = list(fitsol[0])

	plt.clf()
	plt.plot(np.log10(xdata), ydata, 'o')
	plt.plot(np.log10(xdata), fitf(xdata, *fitsol[0]))

json.dump(allpars, open('LDA_3exp_pars.json', 'w'))

