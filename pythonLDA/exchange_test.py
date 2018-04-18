import numpy as np
import scipy.constants as cnst
from scipy.special import sici

import matplotlib.pyplot as plt

V = []

gamma = 0.5772156649

p = 2

for K in range(1, 100):

	Vp = (np.log(p*(2.0*K+1.0)*np.pi) + gamma - sici(p*(2.0*K+1.0)*np.pi)[1])/np.pi
	V.append(Vp)

plt.ion()

plt.plot(V, 'o')
plt.plot( [(np.log(2.0*np.pi*p) + gamma)/np.pi]*len(V), '-')
plt.show()

raw_input()



