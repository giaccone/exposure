import matplotlib.pylab as plt
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)
from exposure.eu_directive import eu_limit

# ----------------
# define limit
# ---------------
freq = np.logspace(0,7,100001)

Blow = eu_limit(freq, 'B', 'low')
Bhigh = eu_limit(freq, 'B', 'high')
Bloc = eu_limit(freq, 'B', 'loc')

Es = eu_limit(freq, 'E', 'sensory')
Eh = eu_limit(freq, 'E', 'health')

# ----------------
# plot
# ---------------
fg1 = plt.figure(facecolor='w')
plt.loglog(freq, Blow, linewidth=2, label='AL low')
plt.loglog(freq, Bhigh, linestyle='--' ,linewidth=2, label='AL high')
plt.loglog(freq, Bloc, linewidth=2, label='AL loc.')
plt.xlabel('frequency (Hz)',fontsize=14)
plt.ylabel(r'Magnetic flux density ($\mu$T)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
#plt.grid()
fg1.tight_layout()

fg2 = plt.figure(facecolor='w')
plt.loglog(freq, Es, linewidth=2, label='E sensory')
plt.loglog(freq, Eh, linestyle='--', linewidth=2, label='E health')
plt.xlabel('frequency (Hz)', fontsize=14)
plt.ylabel(r'Electric field (V/m)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
#plt.grid()
fg2.tight_layout()


plt.show()


