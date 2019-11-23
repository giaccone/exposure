# load modules
import matplotlib.pylab as plt
from exposure.icnirp import icnirp_limit
import numpy as np
plt.ion()

# set frequency range
freq = np.logspace(0,7,1001)
B98o = icnirp_limit(freq, '1998', 'occupational', 'B')
B98p = icnirp_limit(freq, '1998', 'public', 'B')
B10o = icnirp_limit(freq, '2010', 'occupational', 'B')
B10p = icnirp_limit(freq, '2010', 'public', 'B')

E98o = 0.2*icnirp_limit(freq, '1998', 'occupational', 'J')
E98p = 0.2*icnirp_limit(freq, '1998', 'public', 'J')

Ecns10o = icnirp_limit(freq, '2010', 'occupational', 'Ecns')
Epns10o = icnirp_limit(freq, '2010', 'occupational', 'Epns')
Ecns10p = icnirp_limit(freq, '2010', 'public', 'Ecns')
Epns10p = icnirp_limit(freq, '2010', 'public', 'Epns')


fg1 = plt.figure()
plt.loglog(freq, B98o, 'C0', linewidth=2, label='1998, Occ.')
plt.loglog(freq, B98p, 'C0--', linewidth=2, label='1998, Pub.')
plt.loglog(freq, B10o, 'C1',linewidth=2, label='2010, Occ.')
plt.loglog(freq, B10p, 'C1--', linewidth=2, label='2010, Pub.')
plt.xlabel('frequency (Hz)',fontsize=16)
plt.ylabel(r'Magnetic flux density ($\mu$T)', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
plt.grid()
plt.tight_layout()


fg2 = plt.figure()
plt.loglog(freq, E98o, 'C0', linewidth=2, label='1998, Occ.')
plt.loglog(freq, E98p, 'C0--', linewidth=2, label='1998, Pub.')

plt.loglog(freq, Ecns10o, 'C1', linewidth=2, label='2010, Occ., CNS')
plt.loglog(freq, Epns10o, 'C1--', linewidth=2, label='2010, Occ., PNS')
plt.loglog(freq, Ecns10p, 'C3', linewidth=2, label='2010, Pub., CNS')
plt.loglog(freq, Epns10p, 'C3--', linewidth=2, label='2010, Pub., PNS')

plt.xlabel('frequency (Hz)', fontsize=16)
plt.ylabel('Electric field (V/m)', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(bbox_to_anchor=(0.5, 1),fontsize=12)
plt.grid()
plt.tight_layout()