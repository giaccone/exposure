import matplotlib.pylab as plt
import numpy as np
from exposure.eu_directive import eu_filter
from scipy.signal import bode
from matplotlib import rc
rc('font', **{'family': 'serif'})
rc('text', usetex=True)

# Define filter
quantity = 'E'
kind = 'health'

# frequency
domain = 'freq'
f = np.logspace(0, 6, 1000)
H_f, theta_f = eu_filter(quantity, kind, domain, f)

# time
domain = 'time'
num, den = eu_filter(quantity, kind, domain, f)
omega, Hdb, theta_t = bode((num, den), 2 * np.pi * f)
H_t = 10 ** (Hdb / 20)


# plot
plt.figure()
plt.subplot(2,1,1)
plt.title("{} - {})".format(quantity, kind), fontsize=12)
plt.loglog(f, H_f)
plt.loglog(f, H_t)
plt.ylabel('magnitude', fontsize=14)
plt.subplot(2,1,2)
plt.semilogx(f, theta_f)
plt.semilogx(f, theta_t)
plt.ylabel('phase (deg)', fontsize=14)
plt.xlabel('frequency (Hz)', fontsize=14)
plt.legend(('frequency domain', 'time domain'))
plt.tight_layout()

plt.show()