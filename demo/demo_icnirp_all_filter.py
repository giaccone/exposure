# load modules
from exposure.icnirp import icnirp_filter
from scipy.signal import bode
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# all possible parameters
y = ['1998', '2010']
r = ['occupational', 'public']
q = ['B', 'J', 'Ecns', 'Epns']

# cycle on all configurations
for year in y:
    for receptor in r:
        for quantity in q:
            if (year == '1998') and ((quantity == 'Epns') or (quantity == 'Ecns')):
                skip = True
            elif (year == '2010') and (quantity == 'J'):
                skip = True
            else:
                skip = False

            if not skip:
                print("year {}, receptor {}, quantity {}".format(year, receptor, quantity))

                # frequency
                domain = 'freq'
                f = np.logspace(0, 6, 1000)
                H_f, theta_f = icnirp_filter(year, receptor, quantity, domain, f)

                # time
                domain = 'time'
                num, den = icnirp_filter(year, receptor, quantity, domain)
                omega, Hdb, theta_t = bode((num, den), 2 * np.pi * f)
                H_t = 10 ** (Hdb / 20)

                # plot
                plt.figure()
                plt.subplot(2,1,1)
                plt.title("{} - {} (ICNIRP {})".format(receptor, quantity, year), fontsize=12)
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
