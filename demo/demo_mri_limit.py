import matplotlib.pylab as plt
import numpy as np
from exposure.mri import mri_limit
# define a frequency input
esd = np.logspace(-2,1,1000)*1e-3
# define the desired quantity
quantity = 'e'
if quantity.lower() == 'e':
    um = '(V/m)'
elif quantity.lower() == 'dbdt':
    um = '(T/s)'
# elif quantity.lower() == 'b':
#     um = '(T)'

q1 = mri_limit(esd, 'heart', quantity)
q2 = mri_limit(esd, 'PNS', quantity, mode='NOM')
q3 = mri_limit(esd, 'PNS', quantity, mode='FLOM')

hsig = plt.figure(facecolor='w')
plt.loglog(esd*1e3, q1, 'C0', linewidth=2, label='Cardiac stimulation')
plt.loglog(esd*1e3, q2, 'C1', linewidth=2, label='PNS stimulation-NOM')
plt.loglog(esd*1e3, q3, 'C2', linewidth=2, label='PNS stimulation-FLOM')
plt.xlabel('time (ms)', fontsize=14)
plt.ylabel('Threshold' + ' ' + um, fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.grid()
plt.tight_layout()

plt.show()


