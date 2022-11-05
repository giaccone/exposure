from exposure.icnirp import icnirp_limit, icnirp_filter
from exposure.modified.icnirp import icnirp_filter as icnirp_filter_beta
import numpy as np
from fourier.ffttools import fft_analysis
from exposure.pulsed import sum_rule, wpm_time, wpm_freq
import matplotlib.pyplot as plt
plt.ion()

# Define filter
year = '2010'
receptor = 'occupational'
quantity = 'B'

# Define time domain waveform
wave = 'trapz'
if wave == 'sine_pulse':
    N = 9999
    fr = 1000
    T = 1/fr
    winT = 3*T
    dt = winT / N
    t = np.linspace(0, winT - dt, N)
    id1 = int(N/3)
    id2 = int(id1 + N/3)
    Glim = icnirp_limit(fr, year, receptor, quantity)
    B = np.zeros((t.size, 1))
    B[id1:id2, 0] = np.sqrt(2) * Glim * np.sin(2 * np.pi * fr * t[id1:id2])

elif wave == 'trapz':
    N = 9999
    Tr = 2.5e-3
    Tf = 10e-3
    Tw = 100e-3
    T = 300e-3
    t = np.linspace(0, T, N)
    Glim = icnirp_limit(1/2/Tr, year, receptor, quantity)
    B = np.zeros((t.size, 1))
    B[t < Tr, 0] = Glim * np.sqrt(2) / Tr * t[t < Tr]
    B[(t >= Tr) & (t < Tr + Tw)] = Glim * np.sqrt(2)
    Tstart = t[t < Tr + Tw][-1]
    B[(t >= Tr + Tw) & (t < Tr + Tw + Tf), 0] = Glim * np.sqrt(2) - Glim * np.sqrt(2) / Tf * (t[(t >= Tr + Tw) & (t < Tr + Tw + Tf)] - Tstart)



# Evaluate spectrum
f, BFT = fft_analysis(t, B)

# Frequency summation rule
Blim = icnirp_limit(np.abs(f), year, receptor, quantity) * np.sqrt(2)
IS, ISv = sum_rule(np.abs(BFT), Blim)

# WP in frequency domain
WF, phi = icnirp_filter(year, receptor, quantity, 'freq', f)
IWf, WPf = wpm_freq(WF, phi, f, BFT, makeplot='n')

# WP in frequency domain (beta)
WFb, phib = icnirp_filter_beta(year, receptor, quantity, 'freq', f)
IWfb, WPfb = wpm_freq(WFb, phib, f, BFT, makeplot='n')

# WP in time domain
num, den = icnirp_filter(year, receptor, quantity, 'time')
IW, WP = wpm_time(num, den, t, B, makeplot='n')

print("Frequency summation rule: {:.4f}".format(IS))
print("Weighted peak method (freq): {:.4f}".format(IWf))
print("Weighted peak method (freq, new): {:.4f}".format(IWfb))
print("Weighted peak method (time): {:.4f}".format(IW))

# final plot
hsig = plt.figure(facecolor='w')
plt.plot(t, B*1e6, 'C0', linewidth=2, label='B-field ($\mu T$)')
plt.xlabel('time (s)', fontsize=14)
plt.ylabel('signal', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.grid()
plt.tight_layout()


hWP = plt.figure(facecolor='w')
plt.plot(t, np.sqrt(np.sum(WPf ** 2, axis=1)), 'C0', linewidth=2, label='frequency domain')
plt.plot(t, IWf * np.ones(t.shape), 'C0--', linewidth=1)

plt.plot(t, np.sqrt(np.sum(WP ** 2, axis=1)), 'C1', linewidth=2, label='time domain')
plt.plot(t, IW * np.ones(t.shape), 'C1--', linewidth=1)

# freq new
plt.plot(t, np.sqrt(np.sum(WPfb ** 2, axis=1)), 'k--', linewidth=2, label='frequency domain (proposed)')
plt.plot(t, IWfb * np.ones(t.shape), 'k--', linewidth=1)
#
plt.xlabel('time (s)',fontsize=14)
plt.ylabel('weighted signal', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.grid()
#plt.title('$EI_t$ = {:.2f}    -    $EI_f$ = {:.2f}'.format(IW, IWf))
plt.tight_layout()
