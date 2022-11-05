import numpy as np
import scipy.signal as sig
import fourier.ffttools as ft
import matplotlib.pylab as plt


def sum_rule(A, Alim):
    """

    Parameters
    ----------
    A (ndarray):
        spectrum, only the magnitude, of the variable A
    Alim (ndarray):
        limit at each frequency of the spectrum

    Returns
    -------
        EI (float):
            exposure index
        EI_array (float):
            array with the exposure index at each frequency

    """

    # local copy of A
    Aloc = np.copy(A)

    # force Alim and Aloc to be a column vector
    Alim.shape = (-1, 1)
    dim = Aloc.ndim
    if dim == 1:
        Aloc = np.reshape(A, (Aloc.size, -1))
    if dim > 1:
        (r, c) = Aloc.shape
        if r < c:
            Aloc = Aloc.transpose()

    # protection against possible division by zero
    iz = np.where(Alim == 0)
    Aloc[iz] = 0
    Alim[iz] = 1

    ratio = Aloc/Alim
    EI_array = np.sum(ratio, axis=0)
    EI = np.sqrt(np.sum(EI_array ** 2))

    return EI, EI_array


def wpm_time(num, den, t, A, makeplot='y'):
    """
    WPM_TIME performs the Weighted Peak Method in time domain

    Parameters
    ----------
    num (ndarray):
        numerator of the filter
    den (ndarray):
        denominator of the filter
    t (ndarray):
        array with time values. t.shape = (N,)
    A (ndarray):
        array with the quantity ti be weighted. A.shape = (N, Ncomp)
    makeplot (str):
        'y' or 'n' to plot or not plot results, respectively (default 'n')

    Returns
    -------
    EI (float):
        weighted peak index (scalar value)
    WP (ndarray)
        weighted waveforms. WP.shape = (N, Ncomp)
    hfg (figure):
        handle to figures

    AUTHOR: Luca Giaccone (luca.giaccone@polito.it)
    DATE: 02.03.2016
    HISTORY:
    """

    # force A to be a column vector
    dim = A.ndim
    if dim == 1:
        A = np.reshape(A,(A.size,-1))

    if dim > 1:
        (r, c) = A.shape
        if r < c:
            A = A.transpose()

    WP = np.zeros_like(A)

    for k, row in enumerate(A.transpose()):
        _, WP[:, k], _ = sig.lsim((num, den), row, t)

    WPmod = np.sqrt(np.sum(WP ** 2,axis=1))
    EI = np.max(WPmod)

    hfg = [0,0,0]
    if makeplot == 'y':
        nsp = A.shape[1]*100 + 10
        hfg[0] = plt.figure(facecolor='w')
        for k in range(A.shape[1]):
            nsp += 1
            ax = hfg[0].add_subplot(nsp)
            ax.plot(t, A[:, k], 'C0',linewidth=2)
            plt.ylabel('comp #{}'.format(k), fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.grid()
            plt.tight_layout()
        else:
            plt.xlabel('time (s)', fontsize=14)
            plt.tight_layout()

        nsp = A.shape[1] * 100 + 10
        hfg[1] = plt.figure(facecolor='w')
        for k in range(A.shape[1]):
            nsp += 1
            ax = hfg[1].add_subplot(nsp)
            ax.plot(t, WP[:, k], 'C2', linewidth=2, label='comp ' + str(k))
            plt.ylabel('Wcomp #{}'.format(k), fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.grid()
            plt.tight_layout()
        else:
            plt.xlabel('time (s)', fontsize=14)
            plt.tight_layout()


        hfg[2] = plt.figure(facecolor='w')
        plt.plot(t, WPmod, 'C3',linewidth=2, label='module')
        plt.plot(t, EI * np.ones(t.shape), 'k--', linewidth=1)
        plt.title('IW = {:.2f}'.format(EI),fontsize=14)
        plt.xlabel('time (s)', fontsize=14)
        plt.ylabel('WPmod', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(np.linspace(0, EI, 5),fontsize=14)
        plt.legend()
        plt.grid()
        plt.tight_layout()

    if makeplot == 'y':
        return EI, WP, hfg
    else:
        return EI, WP


def wpm_freq(weight_fun, phase, f, A, **kwargs):
    """
    WPM_FREQ performs the Weighted Peak Method in time domain

    Parameters
    ----------
    WF (ndarray):
        module of the weighting function
    phi (ndarray):
        phase of the weighting function
    f (ndarray):
        frequency array (ouput of ft.fftanalysys)
    A (ndarray):
        spectrum of the signal (ouput of ft.fftanalysys)
    makeplot (str):
        'y' or 'n' to plot or not plot results, respectively (default 'n')

    Returns
    -------
    EI (float):
        weighted peak index (scalar value)
    WP (ndarray)
        weighted waveforms. WP.shape = (N, Ncomp)
    hfg (figure):
        handle to figures

    AUTHOR: Luca Giaccone (luca.giaccone@polito.it)
    DATE: 07.03.2016
    HISTORY:
    """

    # defaults
    makeplot = 'n'

    # Check optiona input
    for key in kwargs:
        if key.lower() == 'makeplot':
            makeplot = kwargs[key]

    # force A to be a column vector
    dim = A.ndim
    if dim == 1:
        A = np.reshape(A, (A.size, -1))
    if dim > 1:
        (r, c) = A.shape
        if r < c:
            A = A.transpose()

    weight_fun.shape = (weight_fun.size, -1)
    phase.shape = (phase.size, -1)


    C_WF = np.zeros(weight_fun.shape, dtype=complex)
    C_WF.real = weight_fun * np.cos(phase * np.pi / 180)
    C_WF.imag = weight_fun * np.sin(phase * np.pi / 180)

    #C_WF.shape = A.shape
    WP = ft.inverse_fft(A * C_WF)
    WPmod = np.sqrt(np.sum(WP ** 2, axis=1))
    EI = np.max(WPmod)

    hfg = [0, 0, 0]
    if makeplot == 'y':
        df = f[1] - f[0]
        T = 1/df
        t = np.linspace(0,T,A.shape[0])
        a = ft.inverse_fft(A)

        nsp = A.shape[1] * 100 + 10
        hfg[0] = plt.figure(facecolor='w')
        for k in range(A.shape[1]):
            nsp += 1
            ax = hfg[0].add_subplot(nsp)
            ax.plot(t, a[:, k], 'k', linewidth=3)
            plt.ylabel('comp #{}'.format(k), fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.grid()
            plt.tight_layout()
        else:
            plt.xlabel('time (s)', fontsize=14)
            plt.tight_layout()

        nsp = A.shape[1] * 100 + 10
        hfg[1] = plt.figure(facecolor='w')
        for k in range(A.shape[1]):
            nsp += 1
            ax = hfg[1].add_subplot(nsp)
            ax.plot(t, WP[:, k], 'b', linewidth=3, label='comp ' + str(k))
            plt.ylabel('Wcomp #{}'.format(k), fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.grid()
            plt.tight_layout()
        else:
            plt.xlabel('time (s)', fontsize=14)
            plt.tight_layout()

        hfg[2] = plt.figure(facecolor='w')
        plt.plot(t, WPmod, 'r', linewidth=3, label='module')
        plt.plot(t, EI * np.ones(t.shape), 'k--', linewidth=1)
        plt.title('IW = {:.2f}'.format(EI), fontsize=14)
        plt.xlabel('time (s)', fontsize=14)
        plt.ylabel('WPmod', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(np.linspace(0, EI, 5), fontsize=14)
        plt.legend()
        plt.grid()
        plt.tight_layout()

    if makeplot == 'y':
        return EI, WP, hfg
    else:
        return EI, WP
