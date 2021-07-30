import numpy as np

def eu_limit(f, quantity, kind):
    """
    :param f: 
    :param quantity: 
    :param kind: 
    :return: 
    
    AUTHOR: Luca Giaccone (luca.giaccone@polito.it)
    DATE: 26/05/2017
    HISTORY:

    # EXAMPLE:
    >>> f = 50 # Hz
    >>> quantity = 'B'
    >>> kind = 'low'
    >>> B = EuLimit(f, quantity, kind)
    >>> B*1e6
    1000.0
    """

    # check frequency range
    if np.logical_or(np.any(f < 1), np.any(f > 10e6)):
        print("Frequency out of range [1 10e6] Hz\n")
        return None

    if isinstance(f,(int,float)):
        f = np.array([f])

    # Initialize output
    limit = np.zeros(f.shape)

    # Get limits according to input
    if quantity.upper() == 'B':
        if kind.lower() == 'low':
            from exposure.icnirp import icnirp_limit
            limit = icnirp_limit(f, '2010', 'occupational',quantity)

        elif kind.lower() == 'high':
            limit[(f >= 1) & (f < 3000)] = 0.3 / f[(f >= 1) & (f < 3000)]
            limit[(f >= 3000) & (f <= 10e6)] = 1e-4

        elif kind.lower() == 'loc':
            limit[(f >= 1) & (f < 3000)] = 0.9 / f[(f >= 1) & (f < 3000)]
            limit[(f >= 3000) & (f <= 10e6)] = 3e-4

    elif quantity.upper() == 'E':
        if kind.lower() == 'sensory':
            limit[(f >= 1) & (f < 10)] = 0.7 / f[(f >= 1) & (f < 10)]
            limit[(f >= 10) & (f < 25)] = 0.07
            limit[(f >= 25) & (f < 400)] = 0.0028 * f[(f >= 25) & (f < 400)]
            limit[(f >= 400) & (f < 3e3)] = 1.1
            limit[(f >= 3e3) & (f <= 10e6)] = 3.8e-4 * f[(f >= 3e3) & (f <= 10e6)]

        elif kind.lower() == 'health':
            limit[(f >= 1) & (f < 3e3)] = 1.1
            limit[(f >= 3e3) & (f <= 10e6)] = 3.8e-4 * f[(f >= 3e3) & (f <= 10e6)]


    return limit


def eu_filter(quantity,kind,domain,f=None):
    """
    EuFilter computes filter parameters for the Weighted Peak Method

    :param year: reference guidelines. '1998' or '2010'
    :param receptor: 'occupational' or 'public'
    :param quantity: 'B', 'J', 'Ecns' or 'Epns' (to be defined according to the guidelines)
    :param domain: 'freq' or 'time' to build the frequency of time filter, respectively
    :param domain_var: variable related to the domain chosen. Sampling time (Ts) for time-domain
    or frequency array (f) in freqency-domain
    :param kwargs:
       * 'RCseries': option to build the filters related to the 1998 guidelines according to the 2010
           guidelines approach (series of RC filters)
    :return: parameters of the filter. (WF, phase) in frequency domain. (num, den, num, den) in time-domain

    AUTHOR: Luca Giaccone (luca.giaccone@polito.it)
    DATE: 29.02.2016
    HISTORY:

    EXAMPLE:
    >>> quantity = 'B'
    >>> kind = 'low'
    >>> Ts = 2e-4
    >>> num, den, num, den = EuFilter(quantity, kind, 'time', Ts)
    """


    # Create filter according to input
    if quantity.upper() == 'B':

        if kind.lower() == 'low':

            from exposure.icnirp import icnirp_filter
            if domain.lower() == 'freq':
                weight_fun, phase = icnirp_filter('2010','occupational',quantity, domain, f)
            elif domain.lower() == 'time':
                num, den = icnirp_filter('2010', 'occupational', quantity, domain)

        elif kind.lower() == 'high':

            if domain.lower() == 'freq':
                # get sign
                isgn = np.sign(f)
                f = np.abs(f)
                # Initialize output
                weight_fun = np.zeros(f.shape)
                phase = np.zeros(f.shape)

                weight_fun[(f >= 1) & (f < 3000)] = f[(f >= 1) & (f < 3000)] / (np.sqrt(2) * 0.3)
                weight_fun[(f >= 3000) & (f <= 10e6)] = 1 / (np.sqrt(2) * 1e-4)

                phase[f < 3000] = 90 * isgn[f < 3000]
            elif domain.lower() == 'time':
                a = 2 * np.pi * 3000

                fref = np.array([1e5])
                s = 2j * np.pi * fref
                Href = s / (s + a)
                lim_ref = eu_limit(fref, quantity, kind) * np.sqrt(2)
                k = np.asscalar(1.0 / (lim_ref * np.abs(Href)))

                num = np.array([k, 0])
                den = np.array([1, a])

        elif kind.lower() == 'loc':

            if domain.lower() == 'freq':
                # get sign
                isgn = np.sign(f)
                f = np.abs(f)
                # Initialize output
                weight_fun = np.zeros(f.shape)
                phase = np.zeros(f.shape)

                weight_fun[(f >= 1) & (f < 3000)] = f[(f >= 1) & (f < 3000)] / (np.sqrt(2) * 0.9)
                weight_fun[(f >= 3000) & (f <= 10e6)] = 1 / (np.sqrt(2) * 3e-4)

                phase[f < 3000] = 90 * isgn[f < 3000]

            elif domain.lower() == 'time':
                a = 2 * np.pi * 3000

                fref = np.array([1e5])
                s = 2j * np.pi * fref
                Href = s / (s + a)
                lim_ref = eu_limit(fref, quantity, kind) * np.sqrt(2)
                k = np.asscalar(1.0 / (lim_ref * np.abs(Href)))

                num = np.array([k, 0])
                den = np.array([1, a])

    elif quantity.upper() == 'E':

        if kind.lower() == 'sensory':

            if domain.lower() == 'freq':
                # get sign
                isgn = np.sign(f)
                f = np.abs(f)
                # Initialize output
                weight_fun = np.zeros(f.shape)
                phase = np.zeros(f.shape)

                weight_fun[(f >= 1) & (f < 10)] = f[(f >= 1) & (f < 10)] / 0.7
                weight_fun[(f >= 10) & (f < 25)] = 1 / 0.07
                weight_fun[(f >= 25) & (f < 400)] = 1 / (0.0028 * f[(f >= 25) & (f < 400)])
                weight_fun[(f >= 400) & (f < 3e3)] = 1 / 1.1
                weight_fun[(f >= 3e3) & (f <= 10e6)] = 1 / (3.8e-4 * f[(f >= 3e3) & (f < 10e6)])

                phase[(f >= 1) & (f < 10)] = 90 * isgn[(f >= 1) & (f < 10)]
                phase[(f >= 10) & (f < 25)] = 0
                phase[(f >= 25) & (f < 400)] = -90 * isgn[(f >= 25) & (f < 400)]
                phase[(f >= 400) & (f < 3000)] = 0
                phase[(f >= 3000) & (f < 10e6)] = -90 * isgn[(f >= 3000) & (f < 10e6)]

            if domain.lower() == 'time':
                a = 2 * np.pi * 10
                b = 2 * np.pi * 25
                c = 2 * np.pi * 400
                d = 2 * np.pi * 3000

                fref = np.array([30e3])
                s = 2j * np.pi * fref
                Href = (s * (s + c)) / ((s + a) * (s + b) * (s + d))
                lim_ref = eu_limit(fref, quantity, kind)
                k = np.asscalar(1.0 / (lim_ref * np.abs(Href)))

                num = [k, k * c, 0]
                den = [1, (a + b + d), (a * b + a * d + b * d), a * b * d]

        elif kind.lower() == 'health':

            if domain.lower() == 'freq':
                # get sign
                isgn = np.sign(f)
                f = np.abs(f)
                # Initialize output
                weight_fun = np.zeros(f.shape)
                phase = np.zeros(f.shape)

                weight_fun[(f >= 1) & (f < 3e3)] = 1 / 1.1
                weight_fun[(f >= 3e3) & (f <= 10e6)] = 1 / (3.8e-4 * f[(f >= 3e3) & (f < 10e6)])

                phase[(f >= 1) & (f < 3000)] = 0
                phase[(f >= 3000) & (f < 10e6)] = -90 * isgn[(f >= 3000) & (f < 10e6)]

            elif domain.lower() == 'time':

                a = 2 * np.pi * 3000

                fref = np.array([1])
                s = 2j * np.pi * fref
                Href = 1 / (s + a)
                lim_ref = eu_limit(fref, quantity, kind)
                k = np.asscalar(1.0 / (lim_ref * np.abs(Href)))

                num = [k]
                den = [1, a]

    
    # Assign outputs
    if domain == 'time':
        return num, den
    elif domain == 'freq':
        return weight_fun, phase