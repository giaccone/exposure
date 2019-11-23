# load modules
import numpy as np


def icnirp_limit(f, year, receptor, quantity):
    """
    INCIRP_LIMIT provides basic restrictions or reference levels according to the ICNIRP guidelines.

    Parameters
    ----------
    f (float, ndarray):
        frequency range
    year (str):
        reference guidelines. '1998' or '2010'
    receptor (str):
        'occupational' or 'public'
    quantity (str):
        'B', 'J', 'Ecns' or 'Epns' (to be defined according to the guidelines)

    Returns
    -------
    return (float, ndarray):
        limit (N.B. always in S.I. unit)

    AUTHOR: Luca Giaccone (luca.giaccone@polito.it)
    DATE: 23.02.2016
    HISTORY:
    """

    if isinstance(f,(int,float)):
        f = np.array([f])

    # Initialize output
    limit = np.zeros(f.shape)

    # Get limits according to input
    if year == '1998':
        if receptor == 'occupational':
            if quantity == 'B':
                limit[f < 1] = 0.2
                limit[(f >= 1) & (f < 8)] = 0.2/f[(f >= 1) & (f < 8)]**2
                limit[(f >= 8) & (f < 25)] = 0.025/f[(f >= 8) & (f < 25)]
                limit[(f >= 25) & (f < 820)] = 25e-6/(f[(f >= 25) & (f < 820)]*1e-3)
                limit[(f >= 820) & (f < 65e3)] = 30.7e-6
                limit[(f >= 65e3) & (f < 1e6)] = 2e-6/(f[(f >= 65e3) & (f < 1e6)]*1e-6)
                limit[(f >= 1e6) & (f < 10e6)] = 2e-6/(f[(f >= 1e6) & (f < 10e6)]*1e-6)
                limit[(f >= 10e6) & (f < 400e6)] = 0.2e-6
                limit[(f >= 400e6) & (f < 2000e6)] = 0.01e-6*np.sqrt(f[(f >= 400e6) & (f < 2000e6)]*1e-6)
                limit[(f >= 2000e6) & (f <= 300e9)] = 0.45e-6

            elif quantity == 'J':
                limit[f < 1] = 40 * 1e-3
                limit[(f >= 1) & (f < 4)] = 40*1e-3/f[(f >= 1) & (f < 4)]
                limit[(f >= 4) & (f < 1000)] = 10*1e-3
                limit[(f >= 1000) & (f < 100e3)] = f[(f >= 1000) & (f < 100e3)]/100*1e-3
                limit[(f >= 100e3) & (f <= 10e6)] = f[(f >= 100e3) & (f <= 10e6)]/100 * 1e-3

        elif receptor == 'public':
            if quantity == 'B':
                limit[f < 1] = 4e-2
                limit[(f >= 1) & (f < 8)] = 4e-2 / f[(f >= 1) & (f < 8)] ** 2
                limit[(f >= 8) & (f < 25)] = 5e-3 / f[(f >= 8) & (f < 25)]
                limit[(f >= 25) & (f < 800)] = 5e-6 / (f[(f >= 25) & (f < 800)] * 1e-3)
                limit[(f >= 800) & (f < 3000)] = 6.25e-6
                limit[(f >= 3000) & (f < 150e3)] = 6.25e-6
                limit[(f >= 150e3) & (f < 1e6)] = 0.92e-6 / (f[(f >= 150e3) & (f < 1e6)] * 1e-6)
                limit[(f >= 1e6) & (f < 10e6)] = 0.92e-6 / (f[(f >= 1e6) & (f < 10e6)] * 1e-6)
                limit[(f >= 10e6) & (f < 400e6)] = 0.092e-6
                limit[(f >= 400e6) & (f < 2000e6)] = 0.0046e-6 * np.sqrt(f[(f >= 400e6) & (f < 2000e6)] * 1e-6)
                limit[(f >= 2000e6) & (f <= 300e9)] = 0.2e-6

            elif quantity == 'J':
                limit[f < 1] = 8*1e-3
                limit[(f >= 1) & (f < 4)] = 8*1e-3 / f[(f >= 1) & (f < 4)]
                limit[(f >= 4) & (f < 1000)] = 2*1e-3
                limit[(f >= 1000) & (f < 100e3)] = f[(f >= 1000) & (f < 100e3)] / 500 * 1e-3
                limit[(f >= 100e3) & (f <= 10e6)] = f[(f >= 100e3) & (f <= 10e6)] / 500 * 1e-3

    elif year == '2010':
        if receptor == 'occupational':
            if quantity == 'B':
                limit[(f >= 1) & (f < 8)] = 0.2 / f[(f >= 1) & (f < 8)] ** 2
                limit[(f >= 8) & (f < 25)] = 2.5e-2 / f[(f >= 8) & (f < 25)]
                limit[(f >= 25) & (f < 300)] = 1e-3
                limit[(f >= 300) & (f < 3000)] = 0.3 / f[(f >= 300) & (f < 3000)]
                limit[(f >= 3000) & (f <= 10e6)] = 1e-4

            elif quantity == 'Ecns':
                limit[(f >= 1) & (f < 10)] = 0.5 / f[(f >= 1) & (f < 10)]
                limit[(f >= 10) & (f < 25)] = 0.05
                limit[(f >= 25) & (f < 400)] = 2e-3 * f[(f >= 25) & (f < 400)]
                limit[(f >= 400) & (f < 3000)] = 0.8
                limit[(f >= 3000) & (f <= 10e6)] = 2.7e-4 * f[(f >= 3000) & (f <= 10e6)]

            elif quantity == 'Epns':
                limit[(f >= 1) & (f < 3000)] = 0.8
                limit[(f >= 3000) & (f <= 10e6)] = 2.7e-4 * f[(f >= 3000) & (f <= 10e6)]

        elif receptor == 'public':
            if quantity == 'B':
                limit[(f >= 1) & (f < 8)] = 4e-2 / f[(f >= 1) & (f < 8)] ** 2
                limit[(f >= 8) & (f < 25)] = 5e-3 / f[(f >= 8) & (f < 25)]
                limit[(f >= 25) & (f < 50)] = 2e-4
                limit[(f >= 50) & (f < 400)] = 2e-4
                limit[(f >= 400) & (f < 3000)] = 8e-2 / f[(f >= 400) & (f < 3000)]
                limit[(f >= 3000) & (f <= 10e6)] = 2.7e-5

            elif quantity == 'Ecns':
                limit[(f >= 1) & (f < 10)] = 0.1 / f[(f >= 1) & (f < 10)]
                limit[(f >= 10) & (f < 25)] = 0.01
                limit[(f >= 25) & (f < 1000)] = 4e-4 * f[(f >= 25) & (f < 1000)]
                limit[(f >= 1000) & (f < 3000)] = 0.4
                limit[(f >= 3000) & (f <= 10e6)] = 1.35e-4 * f[(f >= 3000) & (f <= 10e6)]

            elif quantity == 'Epns':
                limit[(f >= 1) & (f < 3000)] = 0.4
                limit[(f >= 3000) & (f <= 10e6)] = 1.35e-4 * f[(f >= 3000) & (f <= 10e6)]

    if limit.size == 1:
        limit = np.ndarray.item(limit)

    return limit


def icnirp_filter(year, receptor, quantity, domain, f=None, rc_series=None):
    """

    Parameters
    ----------
    year (str):
        ICNIRP guidelines publication yeare, '1998' or '2010'
    receptor (str):
        'occupational' or 'public'
    quantity (str):
        string defining the phisical quantity (e.g. 'B', 'J', 'Ecns', 'Epns')
    domain (str):
        'freq' or 'time'
    f (ndarray):
        input required when domain='freq'. It defines the frequency values where
        the filter has to be defined
    rc_series (ndarray):
        optional input that can be used with year='1998'. It takes into
        account also the filter variotion (magnitude/phase) at extremely
        low frequency

    Returns
    -------
    num (ndarray):
        numerator of the filter
    den (ndarray):
        denominator of the filter

    AUTHOR: Luca Giaccone (luca.giaccone@polito.it)
    DATE: 23.11.2019
    HISTORY:
    """

    if year == '1998':
        if receptor == 'occupational':
            if quantity == 'B':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[f < 1] = 1 / (0.2 * np.sqrt(2))
                    weight_fun[(f >= 1) & (f < 8)] = 1 / (0.2 / f[(f >= 1) & (f < 8)] ** 2 * np.sqrt(2))
                    weight_fun[(f >= 8) & (f < 25)] = 1 / (0.025 / f[(f >= 8) & (f < 25)] * np.sqrt(2))
                    weight_fun[(f >= 25) & (f < 820)] = 1 / (25e-6 / (f[(f >= 25) & (f < 820)] * 1e-3) * np.sqrt(2))
                    weight_fun[f >= 820] = 1 / (30.7e-6 * np.sqrt(2))

                    phase[f < 820] = 90 * isgn[f < 820]

                elif domain == 'time':
                    if rc_series == 'y':
                        # angular frequencies
                        a = 2 * np.pi * 8
                        b = 2 * np.pi * 820

                        fref = np.array([1e4])
                        s = 2j * np.pi * fref
                        Href = (s ** 2) / ((s + a) * (s + b))
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k, 0, 0])
                        den = np.array([1, (a + b), (a * b)])

                    else:
                        # angular frequency
                        a = 2 * np.pi * 820

                        # filter parameters
                        fref = np.array([1e4])
                        s = 2j * np.pi * fref
                        Href = s / (s + a)
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k, 0])
                        den = np.array([1, a])

            elif quantity == 'J':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[f < 1] = 1 / (40 * 1e-3 * np.sqrt(2))
                    weight_fun[(f >= 1) & (f < 4)] = 1 / (40 * 1e-3 / f[(f >= 1) & (f < 4)] * np.sqrt(2))
                    weight_fun[(f >= 4) & (f < 1000)] = 1 / (10 * 1e-3 * np.sqrt(2))
                    weight_fun[(f >= 1000) & (f < 100e3)] = 1 / (f[(f >= 1000) & (f < 100e3)] / 100 * 1e-3 * np.sqrt(2))
                    weight_fun[(f >= 100e3) & (f < 10e6)] = 1 / (f[(f >= 100e3) & (f < 10e6)] / 100 * 1e-3 * np.sqrt(2))

                    phase[f > 1000] = -90 * isgn[f > 1000]
                elif domain == 'time':
                    if rc_series == 'y':
                        # angular frequencies
                        a = 2 * np.pi * 1
                        b = 2 * np.pi * 4
                        c = 2 * np.pi * 1000

                        # filter parameters
                        fref = np.array([1e5])
                        s = 2j * np.pi * fref
                        Href = (s + a) / ((s + b) * (s + c))
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k, (k * a)])
                        den = np.array([1, (b + c), (b * c)])

                    else:
                        # angular frequency
                        a = 2 * np.pi * 1000

                        # filter parameters
                        fref = np.array([1e5])
                        s = 2j * np.pi * fref
                        Href = 1 / (s + a)
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k])
                        den = np.array([1, a])

        elif receptor == 'public':
            if quantity == 'B':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[f < 1] = 1 / (4e-2 * np.sqrt(2))
                    weight_fun[(f >= 1) & (f < 8)] = 1 / (4e-2 / f[(f >= 1) & (f < 8)] ** 2 * np.sqrt(2))
                    weight_fun[(f >= 8) & (f < 25)] = 1 / (5e-3 / f[(f >= 8) & (f < 25)] * np.sqrt(2))
                    weight_fun[(f >= 25) & (f < 800)] = 1 / (5e-6 / (f[(f >= 25) & (f < 800)] * 1e-3) * np.sqrt(2))
                    weight_fun[f >= 800] = 1 / (6.25e-6 * np.sqrt(2))

                    phase[f < 800] = 90 * isgn[f < 800]

                elif domain == 'time':
                    if rc_series == 'y':
                        # angular frequencies
                        a = 2 * np.pi * 8
                        b = 2 * np.pi * 800

                        # filter parameters
                        fref = np.array([1e4])
                        s = 2j * np.pi * fref
                        Href = (s ** 2) / ((s + a) * (s + b))
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k, 0, 0])
                        den = np.array([1, (a + b), (a * b)])

                    else:
                        # angular frequency
                        a = 2 * np.pi * 820

                        # filter parameters
                        fref = np.array([1e4])
                        s = 2j * np.pi * fref
                        Href = s / (s + a)
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k, 0])
                        den = np.array([1, a])

            elif quantity == 'J':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[f < 1] = 1 / (8 * 1e-3 * np.sqrt(2))
                    weight_fun[(f >= 1) & (f < 4)] = 1 / (8 * 1e-3 / f[(f >= 1) & (f < 4)] * np.sqrt(2))
                    weight_fun[(f >= 4) & (f < 1000)] = 1 / (2 * 1e-3 * np.sqrt(2))
                    weight_fun[(f >= 1000) & (f < 100e3)] = 1 / (f[(f >= 1000) & (f < 100e3)] / 500 * 1e-3 * np.sqrt(2))
                    weight_fun[(f >= 100e3) & (f < 10e6)] = 1 / (f[(f >= 100e3) & (f < 10e6)] / 500 * 1e-3 * np.sqrt(2))

                    phase[f > 1000] = -90 * isgn[f > 1000]

                elif domain == 'time':
                    if rc_series == 'y':
                        # angular frequencies
                        a = 2 * np.pi * 1
                        b = 2 * np.pi * 4
                        c = 2 * np.pi * 1000

                        # filter parameters
                        fref = np.array([1e5])
                        s = 2j * np.pi * fref
                        Href = (s + a) / ((s + b) * (s + c))
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k, (k * a)])
                        den = np.array([1, (b + c), (b * c)])

                    else:
                        # angular frequencies
                        a = 2 * np.pi * 1000

                        # filter parameters
                        fref = np.array([1e5])
                        s = 2j * np.pi * fref
                        Href = 1 / (s + a)
                        lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                        k = (1.0 / (lim_ref * np.abs(Href))).item()

                        # define numerator and denominator
                        num = np.array([k])
                        den = np.array([1, a])

    elif year == '2010':
        if receptor == 'occupational':
            if quantity == 'B':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[(f >= 1) & (f < 8)] = 1 / (0.2 / f[(f >= 1) & (f < 8)] ** 2 * np.sqrt(2))
                    weight_fun[(f >= 8) & (f < 25)] = 1 / (2.5e-2 / f[(f >= 8) & (f < 25)] * np.sqrt(2))
                    weight_fun[(f >= 25) & (f < 300)] = 1 / (1e-3 * np.sqrt(2))
                    weight_fun[(f >= 300) & (f < 3000)] = 1 / (0.3 / f[(f >= 300) & (f < 3000)] * np.sqrt(2))
                    weight_fun[(f >= 3000) & (f < 10e6)] = 1 / (1e-4 * np.sqrt(2))

                    phase[(f >= 1) & (f < 8)] = isgn[(f >= 1) & (f < 8)] * 180
                    phase[(f >= 8) & (f < 25)] = isgn[(f >= 8) & (f < 25)] * 90
                    phase[(f >= 25) & (f < 300)] = 0
                    phase[(f >= 300) & (f < 3000)] = isgn[(f >= 300) & (f < 3000)] * 90
                    phase[(f >= 3000) & (f < 10e6)] = 0

                elif domain == 'time':
                    # angular frequencies
                    a = 2 * np.pi * 8
                    b = 2 * np.pi * 25
                    c = 2 * np.pi * 300
                    d = 2 * np.pi * 3000

                    # filter parameters
                    fref = np.array([1])
                    s = 2j * np.pi * fref
                    Href = (s ** 2. * (s + c)) / ((s + a) * (s + b) * (s + d))
                    lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                    k = (1.0 / (lim_ref * np.abs(Href))).item()

                    # define numerator and denominator
                    num = np.array([k, k * c, 0, 0])
                    den = np.array([1, (a + b + d), (a * b + a * d + b * d), (a * b * d)])

            elif quantity == 'Ecns':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[(f >= 1) & (f < 10)] = 1 / (0.5 / f[(f >= 1) & (f < 10)] * np.sqrt(2))
                    weight_fun[(f >= 10) & (f < 25)] = 1 / (0.05 * np.sqrt(2))
                    weight_fun[(f >= 25) & (f < 400)] = 1 / (2e-3 * f[(f >= 25) & (f < 400)] * np.sqrt(2))
                    weight_fun[(f >= 400) & (f < 3000)] = 1 / (0.8 * np.sqrt(2))
                    weight_fun[(f >= 3000) & (f < 10e6)] = 1 / (2.7e-4 * f[(f >= 3000) & (f < 10e6)] * np.sqrt(2))

                    phase[(f >= 1) & (f < 10)] = 90 * isgn[(f >= 1) & (f < 10)]
                    phase[(f >= 10) & (f < 25)] = 0
                    phase[(f >= 25) & (f < 1000)] = -90 * isgn[(f >= 25) & (f < 1000)]
                    phase[(f >= 1000) & (f < 3000)] = 0
                    phase[(f >= 3000) & (f < 10e6)] = -90 * isgn[(f >= 3000) & (f < 10e6)]

                elif domain == 'time':
                    # angular frequencies
                    a = 2 * np.pi * 10
                    b = 2 * np.pi * 25
                    c = 2 * np.pi * 400
                    d = 2 * np.pi * 3000

                    # filter parameters
                    fref = np.array([30e3])
                    s = 2j * np.pi * fref
                    Href = (s * (s + c)) / ((s + a) * (s + b) * (s + d))
                    lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                    k = (1.0 / (lim_ref * np.abs(Href))).item()

                    # define numerator and denominator
                    num = [k, k * c, 0]
                    den = [1, (a + b + d), (a * b + a * d + b * d), a * b * d]

            elif quantity == 'Epns':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[(f >= 1) & (f < 3000)] = 1 / (0.8 * np.sqrt(2))
                    weight_fun[(f >= 3000) & (f < 10e6)] = 1 / (2.7e-4 * f[(f >= 3000) & (f < 10e6)] * np.sqrt(2))

                    phase[(f >= 1) & (f < 3000)] = 0
                    phase[(f >= 3000) & (f < 10e6)] = -90 * isgn[(f >= 3000) & (f < 10e6)]

                elif domain == 'time':
                    # angular frequency
                    a = 2 * np.pi * 3000

                    # filter parameters
                    fref = np.array([1])
                    s = 2j * np.pi * fref
                    Href = 1 / (s + a)
                    lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                    k = (1.0 / (lim_ref * np.abs(Href))).item()

                    # define numerator and denominator
                    num = [k]
                    den = [1, a]

        elif receptor == 'public':
            if quantity == 'B':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[(f >= 1) & (f < 8)] = 1 / (4e-2 / f[(f >= 1) & (f < 8)] ** 2 * np.sqrt(2))
                    weight_fun[(f >= 8) & (f < 25)] = 1 / (5e-3 / f[(f >= 8) & (f < 25)] * np.sqrt(2))
                    weight_fun[(f >= 25) & (f < 50)] = 1 / (2e-4 * np.sqrt(2))
                    weight_fun[(f >= 50) & (f < 400)] = 1 / (2e-4 * np.sqrt(2))
                    weight_fun[(f >= 400) & (f < 3000)] = 1 / (8e-2 / f[(f >= 400) & (f < 3000)] * np.sqrt(2))
                    weight_fun[(f >= 3000) & (f < 10e6)] = 1 / (2.7e-5 * np.sqrt(2))

                    phase[(f >= 1) & (f < 8)] = 180 * isgn[(f >= 1) & (f < 8)]
                    phase[(f >= 8) & (f < 25)] = 90 * isgn[(f >= 8) & (f < 25)]
                    phase[(f >= 25) & (f < 50)] = 0
                    phase[(f >= 50) & (f < 400)] = 0
                    phase[(f >= 400) & (f < 3000)] = 90 * isgn[(f >= 400) & (f < 3000)]
                    phase[(f >= 3000) & (f < 10e6)] = 0

                elif domain == 'time':
                    # angular frequency
                    a = 2 * np.pi * 8
                    b = 2 * np.pi * 25
                    c = 2 * np.pi * 300
                    d = 2 * np.pi * 3000

                    # filter parameters
                    fref = np.array([1])
                    s = 2j * np.pi * fref
                    Href = (s ** 2. * (s + c)) / ((s + a) * (s + b) * (s + d))
                    lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                    k = (1.0 / (lim_ref * np.abs(Href))).item()

                    # define numerator and denominator
                    num = np.array([k, k * c, 0, 0])
                    den = np.array([1, (a + b + d), (a * b + a * d + b * d), (a * b * d)])

            elif quantity == 'Ecns':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[(f >= 1) & (f < 10)] = 1 / (0.1 / f[(f >= 1) & (f < 10)] * np.sqrt(2))
                    weight_fun[(f >= 10) & (f < 25)] = 1 / (0.01 * np.sqrt(2))
                    weight_fun[(f >= 25) & (f < 1000)] = 1 / (4e-4 * f[(f >= 25) & (f < 1000)] * np.sqrt(2))
                    weight_fun[(f >= 1000) & (f < 3000)] = 1 / (0.4 * np.sqrt(2))
                    weight_fun[(f >= 3000) & (f < 10e6)] = 1 / (1.35e-4 * f[(f >= 3000) & (f < 10e6)] * np.sqrt(2))

                    phase[(f >= 1) & (f < 10)] = 90 * isgn[(f >= 1) & (f < 10)]
                    phase[(f >= 10) & (f < 25)] = 0
                    phase[(f >= 25) & (f < 1000)] = -90 * isgn[(f >= 25) & (f < 1000)]
                    phase[(f >= 1000) & (f < 3000)] = 0
                    phase[(f >= 3000) & (f < 10e6)] = -90 * isgn[(f >= 3000) & (f < 10e6)]

                elif domain == 'time':
                    # angular frequencies
                    a = 2 * np.pi * 10
                    b = 2 * np.pi * 25
                    c = 2 * np.pi * 1000
                    d = 2 * np.pi * 3000

                    # filter parameters
                    fref = np.array([30e3])
                    s = 2j * np.pi * fref
                    Href = (s * (s + c)) / ((s + a) * (s + b) * (s + d))
                    lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                    k = (1.0 / (lim_ref * np.abs(Href))).item()

                    num = [k, k * c, 0]
                    den = [1, (a + b + d), (a * b + a * d + b * d), a * b * d]

            elif quantity == 'Epns':
                if domain == 'freq':
                    # get sign
                    isgn = np.sign(f)
                    f = np.abs(f)
                    # Initialize output
                    weight_fun = np.zeros(f.shape)
                    phase = np.zeros(f.shape)

                    weight_fun[(f >= 1) & (f < 3000)] = 1 / (0.4 * np.sqrt(2))
                    weight_fun[(f >= 3000) & (f < 10e6)] = 1 / (1.35e-4 * f[(f >= 3000) & (f < 10e6)] * np.sqrt(2))

                    phase[(f >= 1) & (f < 3000)] = 0
                    phase[(f >= 3000) & (f < 10e6)] = -90 * isgn[(f >= 3000) & (f < 10e6)]

                elif domain == 'time':
                    # angular frequency
                    a = 2 * np.pi * 3000

                    # filter parameters
                    fref = np.array([1])
                    s = 2j * np.pi * fref
                    Href = 1 / (s + a)
                    lim_ref = icnirp_limit(fref, year, receptor, quantity) * np.sqrt(2)
                    k = (1.0 / (lim_ref * np.abs(Href))).item()

                    num = [k]
                    den = [1, a]

    # Assign outputs
    if domain == 'time':
        return num, den
    elif domain == 'freq':
        return weight_fun, phase





