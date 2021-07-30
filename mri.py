import numpy as np


def mri_limit(esd, organ, quantity, mode='NOM'):
    """
    MRI_LIMIT returns the threshold (perception level) for MRI application according to IEC 60601-2-33:2010

    :param esd: effective stimulus duration
    :param organ: 'heart' or 'PNS'
    :param quantity: 'E', 'dBdt'
    :param mode: 'NOM' or 'FLOM' (default is NOM)
    :return:
        * limit: threshold (perception level) for the given input

    EXAMPLE: (to be checked)
    >>> esd = 5e-3 #s
    >>> organ = 'heart'
    >>> quantity = 'E'
    >>> limit = mrilimit(esd, organ, quantity)
    >>> print('{:.2f} (V/m)'.format(limit))
    2.47 (V/m)
    """

    # convert to numpy.array if esd is a scalar value
    if isinstance(esd, (int, float)):
        esd = np.array([esd])

    if organ.lower() == 'heart':
        if quantity.lower() ==  'e':
            limit = 2. / (1 - np.exp(-esd / 3e-3))
        elif quantity.lower() == 'dbdt':
            limit = 20 / (1 - np.exp(-esd / 3e-3))
        # elif quantity.lower() ==  'b':
        #     dBdt_max = 20 / (1 - np.exp(-esd / 3e-3))
        #     limit = dBdt_max / (2 * np.pi * f) / np.sqrt(2) # be aware! It is an RMS value
    elif organ.lower() == 'pns':
        # mode selection
        if mode.lower() == 'nom':
            RF = 0.8
        elif  mode.lower() == 'flom':
            RF = 1

        if quantity.lower() == 'e':
            limit = RF * 2.2 * (1 + 0.36e-3 / esd)
        elif quantity.lower() == 'dbdt':
            limit = RF * 20 * (1 + 0.36e-3 / esd)
        # elif quantity.lower() == 'b':
        #     dBdt_max = RF * 20 * (1 + 0.36e-3 / esd)
        #     limit = dBdt_max / (2 * np.pi * f) / np.sqrt(2) # be aware! It is an RMS value

    # convert to scalar if the numpy.array has only one element
    if limit.size == 1:
        limit = np.asscalar(limit)

    return limit