import numpy as np
import mdtraj as md


def calc_density(traj, units='macro'):
    """
    This function is taken directly from mtools (see
    github.com/mattwthompson/mtools).
    """
    vol = np.product(traj.unitcell_lengths, axis=1)
    total_mass = np.sum([x.element.mass for x in traj.topology.atoms])
    if units == 'nano':
        rho = total_mass/vol
    elif units == 'macro':
        rho = total_mass/vol * 1.66054  # Convert from amu/nm^3 to kg/m^3
    else:
        raise ValueError('Unsupported units!')
    return rho


def block_avg(data, block_size):
    """
    Break a 2d numpy array into blocks

    This function is taken directly from Tim Moore's `block_avg` (see
    https://github.com/tcmoore3/block_avg).

    Parameters
    ----------
    data : np.ndarray, shape=(m, n)
        The data to block; must be a 2-dimensional array
    block_size : int
        The size of each block
    Returns
    -------
    blocks : np.ndarray, shape=(m/block_size, n)
        The block averaged data
    stds : np.ndarray, shape=(m/block_size, n)
        The standard deviation of each block
    Notes
    -----
    `m` must not necessarily be divisible by `block_size` ; in the case
    that it isn't, the data is trimmed *from the beginning* so that it is.
    """

    remainder = data.shape[0] % block_size
    if remainder != 0:
        data = data[remainder:]
    n_blocks = int(data.shape[0] / block_size)
    data = data.reshape((n_blocks, block_size, -1))
    blocks = np.mean(data, axis=1)
    stds = np.std(data, axis=1)
    return blocks, stds
