import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import signac

from util.decorators import job_chdir
from util.util import block_avg, calc_density


@job_chdir
def calc_density(job):
    top_file = 'sample.gro'
    top_file = 'sample.trr'
    data_file = 'rho.txt'
    img_file = 'rho.pdf'

    if not all([os.path.isfile(top_file) and os.path.isfile(trh_file)]):
        return

    trj = md.load(trj_file, top=top_file)
    rho = calc_density(trj, units='macro')
    data = np.vstack([trj.time, rho])
    np.savetxt(data_file, np.transpose(data),
               header='# Time (ps)\tDensity (kg/m^3)')

    fig, ax = plt.subplots()
    ax.plot(trj.time, rho)
    ax.set_xlabel('Simulation time (ps)')
    ax.set_ylabel('Density (kg/m^3)')
    ax.set_title('Box of C_{}'.format(job.sp()['C_n']))
    fig.savefig(img_file)


if __name__ == '__main__':
    import flow
    flow.run()
