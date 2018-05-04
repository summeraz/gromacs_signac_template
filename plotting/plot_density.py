import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import signac

fig, ax = plt.subplots()

project = signac.get_project()
for job in project.find_jobs():
    data_file = os.path.join(job.ws(), 'rho.txt')
    if os.path.isfile(data_file):
        data = np.loadtxt(data_file)
        t = data[:, 0]
        rho = data[:, 1]

        t_b, t_std = block_avg(t, 50)
        rho_b, rho_std = block_avg(rho, 50)
        t_b = t_b.reshape(-1)
        t_std = t_std.reshape(-1)
        rho_b = rho_b.reshape(-1)
        rho_std = rho_std.reshape(-1)
        myline, = ax.plot(t_b, rho_b, marker='o', markersize=5,
                          label='C_{}'.format(job.sp()['C_n']))
        ax.fill_between(t_b, rho_b-rho_std, rho_b+rho_std, alpha=0.2,
                        facecolor=myline.get_color())

ax.set_xlabel('Simulation time (ps)')
ax.set_ylabel('Density (kg/m^3)')
ax.legend(title='Alkane length', loc='lower right', fontsize='x-small')
fig.savefig(os.path.join(project.workspace(), 'rho-summary.pdf'))
