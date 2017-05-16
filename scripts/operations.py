import logging
from math import ceil
import os
import subprocess

from pkg_resources import resource_filename
import numpy as np
import mbuild as mb

from util.decorators import job_chdir
    

logger = logging.getLogger(__name__)


@job_chdir
def initialize(job):
    "Inialize the simulation"

    with job:
        water = mb.load(os.path.join(job._project.root_directory(), 'scripts/util/lib/tip3p.mol2'))
        water.name = 'tip3p'
        n_water = job.statepoint()['n_water']
        #seed = job.statepoint()['seed']
        system_box = mb.Box(3*[(n_water/30)**(1/3)]) # Initialize box to 30 waters/nm^3
        system = mb.fill_box(compound=water,
                             n_compounds=n_water,
                             box=system_box,
                             )#seed=seed)
        system.save('init.gro',
            overwrite=True)
        system.save('init.top',
            forcefield_files = os.path.join(job._project.root_directory(), 'scripts/util/lib/tip3p.xml'),
            overwrite=True)

@job_chdir
def minimize(job):
    "Perform an energy minimization"
    grompp = _grompp_str(job, 'em', 'init', 'init')
    grompp_proc = subprocess.Popen(grompp.split())
    grompp_proc.communicate()
    mdrun = _mdrun_str(job, 'em')
    mdrun_proc = subprocess.Popen(mdrun.split())
    mdrun_proc.communicate()


@job_chdir
def equilibrate(job):
    "Perform a equilibration"
    grompp = _grompp_str(job, 'equil', 'em', 'init')
    grompp_proc = subprocess.Popen(grompp.split())
    grompp_proc.communicate()
    mdrun = _mdrun_str(job, 'equil')
    mdrun_proc = subprocess.Popen(mdrun.split())
    mdrun_proc.communicate()


@job_chdir
def sample(job):
    "Perform a production run"
    grompp = _grompp_str(job, 'sample', 'equil', 'init')
    grompp_proc = subprocess.Popen(grompp.split())
    grompp_proc.communicate()
    mdrun = _mdrun_str(job, 'sample')
    mdrun_proc = subprocess.Popen(mdrun.split())
    mdrun_proc.communicate()

def auto(job):
    "This is a meta-operation to execute multiple operations."
    from gromacs_signac_template.project import get_project
    project = get_project()
    logger.info("Running meta operation 'auto' for job '{}'.".format(job))
    for i in range(10):
        next_op = project.next_operation(job)
        if next_op is None:
            logger.info("No next operation, exiting.")
            break
        else:
            logger.info("Running next operation '{}'...".format(next_op))
            func = globals()[next_op.name]
            func(job)
    else:
        logger.warning("auto: Reached max # operations limit!")

def _grompp_str(job, op_name, gro_name, sys_name):
    """Helper function, returns grompp command string for operation """
    grompp_str = ('gmx grompp -f {0}/scripts/util/mdp_files/{1}.mdp -c {2}/{3}.gro '
        '-p {2}/{4}.top -o {2}/{1}.tpr'
        ''.format(
            job._project.root_directory(),
            op_name,
            job.workspace(), 
            gro_name,
            sys_name))
    return grompp_str

def _mdrun_str(job, op_name):
    """Helper function, returns mdrun command string for operation """
    mdrun_str = 'gmx mdrun -v -deffnm {} -nt 12'.format(op_name)
    return mdrun_str
