#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories."""
import logging
import argparse
from hashlib import sha1

import signac
import numpy as np


def main(args):
    project = signac.init_project('Alkanes')
    statepoints_init = []
    for C_n in [6, 8, 10]:
        statepoint = dict(
                    # length of alkane
                    C_n = C_n)
        project.open_job(statepoint).init()
        statepoints_init.append(statepoint)

    # Writing statpoints to hash table as backup
    project.write_statepoints(statepoints_init)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Initialize the data space.")
    parser.add_argument(
        '-n', '--num-replicas',
        type=int,
        default=1,
        help="Initialize multiple replications.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(args)
