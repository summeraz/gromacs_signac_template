Template for MoSDeF- and GROMACS-centric project managed with signac

See [generic signac template](https://github.com/glotzerlab/signac-project-template) and [a fork](https://github.com/summeraz/monolayer_screening)

# Installation

Install these packages:

* [block_avg](https://github.com/tcmoore3/block_avg)
* [mtools](https://github.com/mattwthompson/mtools)

Install other dependencies:

```
python setup.py install
```


# Get started:

Initialize project:

```
python gromacs_signac_template/init.py
```

Build systems and perform energy minimization:

```
python scripts/run.py initialize
python scripts/run.py minimize
```

Submit equilibration and production runs to the cluster
(can also be done locally as above)

```
python submit.py -j equilibrate
python submit.py -j sample
```

At any time you can evaluate the status of each job with:

```
python gromacs_signac_template/status.py -d
```

Analyze system (here, simple timeseries plots of density)

```
python analysis/calc_density.py
```

Look at the PDF file located in the ```workspace``` directory.
