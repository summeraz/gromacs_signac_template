Template for MoSDeF- and GROMACS-centric project managed with signac

See [generic signac template](https://github.com/glotzerlab/signac-project-template) and [a fork](https://github.com/summeraz/monolayer_screening)

# Versions

`signac-flow` underwent significant changes from `0.4` to `0.5`. The master
branch will support versions `0.5` and newer but an older branch (`flow-0.4`)
will be maintained to support version `0.4`.

# Installation

First, verify that each package in `requirements.txt` is installed. If not
already installed, They should all be available via `pip`, `conda`, etc.
The channels `mosdef` and `glotzer` or `conda-forge` are also required.

```
python setup.py install
```

# Get started:

Initialize project:

```
python src/init.py
```

Build systems and perform energy minimization:

```
python src/project.py run initialize
python src/project.py run em  # energy minimization
```

Submit equilibration and production runs to the cluster
(can also be run locally as above)

```
python src/project.py submit -o equil
python src/project.py submit -o sample
```

At any time you can evaluate the status of each job with:

```
python src/project.py status -d
```

Analyze system (here, simple timeseries plots of density)

The status of the project analysis can be examined via:
```
python src/analysis.py status -d
```

To submit the `calc_density` operation (which we will run locally) use:
```
python src/analysis.py run calc_density
```

This operation will write a density timeseries to a text file as well as
a plot in each job's workspace.

To compare the density of all jobs on a single plot, run the
`plot_density.py` script in the `plotting` directory:
```
python plotting/plot_density.py
```

Look at the PDF file located in the `plotting` directory. It should look something like this:

![image](https://user-images.githubusercontent.com/7935382/28077533-a8a43f84-6627-11e7-9370-1206160d185d.png)
