Template for MoSDeF- and GROMACS-centric project managed with signac

See [generic signac template](https://github.com/glotzerlab/signac-project-template) and [a fork](https://github.com/summeraz/monolayer_screening)

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

Sample analysis coming soon ...
