Template for MoSDeF- and GROMACS-centric project managed with signac

See [generic signac template](https://github.com/glotzerlab/signac-project-template) and [a fork](https://github.com/summeraz/monolayer_screening)

# Get started:

Initialize project:

```
python gromacs_signac_template/init.py 654
```

Build systems and perform energy minimization:

```
python scripts/run.py initialize
python scripts/run.py minimize
```

Submit equilibration and production runs to the cluster

```
python gromacs_signac_template/submit.py -j equilibrate
python gromacs_signac_template/submit.py -j sample
```

Sample analysis coming soon ...
