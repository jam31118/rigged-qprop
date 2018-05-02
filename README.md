# Rigged QPROP
Tool set dedicated for users of QPROP: A Schrödinger-solver for intense laser-atom interaction


# Citation
"Qprop is for non-profit use only. If results obtained using Qprop are published, an acknowledgment or citation is compulsory."

An example of citation would be like the following:

D. Bauer and P. Koval, Comput. Phys. Commun. 174, 396 (2006).

V. Mosert and D. Bauer, Comput. Phys. Commun. 207, 452 (2016).


# Prerequisites
## Environment Variable
`export QPROP_HOME=/path/to/home/directory/of/cloned/repository`

`export QPROP_DEP_HOME=/path/in/which/dependencies/are/intalled`

The directory tree pointed by `$QPROP_DEP_HOME` should be like the following:

```bash
$QPROP_DEP_HOME/
  boost_1_65_1/
  gsl-2.4/
  openmpi-1.10.7/
```

where subdirectories are intallation directory of each dependency.

For example, the directory tree of `gsl-2.4/` for GSL library, should be:

```bash
$QPROP_DEP_HOME/
  gsl-2.4/
    bin/
    include/
    lib/
    share/
```

