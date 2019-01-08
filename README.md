# Rigged QPROP
Tool set dedicated for users of [QPROP][qprop-website]: 
A Schrödinger-solver for intense laser-atom interaction


# Installation
Refer to [wiki][wiki-installation] for detailed installation instructions.


# Citation
> "Qprop is for non-profit use only. 
> If results obtained using Qprop are published, 
> an acknowledgment or citation is compulsory."

An example of citation would be like the following:

D. Bauer and P. Koval, Comput. Phys. Commun. 174, 396 (2006).

V. Mosert and D. Bauer, Comput. Phys. Commun. 207, 452 (2016).


# Visualization
For visualizing the calculation result by `rigged-qprop`, have a look at [visual-qprop](https://github.com/jam31118/qprop), a toolset dedicated to visualize and analyze the calculation results such as animated time evolution, plotting photoelectric momentum spectra etc.


# PPP : Paralleled Post Propagation
Implemented for parallel computing, PPP routine boosts the post propagation, which is often essential for evaluating low-energy structure in momentum spectra through t-SURFF routine.

After compilation along with this PPP support, check out the example using this routine, for example:
``` bash
cd $QPROP_HOME/src/example/ppp/attoclock-ppp
bash ../run.sh
```
Noting that, for the `attoclock` example, more than 90% of the timesteps for propagation is the field-free post-propagation, it can benefit from `PPP` quite a lot.


[qprop-website]: http://qprop.de
[wiki-installation]: https://github.com/jam31118/rigged-qprop/wiki/Installation
