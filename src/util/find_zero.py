import subprocess as sp
from qprop.parameter import update_param_in_file

from os.path import isfile
from sys import argv, stderr

if len(argv) != 3:
    print("[ERROR] Enter (1) imaginary propagation binary " 
            + "/ (2) first ioniation potential in atomic unit", file=stderr)
    exit(1)
# imag_prop_bin = './hydrogen_im'
imag_prop_bin = argv[1]
assert isfile(imag_prop_bin)
Ip = float(argv[2])


#exit(1)

def get_energy(effpot_alpha, Ip, imag_prop_bin, energy_filepath='initial-energy.dat'):
    assert isfile(imag_prop_bin)
    if hasattr(effpot_alpha, '__getitem__') and hasattr(effpot_alpha, '__len__'):
#    if type(effpot_alpha) in [list, tuple]:
        assert len(effpot_alpha) == 1
        effpot_alpha = effpot_alpha[0]
    if effpot_alpha < 0: effpot_alpha = 0.0
    update_param_in_file('initial.param','effpot-alpha','double',effpot_alpha)
    command_words = [imag_prop_bin]
    run_result = sp.run(command_words, check=True, stdout=sp.PIPE, stderr=sp.PIPE)
    run_result.check_returncode()
    assert isfile(energy_filepath)
    energy = read_scalar_from_file(energy_filepath)
    diff = energy - Ip
    print("effpot_alpha: {0} / diff: {1} / energy: {2}".format(effpot_alpha, diff, energy), flush=True)
    return diff


from os.path import isfile
def read_scalar_from_file(filepath):
    assert isfile(filepath)
    scalar = None
    content = None
    with open(filepath, 'r') as f: content = f.read()
    not_empty_lines = [line for line in content.split('\n') if line != '']
    if len(not_empty_lines) == 1:
        line = not_empty_lines[0]
        try: scalar = float(line)
        except: raise Exception("Failed to convert '{0}' into float".format(line))
    else: 
        err_mesg = "The number of non-empty line in file {0} should be one. However, it was {1}"
        raise IOError(err_mesg.format(filepath, len(not_empty_lines)))
    return scalar

# Ip = -4.4576e-01
initial_guess = 3.0
print("Starting with initial guess {0:.4}".format(initial_guess), flush=True)

fargs = (Ip, imag_prop_bin)
from scipy.optimize import fsolve
alpha_at_zero = fsolve(get_energy, initial_guess, args=fargs)
print("alpha_at_zero: {0}".format(alpha_at_zero), flush=True)


## generate global values
#import numpy as np
#alpha_array = np.linspace(0,4,200)
#energy_array = np.empty_like(alpha_array)
#for idx, alpha in enumerate(alpha_array):
#    energy_array[idx] = get_energy(alpha, *fargs)
#np.savetxt('alpha-array.txt',alpha_array)
#np.savetxt('energy-array.txt',energy_array)


## [TODO] From name of species e.g. 'xe' to set of parameters such as nuclear-charge, etc.
## then find zero (but avoiding .. a strange one. - fsolve doesn't seems to work for some initial guess) 

## Plot along different, equidistanced and along long wide range, energy

#for effpot_alpha in [0.01 * i for i in range(4)]:
#    energy = get_energy(effpot_alpha, Ip, imag_prop_bin)
#    print("effpot-alpha: {0} / energy: {1}".format(effpot_alpha, energy))

# alpha double 2.6442958 (previous, 0.05 case)

## for pot-cutoff = 40.0
# 2.40694195 (delta-r = 0.1)
# 2.64429747 (delta-r = 0.05)
# 2.567087   (delta-r = 0.025)

## for pot-cutoff = 100.0
# 1.85427089 (delta-r 0.2)
# 2.40694195 (delta-r 0.1)
# 2.64429747 (delta-r 0.05)
# 2.567087   (delta-r 0.025)

