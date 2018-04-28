import subprocess as sp
from qprop.parameter import update_param_in_file

def get_energy(effpot_alpha):
    update_param_in_file('initial.param','effpot-alpha','double',effpot_alpha)
    command_words = ['./hydrogen_im']
    run_result = sp.run(command_words, check=True, stdout=sp.PIPE, stderr=sp.PIPE)
    energy = read_scalar_from_file('initial-energy.dat')
    return energy


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

for effpot_alpha in [0.01 * i for i in range(4)]:
    energy = get_energy(effpot_alpha)
    print("effpot-alpha: {0} / energy: {1}".format(effpot_alpha, energy))

