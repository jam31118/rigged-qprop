import numpy as np
psi_new, psi_qprop = "../new-prop-calc-home/tsurffpsi.raw", "../qprop-calc-home/tsurffpsi.raw"
psi_arr_new, psi_arr_qprop = (np.fromfile(file_path, dtype=complex) for file_path in [psi_new, psi_qprop])

are_close_arr = np.isclose(psi_arr_new, psi_arr_qprop)
all_close = np.all(are_close_arr)

print(psi_new)
print(psi_qprop)

print(are_close_arr)
print("all close: {}".format(all_close))

