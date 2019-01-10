import numpy as np
psi_new, psi_qprop = "./ati-ppp/tsurffpsi.raw", "./ati/tsurffpsi.raw"
psi_arr_new, psi_arr_qprop = (np.fromfile(file_path, dtype=complex) for file_path in [psi_new, psi_qprop])

are_close_arr = np.isclose(psi_arr_new, psi_arr_qprop)
all_close = np.all(are_close_arr)

print("comparing `{}` and `{}`".format(psi_new, psi_qprop))

print("new: ", psi_arr_new)
print("original: ", psi_arr_qprop)
print("diff: ", psi_arr_new - psi_arr_qprop)

print(are_close_arr)
print("all close: {}".format(all_close))
print("`{}` out of `{}` elements are close enough".format(are_close_arr.sum(), are_close_arr.size))

print("indice of not close elements: {}".format(np.arange(are_close_arr.size)[~are_close_arr]))

