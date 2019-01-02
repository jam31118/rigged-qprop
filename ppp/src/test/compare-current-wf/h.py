import numpy as np
wf_file_new, wf_file_qprop = "../new-prop-calc-home/current-wf.bin", "../qprop-calc-home/current-wf.bin"
wf_arr_new, wf_arr_qprop = (np.fromfile(file_path, dtype=complex) for file_path in [wf_file_new, wf_file_qprop])

are_close_arr = np.isclose(wf_arr_new, wf_arr_qprop)
all_close = np.all(are_close_arr)

print(wf_arr_new)
print(wf_arr_qprop)

print(are_close_arr)
print("all close: {}".format(all_close))

