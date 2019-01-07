import numpy as np


type_name = 'float'
file_name = "data-{}.bin".format(type_name)
N = 13


a = np.arange(N, dtype=type_name) * 0.1

with open(file_name, "wb") as f:
    a.tofile(f)

