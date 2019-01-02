from qprop.core import Qprop20
dir_paths = ['../new-prop-calc-home','../qprop-calc-home']
qs = [Qprop20(dir_path) for dir_path in dir_paths]
msps = [q.momentum_spectrum_polar for q in qs]
for msp, q in zip(msps, qs): msp.readPolarSpectrumData(q)
kProbGrids = [msp.kProbGrid for msp in msps]
import numpy as np
close_arr = np.isclose(*kProbGrids)
print("kProbGrid close array: {}".format(close_arr))
all_close = np.all(close_arr)
print("kProbGrid are close enough: {}".format(all_close))

