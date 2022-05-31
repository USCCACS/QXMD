import numpy as np
from qxREAD import qxmd_phonon_trans_prob
LUMOtoHOMO=qxmd_phonon_trans_prob('../../data',5,4)
#print(LUMOtoHOMO.step)
#print(LUMOtoHOMO.prob)
time_array=np.array(LUMOtoHOMO.step)*0.2412*1e-15 #convert to s
prob_array=np.array(LUMOtoHOMO.prob)
fit=np.polyfit(time_array,prob_array,1)
rate=fit[0]
print("%20s %10.3E" % ('Recombination Rate (s^-1) :' ,rate))
