from qxREAD import qxmd_energy_NVE
myeng=qxmd_energy_NVE('../../data/')
print(myeng.steps)
print(myeng.energy)
print(myeng.PE)
print(myeng.KE)
print(myeng.T)
