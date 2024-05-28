from qxREAD import qxmd_eigs
myeigs=qxmd_eigs()
myeigs.read_td_eigs('../../data')
print((myeigs.td_occs[3]))
print((myeigs.td_eigs[3]))
print((myeigs.td_occs[4]))
print((myeigs.td_eigs[4]))
