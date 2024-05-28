from qxREAD import qxmd_eigs
myeigs=qxmd_eigs()
myeigs.read_eigs('../../data')
print((myeigs.occs[0]))
print((myeigs.eigs[0]))
frame=len(myeigs.occs[0])-1 #get last frame, python 0 indexed so -1
sigma=0.1 #eV smearing for DOS
NEDOS=500 #Number of DOS pts
EMIN=-10  #min energy for DOS
EMAX=10   #max energy for DOS
dos=myeigs.compute_dos(frame,sigma,NEDOS,EMIN,EMAX)
print(dos)
