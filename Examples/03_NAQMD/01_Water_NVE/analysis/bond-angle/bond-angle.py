import numpy as np
from qxREAD import qxmd_box,qxmd_ions
mybox=qxmd_box('../../data/')
mytraj=qxmd_ions('../../data/')
mytraj.to_real(mybox)
print(mytraj.ions[0])
OH1=np.zeros(mytraj.nframes)
OH2=np.zeros(mytraj.nframes)
HOHangle=np.zeros(mytraj.nframes)
for i in range(mytraj.nframes) :
	frame=mytraj.ions[i].copy() #mytraj.ions is a list, if we don't use copy frame will act as pointer
        #			     and we may accidently overwrite imported data
	drij=frame[0,:]-frame[1,:]
	drik=frame[0,:]-frame[2,:]
	rij=np.sqrt(np.sum(drij*drij))
	rik=np.sqrt(np.sum(drik*drik))
	OH1[i]=rij
	OH2[i]=rik
	drij=drij/rij
	drik=drik/rik
	cosine=np.dot(drij,drik)
	HOHangle[i]=np.arccos(cosine)*180.0/np.pi

print(OH1)
print(OH2)
print(HOHangle)

