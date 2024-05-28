import numpy as np
'''
--------------------------------------------------------
This is module for Importing data from QXMD simulations
May add VASP support in the Future 


Current tools
	qxmd_eigs -> Object that stores QXMD eigen values
	qxmd_phonon_trans_prob -> Stores tranistion proability from state i->j 
				  for NAQMD simmulation
	qxmd_ions -> Object that stores QXMD atoms along trajectory
	qxmd_box  -> Object that stores QXMD box
	qxmd_energy_OPT --> get Energy for optimization run
	qxmd_energy_NVE --> get Energy for NVE run
	qxmd_energy_NVT --> get Energy for NVT run

--------------------------------------------------------
'''
class qxmd_eigs :
	def __init__(self):
		self.eigs=[]
		self.occs=[]
		self.td_eigs=[]
		self.td_occs=[]
	def read_eigs(self,data_dir) :
		filename=data_dir+'/qm_eig.d'
		eig_file=open(filename,'r')
		print('File opened ' + filename)
		filename=data_dir+'/qm_fer.d'
		fer_file=open(filename,'r')
		print('File opened ' + filename)

		line=eig_file.readline()
		line=eig_file.readline()
		line=line.strip().split()
		step=int(line[0])
		scf=int(line[1])
		nbands=int(line[2])
		self.nbands=nbands
		print('nbands : '+str(nbands))

		line=fer_file.readline()
		line=fer_file.readline()
		line=line.strip().split()
		step=int(line[0])
		scf=int(line[1])
		EFER=float(line[2])
		print('step : ' +str(step))
		print('Fermi Energy (Ry) : ' +str(EFER))
		print('Fermi Energy (eV) : ' +str(EFER*13.6))
		
		for i in range(nbands) :
			self.eigs.append([])
			self.occs.append([])
		for i in range(nbands) :
			line=eig_file.readline()
			line=line.strip().split()
			#print(line)
			idx=int(line[0])
			E=float(line[1])
			OCC=float(line[2])
			self.eigs[i].append((E-EFER)*13.6)
			self.occs[i].append(OCC)
		while True :
		
			line=fer_file.readline()
			if not line :
				break
			line=line.strip().split()
			step=int(line[0])
			scf=int(line[1])
			EFER=float(line[2])
			print('step : ' +str(step))
			print('Fermi Energy (Ry) : ' +str(EFER))
			print('Fermi Energy (eV) : ' +str(EFER*13.6))
			line=eig_file.readline()
			line=line.strip().split()
			step=int(line[0])
			scf=int(line[1])
			nbands=int(line[2])
			for i in range(nbands) :
				line=eig_file.readline()
				line=line.strip().split()
				idx=int(line[0])
				E=float(line[1])
				OCC=float(line[2])
				self.eigs[i].append((E-EFER)*13.6)
				self.occs[i].append(OCC)
			
	
		eig_file.close()
		fer_file.close()
			
	def read_td_eigs(self,data_dir) :
		filename=data_dir+'/qm_td_eig.d'
		eig_file=open(filename,'r')
		print('File opened ' + filename)
		filename=data_dir+'/qm_fer.d'
		fer_file=open(filename,'r')
		print('File opened ' + filename)

		line=eig_file.readline()
		line=eig_file.readline()
		line=line.strip().split()
		step=int(line[0])
		scf=int(line[1])
		nbands=int(line[2])
		self.nbands=nbands
		print('nbands : '+str(nbands))

		line=fer_file.readline()
		line=fer_file.readline()
		line=line.strip().split()
		step=int(line[0])
		scf=int(line[1])
		EFER=float(line[2])
		print('step : ' +str(step))
		print('Fermi Energy (Ry) : ' +str(EFER))
		print('Fermi Energy (eV) : ' +str(EFER*13.6))
		
		for i in range(nbands) :
			self.td_eigs.append([])
			self.td_occs.append([])
		for i in range(nbands) :
			line=eig_file.readline()
			line=line.strip().split()
			#print(line)
			idx=int(line[0])
			E=float(line[1])
			OCC=float(line[2])
			self.td_eigs[i].append((E-EFER)*13.6)
			self.td_occs[i].append(OCC)
		while True :
		
			line=fer_file.readline()
			if not line :
				break
			line=line.strip().split()
			step=int(line[0])
			scf=int(line[1])
			EFER=float(line[2])
			print('step : ' +str(step))
			print('Fermi Energy (Ry) : ' +str(EFER))
			print('Fermi Energy (eV) : ' +str(EFER*13.6))
			line=eig_file.readline()
			line=line.strip().split()
			step=int(line[0])
			scf=int(line[1])
			nbands=int(line[2])
			for i in range(nbands) :
				line=eig_file.readline()
				line=line.strip().split()
				idx=int(line[0])
				E=float(line[1])
				OCC=float(line[2])
				self.td_eigs[i].append((E-EFER)*13.6)
				self.td_occs[i].append(OCC)
			
	
		eig_file.close()
		fer_file.close()
	def compute_dos(self,frame,sigma,NEDOS,EMIN,EMAX) :
		deltaE=(EMAX-EMIN)/NEDOS
		myE=EMIN
		dos=np.zeros((NEDOS,2))
		prefac=1/(sigma*np.sqrt(2*np.pi))
		for i in range(NEDOS) :
			for j in range(self.nbands) :
				eig=self.eigs[j][frame]
				dos[i,0]=myE
				exponent=-0.5*(eig-myE)**2/(sigma**2)
				if (exponent<-30) :
					dos[i,1]=dos[i,1]+0.0
				else :
					dos[i,1]=dos[i,1]+ prefac*np.exp(exponent)
			myE=myE+deltaE
		return(dos)
class qxmd_phonon_trans_prob :
	def __init__(self,directory,i,j):
		self.step=[]
		self.prob=[]
		filename=directory+'/qm_fsshprob_'+str(i)+'to'+str(j)+'-u.d' 
		read_file=open(filename,'r')
		line = read_file.readline() # Comment
		while True :
			line = read_file.readline() 
			if not line :
				break
			line=line.strip().split()
			self.step.append(int(line[0]))
			self.prob.append(float(line[2]))
			
		
class qxmd_box :
	def __init__(self,directory):
		self.La=None
		self.Lb=None
		self.Lc=None
		self.alpha=None
		self.beta=None
		self.gamma=None
		self.Hmat=np.zeros((3,3))
		H=np.zeros((3,3))
		bohrtoang=0.529177
		filename=directory+'/qm_box.d'
		read_file=open(filename,'r')
		print('File Opened : '+filename)
		line = read_file.readline() # Comment
		line = read_file.readline() # Header
		line = read_file.readline()
		line = line.strip().split()
		self.steps=int(line[0])
		la=float(line[1])*bohrtoang
		lb=float(line[2])*bohrtoang
		lc=float(line[3])*bohrtoang
		angle1=float(line[4])
		angle2=float(line[5])
		angle3=float(line[6])
		lal=angle1*np.pi/180.00
		lbe=angle2*np.pi/180.00
		lga=angle3*np.pi/180.00
		#Box Conversion stolen from RXMD
		hh1=lc*(np.cos(lal)-np.cos(lbe)*np.cos(lga))/np.sin(lga)
		hh2=lc*np.sqrt(1.0-np.cos(lal)**2-np.cos(lbe)**2-np.cos(lga)**2 + 2*np.cos(lal)*np.cos(lbe)*np.cos(lga) )/np.sin(lga)
		H[0,0]=la;             H[1,0]=0.00;        H[2,0]=0.00
		H[0,1]=lb*np.cos(lga); H[1,1]=lb*np.sin(lga); H[2,1]=0.00
		H[0,2]=lc*np.cos(lbe); H[1,2]=hh1;         H[2,2]=hh2


		self.La=la
		self.Lb=lb
		self.Lc=lc
		self.alpha=angle1
		self.beta=angle2
		self.gamma=angle3
		self.Hmat=H
class qxmd_ions :
	def __init__(self,directory):
		filename=directory+'/qm_ion.d'
		read_file=open(filename,'r')
		print('File Opened : '+filename)
		line = read_file.readline()
		vals=list()
		nstp=0
		while True:
			line = read_file.readline()
			if not line:
				break
				#end if
			line = line.strip().split()
			#print(line)
			#print(nstp)
			nstp=nstp+1
			nspc=int(line[1])
			spc=list()
			for i in range(nspc):
				spc.append(int(line[2+i]))
				# end for
			#print(spc)
			self.natom=sum(spc)

			line = read_file.readline()
			line = line.strip().split()
			scale=float(line[0])
			new=list()
			for i in range(int(self.natom)//int(3)):
				line = read_file.readline()
				line = line.strip().split()
				new.extend(line)
				#end for
			count=1
			count2=0
			mat=np.zeros(shape=(self.natom,3))
			for i in new:
				if ( count % 3 ==1):
					mat[count2,0]=float(i)* scale
				# end if
				if ( count  % 3== 2):
					mat[count2,1]=float(i) *scale
				#end if
				if (count % 3 ==0):
					mat[count2,2]=float(i) * scale
					count2=count2+1
				#end if
				count=count+1
				#end for
			vals.append(mat)
			#end while
		read_file.close()
		self.ions=vals
		self.spcs=spc
		self.nframes=len(vals)
		
	def to_real(self,mybox):
		for i in range(self.nframes) :
			frame=self.ions[i]
			end=self.natom
			#print(self.natom)
			frame[0:end,0]=(mybox.Hmat[0,0]*frame[0:end,0]+mybox.Hmat[0,1]*frame[0:end,1]+
				mybox.Hmat[0,2]*frame[0:end,1] )
			
			frame[0:end,1]=(mybox.Hmat[1,0]*frame[0:end,0]+mybox.Hmat[1,1]*frame[0:end,1]+
				mybox.Hmat[1,2]*frame[0:end,2] )

			frame[0:end,2]=(mybox.Hmat[2,0]*frame[0:end,0]+mybox.Hmat[2,1]*frame[0:end,1]+
				mybox.Hmat[2,2]*frame[0:end,2] )


class qxmd_energy_OPT :
	def __init__(self,directory):
		self.energy=[]
		self.steps=[]
		filename=directory+'/md_eng.d'
		read_file=open(filename,'r')
		print('File Opened : '+filename)
		line = read_file.readline()
		while True:
			line = read_file.readline()
			if not line:
				break
			line = line.strip().split()
			self.steps.append(int(line[0]))
			self.energy.append(27.2114*float(line[1]))
		read_file.close()
		
class qxmd_energy_NVE :
	def __init__(self,directory):
		self.energy=[]
		self.PE=[]
		self.KE=[]
		self.T=[]
		self.steps=[]
		filename=directory+'/md_eng.d'
		read_file=open(filename,'r')
		print('File Opened : '+filename)
		line = read_file.readline()
		while True:
			line = read_file.readline()
			if not line:
				break
			line = line.strip().split()
			self.steps.append(int(line[0]))
			self.energy.append(27.2114*float(line[1]))
			self.PE.append(27.2114*float(line[2]))
			self.KE.append(27.2114*float(line[3]))
			self.T.append(float(line[4]))
		read_file.close()
		
class qxmd_energy_NVT :
	def __init__(self,directory):
		self.energy=[]
		self.PE=[]
		self.KE=[]
		self.T=[]
		self.steps=[]
		self.energyplusBath=[]
		self.BathKE=[]
		self.BathPE=[]
		filename=directory+'/md_eng.d'
		read_file=open(filename,'r')
		print('File Opened : '+filename)
		line = read_file.readline()
		while True:
			line = read_file.readline()
			if not line:
				break
			line = line.strip().split()
			self.steps.append(int(line[0]))
			self.energy.append(27.2114*float(line[1]))
			self.PE.append(27.2114*float(line[2]))
			self.KE.append(27.2114*float(line[3]))
			self.T.append(float(line[4]))
			self.energyplusBath.append(27.2114*float(line[5]))
			self.BathKE.append(27.2114*float(line[6]))
			self.BathPE.append(27.2114*float(line[7]))
		read_file.close()
		
