#################################################################
#                                                               #
#  Makefile for electronic structure calculation by PW method   #
#                                                               #
#                                                               #
#    serial           : make qxmd                               #
#    MPI              : make qxmdmpi                            #
#                                                               #
#################################################################
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES:	.F90 .o .l

CPPDEF    = -DPCHOLESKY=0 -DCCREAL4=1 -DKGBOXR=0
FFTWLIB   = -lrfftw -lfftw -lm
FFTW3LIB  = -lfftw3 -lm
LDFLAG    = -L$(HOME)/FFTW/lib -L$(HOME)/FFTW3/lib
EXTF90    = f90
########################################################FOR DEC ALPHA
#DEC#     LINKER        = f90
#DEC#     LDFLAGS       = -fast $(LDFLAG)
#DEC#     FFLAGS        = -c -fast
#DEC#     MPILIBS       = -L/usr/local/mpi/lib/alpha/ch_p4 -lmpi
#DEC#     PVMLIBS       = -L/usr/local/pvm3/lib/ALPHA -lfpvm3 -lgpvm3 -lpvm3
#DEC#     CPPDEFS       = $(CPPDEF) -DLIBFFTW -DPOINTER64 -DDEC_ALPHA
#DEC#     FFTLIB        = $(FFTWLIB)
#######################################################FOR LINUX (GNU)
##GNU#     LINKER        = mpif90
##GNU#     LDFLAGS       = -O3 -malign-double -funroll-all-loops \
##GNU#                    -fomit-frame-pointer -ffast-math -fstrength-reduce \
##GNU#                    -fthread-jumps -m486 $(LDFLAG)
##GNU#     FFLAGS        = -c -O3 -malign-double -funroll-all-loops \
##GNU#                    -fomit-frame-pointer -ffast-math -fstrength-reduce \
##GNU#                    -fthread-jumps -m486
##GNU#     CPPDEFS       = $(CPPDEF) -DLIBFFTW -DGNU
##GNU#     FFTLIB        = $(FFTWLIB)
#######################################################FOR LINUX (GNU)
#GNU#     LINKER        = gfortran
#GNU#     LDFLAGS       = 
#GNU#     FFLAGS        = -c
#GNU#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DGNU
#GNU#     FFTLIB        = 
#GNU#     EXTF90        = f
##########################################FOR PC LINUX (Fujitsu Fortran&C)
## -X9 : compile sources as Fortran95 source
## -Am : create *.mod (module)
## -Eg : check procedure diffinition and reference, may take long execution time
#FFC#     LINKER        = mpif90
##FFC#     LDFLAGS       = -Kfast -static-flib $(LDFLAG)
##FFC#     FFLAGS        = -c -Kfast -Am -X9 -f 2008 -f 2004 -f 2006
#FFC#     LDFLAGS       = -static-flib $(LDFLAG)
#FFC#     FFLAGS        = -c -Am -X9 -f 2008 -f 2004 -f 2006
#FFC#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DFFC -DVECTOR
#FFC#     FFTLIB        = $(FFTWLIB)
##########################################FOR PC LINUX (Fujitsu Fortran&C) on VT-64
## -X9 : compile sources as Fortran95 source
## -Am : create *.mod (module)
## -Eg : check procedure diffinition and reference, may take long execution time
#FFCV64#     LINKER        = mpif90
##FFCV64#     LDFLAGS       = -Kfast $(LDFLAG)
##FFCV64#     FFLAGS        = -c -Kfast -Am -X9 -f 2008 -f 2004 -f 2006 -f 8340
#FFCV64#     LDFLAGS       = $(LDFLAG)
#FFCV64#     FFLAGS        = -c -Am -X9 -f 2008 -f 2004 -f 2006 -f 8340
#FFCV64#     CPPDEFS       = -traditional $(CPPDEF) -DSSL2VP -DFFC -DPOINTER64 -DVECTOR
##FFCV64#     FFTLIB        = $(FFTWLIB)
#FFCV64#     FFTLIB        = -SSL2 -static
##########################################FOR PC LINUX (Intel Fortran)
#IFC#     LINKER        = mpif90
#IFC##     LDFLAGS       = -axW -static $(LDFLAG)
#IFC##     FFLAGS        = -c -axW -w90 -w95
#IFC##     LDFLAGS       = -static $(LDFLAG)
#IFC#     LDFLAGS       = $(LDFLAG)
#IFC#     FFLAGS        = -c -w90 -w95
#IFC#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DIFC -DVECTOR
#IFC#     FFTLIB        = $(FFTWLIB)
##########################################FOR PC LINUX (Intel Fortran) on x86-64
#IFCX86#     LINKER        = mpif90
#IFCX86##     LDFLAGS       = -axW -static $(LDFLAG)
#IFCX86##     FFLAGS        = -c -axW -w90 -w95
#IFCX86##     LDFLAGS       = -static $(LDFLAG)
#IFCX86#     LDFLAGS       = $(LDFLAG)
#IFCX86#     FFLAGS        = -c -w90 -w95
#IFCX86#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DPOINTER64 -DIFC -DVECTOR
#IFCX86#     FFTLIB        = $(FFTWLIB)
##########################################FOR PC LINUX (Intel Fortran) on VT-64
#IFCV64#     LINKER        = mpif90
#IFCV64##     LDFLAGS       = -axW -static $(LDFLAG)
#IFCV64##     FFLAGS        = -c -axW -warn none
#IFCV64##     LDFLAGS       = -static $(LDFLAG)
#IFCV64##     LDFLAGS       = -L/usr/local/fftw-2.1.5-openmpi-1.3.3-gnu64-4.1.2/lib
#IFCV64#     LDFLAGS       = $(LDFLAG)
#IFCV64#     FFLAGS        = -c -warn none
#IFCV64##     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DPOINTER64 -DIFC -DVECTOR
#IFCV64##     FFTLIB        = $(FFTWLIB)
#IFCV64#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DIFC -DVECTOR
#IFCV64#     FFTLIB        = $(FFTW3LIB)
##########################################FOR PC LINUX (Intel Fortran) on x86-64 at USC-HPC
#HPC#     LINKER        = mpif90
#HPC#     LDFLAGS       = $(LDFLAG)
#HPC#     FFLAGS        = -c -warn none
#HPC##     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DPOINTER64 -DIFC -DHPC -DVECTOR
#HPC##     FFTLIB        = $(FFTWLIB)
#HPC#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DIFC -DHPC -DVECTOR
#HPC#     FFTLIB        = $(FFTW3LIB)
########################################################FOR HP
#HP#      LINKER        = f90
#HP#      LDFLAGS       = +O2 +U77 -Wl,$(LDFLAG)
#HP#      FFLAGS        = -c +O2 +U77
#HP#      CPPDEFS       = $(CPPDEF) -DLIBFFTW -DHP
#HP#      FFTLIB        = $(FFTWLIB)
########################################################FOR IBM SP
#SP#      LINKER        = mpxlf
#SP#      LDFLAGS       = -O3 $(LDFLAG)
#SP#      FFLAGS        = -c -O3
#SP#      CPPDEFS       = $(CPPDEF) -DLIBFFTW -DIBM_SP
#SP#      FFTLIB        = $(FFTWLIB)
########################################################FOR SGI ORIGIN 2000
#ORIGIN#  LINKER        = f90
#ORIGIN#  LDFLAGS       = -static -O3 -mips4 -r10000 -OPT:Olimit=0 $(LDFLAG)
#ORIGIN#  FFLAGS        = -c -static -O3 -mips4 -r10000 \
#ORIGIN#                  -OPT:Olimit=0
#ORIGIN#  MPILIBS       = -lmpi
#ORIGIN#  CPPDEFS       = $(CPPDEF) -DLIBFFTW -DSGI_ORIGIN
#ORIGIN#  FFTLIB        = $(FFTWLIB)
########################################################FOR SGI ORIGIN 2600
#ORIGIN26# LINKER        = f90
##ORIGIN26#     LDFLAGS       = -n32 -mips4 -Ofast=ip27 -OPT:IEEE_arithmetic=3 $(LDFLAG)
#ORIGIN26# LDFLAGS       = -n32 -mips4 -O3 -OPT:Olimit=0 $(LDFLAG)
#ORIGIN26# FFLAGS        = -c $(LDFLAGS)
#ORIGIN26# MPILIBS       = -lmpi
#ORIGIN26# CPPDEFS       = $(CPPDEF) -DLIBFFTW -DSGI_ORIGIN
#ORIGIN26# FFTLIB        = $(FFTWLIB)
########################################################FOR CRAY T3E
#T3E#     LINKER        = f90
#T3E#     LDFLAGS       = -O3 -e v $(LDFLAG)
#T3E#     FFLAGS        = -c -O3 -M 1110 -e v
#T3E#     CPPDEFS       = $(CPPDEF) -DLIBFFTW -DCRAY_T3E
#T3E#     FFTLIB        = $(FFTWLIB)
########################################################FOR HITACHI SR2201
#SR2201#  LINKER        = f90
#SR2201#  LDFLAGS       = -L$(MPIR_HOME)/lib/hmpp2/cml -lfmpi -lmpi $(LDFLAG)
##SR2201#  FFLAGS        = -c -i,U -W0,'list(e(4))' \
##SR2201#                  -I$(MPIR_HOME)/include
#SR2201#  FFLAGS        = -c -i,U -W0,'list(e(4)),OPT(O(SS))' \
#SR2201#                  -I$(MPIR_HOME)/include
#SR2201#  CPPDEFS       = $(CPPDEF) -DLIBFFTW -DSR2201
#SR2201#  FFTLIB        = $(FFTWLIB)
########################################################FOR HITACHI SR8000
#SR8000#  LINKER        = xmpif90
##SR8000#  LDFLAGS       = -parallel $(LDFLAG)
#SR8000#  LDFLAGS       = $(LDFLAG)
##SR8000#  FFLAGS        = -c "-W0,'list(e(4)),MP(P(0)),OPT(O(SS))'"
#SR8000#  FFLAGS        = -c -Oss -noparallel "-W0,'list(e(4))'"
#SR8000#  CPPDEFS       = $(CPPDEF) -DLIBFFTW -DSR8000 -DVECTOR \
#SR8000#                  -DFAST_NONLOCAL_FORCE
#SR8000#  FFTLIB        = $(FFTWLIB)
########################################################FOR VPP 5000
## -X9 : compile sources as Fortran95 source
## -Am : create *.mod (module)
#VPP#     LINKER        = mpifrt
#VPP#     LDFLAGS       = 
#VPP#     FFLAGS        = -c -Sw -Aab -Kfast -X9 -Am
#VPP#     CPPDEFS       = $(CPPDEF) -DSSL2VP -DVPP -DVECTOR \
#VPP#                     -DFAST_NONLOCAL_FORCE
#VPP#     FFTLIB        = -lssl2vp
########################################################FOR SX7 (with ASL/SX library)
## -R2     : Specifies that transformation listing and formatted
##           listing are output.
## -ftrace : Specifies that the object file in execution format
##           corresponding to the ftrace command be generated
#SX7#     LINKER        = sxmpif90
#SX7#     LDFLAGS       = -ftrace
#SX7#     FFLAGS        = -c -R2 -ftrace -Wf,-pvctl loopcnt=9000
#SX7#     CPPDEFS       = -traditional $(CPPDEF) -DASLSX -DSX7 -DVECTOR \
#SX7#                     -DFAST_NONLOCAL_FORCE
#SX7#     FFTLIB        = -lasl
########################################################FOR SGI Altix3700/1280 at ISSP
######################################################## with Itanium & Intel Fortran
#ALTIX#     LINKER        = ifort
#ALTIX##     LDFLAGS       = -axW -static $(LDFLAG)
#ALTIX##     FFLAGS        = -c -axW -w90 -w95
#ALTIX#     LDFLAGS       = $(LDFLAG)
#ALTIX#     FFLAGS        = -c -w90 -w95
#ALTIX#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW -DIFC -DVECTOR -DPOINTER64
#ALTIX#     FFTLIB        = $(FFTWLIB)
#ALTIX#     MPILIBS       = -lmpi
########################################################FOR HITACHI SR11000 at ISSP
#SR11000#  LINKER        = mpixf90_r
#SR11000#  LDFLAGS       = $(LDFLAG)
#SR11000#  FFLAGS        = -c -O3 -noparallel -i,L -i,EU
#SR11000#  CPPDEFS       = $(CPPDEF) -DLIBFFTW -DSR11000 -DVECTOR \
#SR11000#                  -DFAST_NONLOCAL_FORCE
#SR11000#  FFTLIB        = $(FFTWLIB)
########################################################FOR HITACHI SR16000 at KEK
##SR16000#  LINKER        = mpif90
#SR16000#  LINKER        = mpxlf
##SR11000#  LDFLAGS       = $(LDFLAG)
#SR16000#  LDFLAGS       = -I/srhome/.tools/usr/local/fftw-3.3.3/include \
#SR16000#                  -L/srhome/.tools/usr/local/fftw-3.3.3/lib -O4 -qstrict
##SR16000#  FFLAGS        = -c -Os -64 -model=M1 -noparallel
#SR16000#  FFLAGS        = -c -O4 -qstrict
#SR16000#  CPPDEFS       = $(CPPDEF) -DLIBFFTW3 -DSR16000 -DVECTOR \
#SR16000#                  -DFAST_NONLOCAL_FORCE
#SR16000#  FFTLIB        = $(FFTW3LIB)
########################################################FOR IBM eServer model p5
#P5#  LINKER        = mpf90
#P5#  LDFLAGS       = $(LDFLAG) -O3 -qarch=pwr5 -qtune=pwr5
#P5#  FFLAGS        = -c -O3 -qarch=pwr5 -qtune=pwr5 -qfixed
#P5#  CPPDEFS       = $(CPPDEF) -DLIBFFTW -DPOINTER64 -DP5 -DVECTOR \
#P5#                  -DFAST_NONLOCAL_FORCE
#P5#  FFTLIB        = $(FFTWLIB)
##########################################FOR FUJITSU PRIMEQUEST 580
## -X9 : compile sources as Fortran95 source
## -Am : create *.mod (module)
#PRIMEQUEST#     LINKER        = mpifrt
#PRIMEQUEST#     LDFLAGS       = -Kfast
#PRIMEQUEST#     FFLAGS        = -c -Kfast -Am -X9 -f w,2004,2006,2008
#PRIMEQUEST#     CPPDEFS       = -traditional $(CPPDEF) -DSSL2VP -DPRIMEQUEST -DVECTOR
#PRIMEQUEST#     FFTLIB        = -SSL2
##########################################FOR FUJITSU PRIMERGY RX200S3
## -X9 : compile sources as Fortran95 source
## -Am : create *.mod (module)
## -pg : option only for PRIMERGY
#PRIMERGY#     LINKER        = mpifrt
#PRIMERGY#     LDFLAGS       = -pg -Kfast
#PRIMERGY#     FFLAGS        = -pg -c -Kfast -Am -X9 -f w,2004,2006,2008
#PRIMERGY#     CPPDEFS       = -traditional $(CPPDEF) -DSSL2VP -DPRIMERGY -DVECTOR
#PRIMERGY#     FFTLIB        = -SSL2
##########################################FOR FUJITSU PRIMEHPC FX10
## -X9 : compile sources as Fortran95 source
#PRIMEHPC#     LINKER        = mpifrtpx
#PRIMEHPC#     LDFLAGS       = -Kfast
#PRIMEHPC#     FFLAGS        = -c -Kfast -X9 -f w,2004,2006,2008
#PRIMEHPC#     CPPDEFS       = -traditional $(CPPDEF) -DSSL2VP -DPRIMEHPC -DVECTOR
#PRIMEHPC#     FFTLIB        = -SSL2
##########################################FOR FUJITSU PRIMERGY CX400
## -X9 : compile sources as Fortran95 source
#CX400#     LINKER        = mpifrt
#CX400#     LDFLAGS       = -Kfast
#CX400#     FFLAGS        = -c -Kfast -X9 -f w,2004,2006,2008
#CX400#     CPPDEFS       = -traditional $(CPPDEF) -DSSL2VP -DCX400 -DVECTOR
#CX400#     FFTLIB        = -SSL2
##########################################FOR FUJITSU PRIMERGY (ITO) at Kyushu Univ.
#CXITO#     LINKER        = mpifrt
#CXITO#     LDFLAGS       = -Kfast
#CXITO#     FFLAGS        = -c -Kfast -f w,2004,2006,2008
#CXITO#     CPPDEFS       = -traditional $(CPPDEF) -DSSL2VP -DVECTOR
#CXITO#     FFTLIB        = -SSL2
##########################################FOR FUJITSU PRIMERGY (ITO) at Kyushu Univ.
#CXITOFFTW#     LINKER        = mpifrt
#CXITOFFTW#     LDFLAGS       = -Kfast $(LDFLAG)
#CXITOFFTW#     FFLAGS        = -c -Kfast -f w,2004,2006,2008
#CXITOFFTW#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DVECTOR
#CXITOFFTW#     FFTLIB        = $(FFTW3LIB)
##########################################FOR HITACHI HA8000-tc/HT210
#HA8000#     LINKER        = mpiifort
##HA8000#     LDFLAGS       = -fast $(LDFLAG)
##HA8000#     FFLAGS        = -c -fast -warn none
##HA8000#     LDFLAGS       = -O3 -ipo -no-prec-div -xHOST $(LDFLAG)
##HA8000#     FFLAGS        = -c -O3 -ipo -no-prec-div -xHOST -warn none
#HA8000#     LDFLAGS       = $(LDFLAG)
#HA8000#     FFLAGS        = -c -warn none
#HA8000#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DIFC -DVECTOR
#HA8000#     FFTLIB        = $(FFTW3LIB)
##############################FOR Intel Fortran on SGI ICE XA/UV (systemB) at ISSP
#SEKIREI#     LINKER        = mpif90
#SEKIREI#     LDFLAGS       = -mkl=sequential -parallel
#SEKIREI#     FFLAGS        = -c -warn none -parallel -I $MKLROOT/include/fftw 
#SEKIREI#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DIFC -DVECTOR
##SEKIREI#     FFTLIB        = $(FFTW3LIB)
###########################################################################


FC		= $(LINKER)


SDIR		= ../Sources

IDIR		= ../Include


OBJS		= \
		allinc.o \
		modules.o \
		ftmain.o \
		constant.o \
		input.o \
		pwlda.o \
		pwreset.o \
		kpoints.o \
		bulk.o \
		stress.o \
		openfile.o \
		reclat.o \
		ba.o \
		eigen.o \
		schmidt.o \
		fermi.o \
		pcc.o \
		vxc.o \
		chgdns.o \
		funcs.o \
		ppkb.o \
		ppkbr.o \
		ppvd.o \
		nlkbpp.o \
		nlkbppr.o \
		force.o \
		tddft_fssh.o \
		sbalin.o \
		savedd.o \
		engrad.o \
		commwv.o \
		fft3.o \
		datadump.o \
		remd_modules.o \
		remd.o \
		remd_in.o


PROGRAM		= qxmd
PROGMPI		= qxmdmpi
PROGPVM		= qxmdpvm


## For the free-source form, comment the next line.
#SEDDEF		= | sed 's/&$$//' >

#.F90.o:
#	$(FC) $(FFLAGS) $(CPPDEFS) $*.F90
#
#.F90.o:
#	@rm -f $*.f90
#	/lib/cpp -P -C $(CPPDEFS) $*.F90 $(SEDDEF) $*.f90
#	$(FC) $(FFLAGS) $*.f90
#	rm  $*.f90
.F90.o:
	@rm -f $*.$(EXTF90)
#	/lib/cpp -P -C $(CPPDEFS) $*.F90 $(SEDDEF) $*.$(EXTF90)
	/lib/cpp -P $(CPPDEFS) $*.F90 $(SEDDEF) $*.$(EXTF90)
#GNU#	mv $*.$(EXTF90) $*.f90
	$(FC) $(FFLAGS) $*.f90
	rm -f $*.$(EXTF90)




$(PROGRAM):	$(OBJS) no_comm.o
	@echo "Loading $(PROGRAM) ... "
	$(LINKER) $(LDFLAGS) $(OBJS) no_comm.o -o $(PROGRAM) \
		$(LIBS) $(FFTLIB)
	mv $(PROGRAM) ..
	@echo "done"

$(PROGMPI):	$(OBJS) mpi_comm.o
	@echo "Loading $(PROGMPI) ... "
	$(LINKER) $(LDFLAGS) $(OBJS) mpi_comm.o -o $(PROGMPI) \
		$(MPILIBS) $(FFTLIB)
	mv $(PROGMPI) ..
	@echo "done"

$(PROGPVM):	$(OBJS) pvm_comm.o
	@echo "Loading $(PROGMPI) ... "
	$(LINKER) $(LDFLAGS) $(OBJS) pvm_comm.o -o $(PROGPVM) \
		$(PVMLIBS) $(FFTLIB)
	mv $(PROGPVM) ..
	@echo "done"


###############
# maintenance #
###############
clean:
		@rm -f $(OBJS)
		@rm -f *.o *.mod *.MOD

allclean:
		@rm -f $(OBJS)
		@rm -f *.o *.mod *.MOD
		@rm -f *.F90 *.h *.i


modules.o:     allinc.o
openfile.o:    allinc.o modules.o remd_modules.o
input.o:       allinc.o constant.o modules.o
pwlda.o:       allinc.o modules.o constant.o
pwreset.o:     allinc.o modules.o constant.o
kpoints.o:     allinc.o ftmain.o modules.o
bulk.o:        allinc.o
stress.o:      allinc.o modules.o
reclat.o:      allinc.o
ba.o:          allinc.o
eigen.o:       allinc.o modules.o
schmidt.o:     allinc.o modules.o
fermi.o:       allinc.o
pcc.o:         allinc.o commwv.o ppkb.o
vxc.o:         allinc.o modules.o
chgdns.o:      allinc.o modules.o
ppkb.o:        allinc.o commwv.o modules.o
ppkbr.o:       allinc.o ppkb.o
ppvd.o:        allinc.o ppkb.o
nlkbpp.o:      allinc.o commwv.o ppkb.o
nlkbppr.o:     allinc.o ppkbr.o nlkbpp.o
force.o:       allinc.o modules.o
tddft_fssh.o:  allinc.o constant.o
sbalin.o:      allinc.o modules.o
savedd.o:      allinc.o modules.o
engrad.o:      allinc.o modules.o
commwv.o:      allinc.o
fft3.o:        fftw_f77.i fftw3.f03
datadump.o:    allinc.o modules.o constant.o
remd_modules.o: allinc.o
remd_in.o:      constant.o
remd.o:         remd_modules.o constant.o


allinc.F90:       $(SDIR)/allinc.f90       ; cp $(SDIR)/allinc.f90 ./allinc.F90
modules.F90:      $(SDIR)/modules.f90      ; cp $(SDIR)/modules.f90 ./modules.F90
ftmain.F90:       $(SDIR)/ftmain.f90       ; cp $(SDIR)/ftmain.f90 ./ftmain.F90
constant.F90:     $(SDIR)/constant.f90     ; cp $(SDIR)/constant.f90 ./constant.F90
input.F90:        $(SDIR)/input.f90        ; cp $(SDIR)/input.f90 ./input.F90
pwlda.F90:        $(SDIR)/pwlda.f90        ; cp $(SDIR)/pwlda.f90 ./pwlda.F90
pwreset.F90:      $(SDIR)/pwreset.f90      ; cp $(SDIR)/pwreset.f90 ./pwreset.F90
kpoints.F90:      $(SDIR)/kpoints.f90      ; cp $(SDIR)/kpoints.f90 ./kpoints.F90
tddft_fssh.F90:   $(SDIR)/tddft_fssh.f90   ; cp $(SDIR)/tddft_fssh.f90 ./tddft_fssh.F90
openfile.F90:     $(SDIR)/openfile.f90     ; cp $(SDIR)/openfile.f90 ./openfile.F90
reclat.F90:       $(SDIR)/reclat.f90       ; cp $(SDIR)/reclat.f90 ./reclat.F90
ba.F90:           $(SDIR)/ba.f90           ; cp $(SDIR)/ba.f90 ./ba.F90
eigen.F90:        $(SDIR)/eigen.f90        ; cp $(SDIR)/eigen.f90 ./eigen.F90
schmidt.F90:      $(SDIR)/schmidt.f90      ; cp $(SDIR)/schmidt.f90 ./schmidt.F90
fermi.F90:        $(SDIR)/fermi.f90        ; cp $(SDIR)/fermi.f90 ./fermi.F90
pcc.F90:          $(SDIR)/pcc.f90          ; cp $(SDIR)/pcc.f90 ./pcc.F90
vxc.F90:          $(SDIR)/vxc.f90          ; cp $(SDIR)/vxc.f90 ./vxc.F90
chgdns.F90:       $(SDIR)/chgdns.f90       ; cp $(SDIR)/chgdns.f90 ./chgdns.F90
funcs.F90:        $(SDIR)/funcs.f90        ; cp $(SDIR)/funcs.f90 ./funcs.F90
ppkb.F90:         $(SDIR)/ppkb.f90         ; cp $(SDIR)/ppkb.f90 ./ppkb.F90
ppkbr.F90:        $(SDIR)/ppkbr.f90        ; cp $(SDIR)/ppkbr.f90 ./ppkbr.F90
ppvd.F90:         $(SDIR)/ppvd.f90         ; cp $(SDIR)/ppvd.f90 ./ppvd.F90
nlkbpp.F90:       $(SDIR)/nlkbpp.f90       ; cp $(SDIR)/nlkbpp.f90 ./nlkbpp.F90
nlkbppr.F90:      $(SDIR)/nlkbppr.f90      ; cp $(SDIR)/nlkbppr.f90 ./nlkbppr.F90
force.F90:        $(SDIR)/force.f90        ; cp $(SDIR)/force.f90 ./force.F90
bulk.F90:         $(SDIR)/bulk.f90         ; cp $(SDIR)/bulk.f90 ./bulk.F90
stress.F90:       $(SDIR)/stress.f90       ; cp $(SDIR)/stress.f90 ./stress.F90
sbalin.F90:       $(SDIR)/sbalin.f90       ; cp $(SDIR)/sbalin.f90 ./sbalin.F90
savedd.F90:       $(SDIR)/savedd.f90       ; cp $(SDIR)/savedd.f90 ./savedd.F90
engrad.F90:       $(SDIR)/engrad.f90       ; cp $(SDIR)/engrad.f90 ./engrad.F90
commwv.F90:       $(SDIR)/commwv.f90       ; cp $(SDIR)/commwv.f90 ./commwv.F90
fft3.F90:         $(SDIR)/fft3.f90         ; cp $(SDIR)/fft3.f90 ./fft3.F90
datadump.F90:     $(SDIR)/datadump.f90     ; cp $(SDIR)/datadump.f90 ./datadump.F90
no_comm.F90:      $(SDIR)/no_comm.f90      ; cp $(SDIR)/no_comm.f90 ./no_comm.F90
mpi_comm.F90:     $(SDIR)/mpi_comm.f90     ; cp $(SDIR)/mpi_comm.f90 ./mpi_comm.F90
pvm_comm.F90:     $(SDIR)/pvm_comm.f90     ; cp $(SDIR)/pvm_comm.f90 ./pvm_comm.F90
fftw_f77.i:       $(IDIR)/fftw_f77.i       ; cp $(IDIR)/fftw_f77.i ./fftw_f77.i
fftw3.f03:        $(IDIR)/fftw3.f03        ; cp $(IDIR)/fftw3.f03 ./fftw3.f03
remd_modules.F90: $(SDIR)/remd_modules.f90 ; cp $(SDIR)/remd_modules.f90 ./remd_modules.F90
remd.F90:         $(SDIR)/remd.f90         ; cp $(SDIR)/remd.f90 ./remd.F90
remd_in.F90:      $(SDIR)/remd_in.f90      ; cp $(SDIR)/remd_in.f90 ./remd_in.F90
