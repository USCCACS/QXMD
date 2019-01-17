#################################################################
#                                                               #
#  Makefile for electronic structure calculation by PW method   #
#                                                               #
#                     : make machine                            #
#                                                               #
#                                            work directory     #
#    serial           : make qxmd              QXMD             #
#    MPI              : make qxmdmpi           QXMD             #
#                                                               #
#################################################################
SHELL = /bin/sh


MAKEFILE	= Makefile

SDIR		= Sources

WDIR		= QXMD
DDIR		= data



####################################################################
PROG		= qxmd
PROGMPI		= qxmdmpi


$(PROG):	$(WDIR)/$(MAKEFILE) $(WDIR)/$(PROG)

$(PROGMPI):	$(WDIR)/$(MAKEFILE) $(WDIR)/$(PROGMPI)


$(WDIR)/$(PROG):
	cd $(WDIR); make $(PROG)

$(WDIR)/$(PROGMPI):
	cd $(WDIR); make $(PROGMPI)


$(WDIR)/$(MAKEFILE):	$(SDIR)/$(MAKEFILE)
	@echo Execute: make machine
	@make help
	@exit 1


####################################################################
# Machine dependent 
SOURCEFILE = $(SDIR)/$(MAKEFILE)
WORKFILE   = $(WDIR)/$(MAKEFILE)

dec: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#DEC#//" $(SOURCEFILE) > $(WORKFILE)

gnu: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#GNU#//" $(SOURCEFILE) > $(WORKFILE)

ffc: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#FFC#//" $(SOURCEFILE) > $(WORKFILE)

ffc-v64: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#FFCV64#//" $(SOURCEFILE) > $(WORKFILE)

ifc: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#IFC#//" $(SOURCEFILE) > $(WORKFILE)

ifc-x86: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#IFCX86#//" $(SOURCEFILE) > $(WORKFILE)

ifc-v64: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#IFCV64#//" $(SOURCEFILE) > $(WORKFILE)

hpc: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#HPC#//" $(SOURCEFILE) > $(WORKFILE)

hp: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#HP#//" $(SOURCEFILE) > $(WORKFILE)

sp: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SP#//" $(SOURCEFILE) > $(WORKFILE)

origin: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#ORIGIN#//" $(SOURCEFILE) > $(WORKFILE)

origin26: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#ORIGIN26#//" $(SOURCEFILE) > $(WORKFILE)

t3e: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#T3E#//" $(SOURCEFILE) > $(WORKFILE)

sr2201: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SR2201#//" $(SOURCEFILE) > $(WORKFILE)

sr8000: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SR8000#//" $(SOURCEFILE) > $(WORKFILE)

sr11000: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SR11000#//" $(SOURCEFILE) > $(WORKFILE)

sr16000: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SR16000#//" $(SOURCEFILE) > $(WORKFILE)

vpp: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#VPP#//" $(SOURCEFILE) > $(WORKFILE)

sx7: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SX7#//" $(SOURCEFILE) > $(WORKFILE)

altix: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#ALTIX#//" $(SOURCEFILE) > $(WORKFILE)

p5: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#P5#//" $(SOURCEFILE) > $(WORKFILE)

primequest: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#PRIMEQUEST#//" $(SOURCEFILE) > $(WORKFILE)

primergy: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#PRIMERGY#//" $(SOURCEFILE) > $(WORKFILE)

primehpc: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#PRIMEHPC#//" $(SOURCEFILE) > $(WORKFILE)

cx400: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#CX400#//" $(SOURCEFILE) > $(WORKFILE)

cxito: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#CXITO#//" $(SOURCEFILE) > $(WORKFILE)

ha8000: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#HA8000#//" $(SOURCEFILE) > $(WORKFILE)

sekirei: $(WDIR) $(DDIR) $(SOURCEFILE)
	sed "s/^#SEKIREI#//" $(SOURCEFILE) > $(WORKFILE)


$(WDIR):
	mkdir $(WDIR)

$(DDIR):
	mkdir $(DDIR)

$(SDIR)/$(MAKEFILE):
	@echo cannot find $(SDIR)/$(MAKEFILE)
	@exit 1


####################################################################
help:
	@echo 'make dec        : DEC Alpha Degital Fortran'
	@echo 'make gnu        : GNU Fortran (only for classical MD)'
	@echo 'make ffc        : PC LINUX (Fujitsu Fortran&C)'
	@echo 'make ffc-v64    : PC LINUX (Fujitsu Fortran&C) on VT-64'
	@echo 'make ifc        : PC LINUX (Intel Fortran)'
	@echo 'make ifc-x86    : PC LINUX (Intel Fortran) on x86-64'
	@echo 'make ifc-v64    : PC LINUX (Intel Fortran) on VT-64'
	@echo 'make hpc        : PC LINUX (Intel Fortran) on x86-64 at USC-HPC'
	@echo 'make hp         : HP FORTRAN 90'
	@echo 'make sp         : IBM SP'
	@echo 'make origin     : SGI ORIGIN 2000'
	@echo 'make origin26   : SGI ORIGIN 2600'
	@echo 'make t3e        : CRAY T3E'
	@echo 'make sr2201     : HITACHI SR2201'
	@echo 'make sr8000     : HITACHI SR8000'
	@echo 'make sr11000    : HITACHI SR11000 at ISSP'
	@echo 'make sr16000    : HITACHI SR16000 at KEK'
	@echo 'make ha8000     : HITACHI HA8000-tc/HT210'
	@echo 'make vpp        : Fujitsu VPP5000'
	@echo 'make sx7        : NEC SX-7'
	@echo 'make altix      : SGI Altix3700/1280 at ISSP'
	@echo 'make p5         : IBM eServer model p5'
	@echo 'make primequest : FUJITSU PRIMEQUEST 580'
	@echo 'make primergy   : FUJITSU PRIMERGY RX200S3'
	@echo 'make cx400      : FUJITSU PRIMERGY CX400'
	@echo 'make cxito      : FUJITSU PRIMERGY CX (ITO) at Kyushu Univ.'
	@echo 'make primehpc   : FUJITSU PRIMEHPC FX10'
	@echo 'make sekirei    : Intel Fortran on SGI ICE XA/UV (systemB) at ISSP'

