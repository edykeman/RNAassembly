# Makefile for VAssemble
#
# Written: April 1st 2016
# Last Update: April 1st 2016
# Author: Eric C. Dykeman

F90 = f95
FFLAGS = -O4

# ---------------------------------------------------------------------
# OBJECTS
# ---------------------------------------------------------------------

MAIN = vassembly.o

MODULES = systemvar.o

CLASSES = class_vrna.o

INOUT = readdata.o

ENERGY = ebulge.o  ehair.o  estack.o  tint11.o  tint12.o  tint22.o \
         tloop.o  tstack.o  tstackh.o  tstacki.o rnarates.o

MISC = snreaction.o random.o rkappa.o

OBJECTS = $(MAIN) $(MODULES) $(CLASSES) $(INOUT) $(ENERGY) $(MISC)

# ---------------------------------------------------------------------
# MAKE COMMANDS
# ---------------------------------------------------------------------

all: vlc.x clean

vlc.x: $(OBJECTS)
	$(F90) $(FFLAGS) -o VLC.x $(OBJECTS)

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod *.x

# ---------------------------------------------------------------------
# COMPILE COMMANDS
# ---------------------------------------------------------------------

# --------- MAIN ------------

vassembly.o: vassembly.f90 $(OBJECTS)
	$(F90) $(FFLAGS) -c vassembly.f90

# -------- MODULES ----------

systemvar.o: systemvar.f90
	$(F90) $(FFLAGS) -c systemvar.f90

# -------- CLASSES ----------

class_vrna.o: class_vrna.f90 $(MODULES)
	$(F90) $(FFLAGS) -c class_vrna.f90

# --------- INOUT -----------

readdata.o: readdata.f90 $(MODULES)
	$(F90) $(FFLAGS) -c readdata.f90

# -------- ENERGY -----------

ebulge.o: ENERGY/ebulge.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ENERGY/ebulge.f90
ehair.o: ENERGY/ehair.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ENERGY/ehair.f90
estack.o: ENERGY/estack.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ENERGY/estack.f90
tint11.o: ENERGY/tint11.f90
	$(F90) $(FFLAGS) -c ENERGY/tint11.f90
tint12.o: ENERGY/tint12.f90
	$(F90) $(FFLAGS) -c ENERGY/tint12.f90
tint22.o: ENERGY/tint22.f90
	$(F90) $(FFLAGS) -c ENERGY/tint22.f90
tloop.o: ENERGY/tloop.f90
	$(F90) $(FFLAGS) -c ENERGY/tloop.f90
tstack.o: ENERGY/tstack.f90
	$(F90) $(FFLAGS) -c ENERGY/tstack.f90
tstackh.o: ENERGY/tstackh.f90
	$(F90) $(FFLAGS) -c ENERGY/tstackh.f90
tstacki.o: ENERGY/tstacki.f90
	$(F90) $(FFLAGS) -c ENERGY/tstacki.f90
rnarates.o: ENERGY/rnarates.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ENERGY/rnarates.f90

# --------- MISC ------------

snreaction.o: snreaction.f90 $(MODULES)
	$(F90) $(FFLAGS) -c snreaction.f90
rkappa.o: rkappa.f90
	$(F90) $(FFLAGS) -c rkappa.f90
random.o: random.f90
	$(F90) $(FFLAGS) -c random.f90
