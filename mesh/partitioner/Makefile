#-----------------------------------------------------------------------
#   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
#   for Information Science & Technology (RIST)
#-----------------------------------------------------------------------

# 
# makefile for GeoFEM partitioner
#

MPICHBINDIR = /opt/local/mpich2/bin
MPICHLIBDIR = /opt/local/mpich2/lib
MPICHINCDIR = /opt/local/mpich2/include
INSTDIR     = 
CC          = /opt/local/mpich2/bin/mpicc
F90         = gfortran
MPILIBS     = -lfmpich -lmpich
RANLIB      = ranlib
AR          = ar
F90OPTFLAGS = -O
OPTFLAGS    = -O

UTILDIR     = ../geoutil

#
# Compile flags
#
CFLAGS   =$(OPTFLAGS) -I$(UTILDIR) -I$(UTILDIR)/file -I$(MPICHINCDIR)
F90FLAGS =$(F90OPTFLAGS) -I. -I$(UTILDIR) -I$(UTILDIR)/file -I$(MPICHINCDIR)
F90LFLAGS=-L$(MPICHLIBDIR) $(MPILIBS)

.SUFFIXES:
.SUFFIXES: .o .f90 .c
#
.f90.o:; $(F90) -c $(F90FLAGS)  $(F90OPTFLAG) $<
#
OBJS = \
partitioner.o analyzer.o geofem_tiger.o \
proc_local.o local_data.o calc_edgcut.o cre_local_data.o \
interface_nodes.o create_metis_input.o\
neib_pe.o paraset.o rcm_gibbs.o util.o

TARGET = partGeoFEM

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(F90) $(F90FLAGS) -o $(TARGET) \
	$(OBJS) \
	-L$(UTILDIR)/file -lGeoUtil \
	$(F90LFLAGS)

analyzer.o: analyzer.f90  partitioner.o partitioner.o
partitioner.o: partitioner.f90
local_data.o: local_data.f90 partitioner.o
create_metis_input.o: create_metis_input.f90 partitioner.o
interface_nodes.o: interface_nodes.f90 partitioner.o partitioner.o
multicolor.o: multicolor.f90
geofem_tiger.o: geofem_tiger.f90  analyzer.o partitioner.o partitioner.o
prep_next_partition.o: prep_next_partition.f90 partitioner.o
calc_edgcut.o: calc_edgcut.f90 partitioner.o
paraset.o: paraset.f90 partitioner.o
neib_pe.o: neib_pe.f90 partitioner.o partitioner.o
reduce_edgcut.o: reduce_edgcut.f90 partitioner.o
rhoset.o: rhoset.f90 partitioner.o
ordering.o: ordering.f90 partitioner.o
partition.o: partition.f90 partitioner.o
cre_local_data.o: cre_local_data.f90 partitioner.o
proc_local.o: proc_local.f90 partitioner.o

clean:
	rm -f *.o $(TARGET) *.mod *~ PI*

install:
	cp $(TARGET) $(INSTDIR)/bin
