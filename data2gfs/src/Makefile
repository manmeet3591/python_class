SH=/bin/sh
FC=ifort 
FFLAGS= -O3 -openmp -xhost -fp-model strict -assume byterecl -CB -traceback
INCDIR=/home/miyachi/ncepgfs/lib/incmod
LIBDIR=/home/miyachi/ncepgfs/lib
INCS=-I$(INCDIR)/sigio_4 -I$(INCDIR)/w3lib-2.0_4 -I$(INCDIR)/g2_4
LDLIBS=-L$(LIBDIR) -lip_4 -lsp_4 -lsigio_4 -lbacio_4 -lg2_4 -lw3lib-2.0_4 -lpng -ljasper

OBJ = type_module.o vcoord_def_module.o constant_module.o read_module.o calcmet_module.o ip_module.o p2p_module.o gfs_module.o param_module.o #  debug_module.o

TARGET = ../exec/data2gfs

.SUFFIXES:
.SUFFIXES: .f90 .o

$(TARGET):	$(OBJ) main.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) 
.f90.o:
	$(FC) $(INCS) $(FFLAGS) -c $< 

clean:
	rm -f *.mod *.o $(TARGET)

param_module.o : type_module.o
#debug_module.o : type_module.o
p2p_module.o:	type_module.o calcmet_module.o constant_module.o
ip_module.o:	type_module.o
gfs_module.o:	type_module.o 
extrapolate_module.o:	type_module.o
calcmet_module.o:	type_module.o constant_module.o
read_module.o:	type_module.o vcoord_def_module.o constant_module.o param_module.o ip_module.o
constant_module.o:	type_module.o
vcoord_def_module.o:	type_module.o
main.o:	$(OBJ)

