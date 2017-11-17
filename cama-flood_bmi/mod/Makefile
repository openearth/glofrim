include   ../adm/Mkinclude

OBJECTS=\
parkind1.o \
mod_input.o \
mod_map.o \
mod_prog.o \
mod_diag.o \
mod_time.o \
mod_output.o \
iso_c_utils.o \
logging.o 

.SUFFIXES : .o .f90 .F
.F.f:
	$(CPP) -P -I ./ $*.F > $*.f90
.F.o:
	$(CPP) -P -I ./ $*.F > $*.f90
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
.f90.o:
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
all:    $(OBJECTS)
	ar -rv mods.a $(OBJECTS)

clean:
	 ${RM} -f *.o core *~ *trace temp* *.f90 *.s *.mod *.a
