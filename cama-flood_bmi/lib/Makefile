include   ../adm/Mkinclude

OBJECTS=\
lib_dates.o \
lib_netcdf_util.o \
len_trim.o\
endian4.o\
igetday.o

MODS = -I../mod

.SUFFIXES : .o .f90 .F
.F.f:
	$(CPP) -P -I ./ $*.F > $*.f90
.F.o:
	$(CPP) -P -I ./ $*.F > $*.f90
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
.f90.o:
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
all:    $(OBJECTS)
	ar -rv libs.a $(OBJECTS)

clean:
	 ${RM} -f *.o core *~ *trace temp* *.f90 *.s *.mod *.a 

