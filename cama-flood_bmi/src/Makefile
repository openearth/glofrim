include   ../adm/Mkinclude
MODS = -I../mod -I../lib
LIBS = ../mod/mods.a ../lib/libs.a

.DEFAULT_GOAL := all

OBJECTS=\
calc_fldstg.o \
calc_dikstg.o \
calc_rivout.o \
calc_fldout.o \
calc_pthout.o \
calc_damout.o \
calc_outpre.o \
calc_rivout_kine.o \
calc_stonxt.o \
calc_watbal.o \
create_outcdf.o \
create_outbin.o \
control_inp.o \
control_inp_bmi.o \
control_bound.o\
control_phy.o \
control_out.o \
control_rest.o \
control_tstp.o \
init_inputnam.o \
init_map.o \
init_topo.o \
init_cond.o \
init_time.o \
control0.o

# these are the objects that should go into the shared library
SHARED_OBJECTS=\
	external_bmi.o

STATIC = srcs.a
$(STATIC): $(OBJECTS)
	ar -rv srcs.a $(OBJECTS)

.SUFFIXES : .o .f90 .F
.F.f:
	$(CPP) -P -I ./ $*.F > $*.f90
.F.o:
	$(CPP) -P -I ./ $*.F > $*.f90
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
.f90.o:
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o


# The shared library
LIBCAMO = libcama.so
$(LIBCAMO): $(OBJECTS) $(SHARED_OBJECTS) $(STATIC)
	$(FC) $(FFLAGS) $(LFLAGS) -shared $(SHARED_OBJECTS) $(STATIC) $(LIBS) -o $@ $(LIB)

# The executable
TARGET = MAIN_day
$(TARGET): $(OBJECTS) $(STATIC) $(TARGET).o
	$(FC) $(FFLAGS) $(LFLAGS) $@.o $(STATIC) $(LIBS) -o $@ $(LIB)

all: $(TARGET) $(LIBCAMO)

clean:
	 ${RM} -f *.o core *~ *trace temp* *.mod *.a *.f90 *.s MAIN_day *.so
