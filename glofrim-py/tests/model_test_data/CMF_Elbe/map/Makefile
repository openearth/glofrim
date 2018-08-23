##################################################
include   
#
#       Macro (Settings)
##################################################
.SUFFIXES : .o .f90 .F
.F.f:
	$(CPP) -P -I ./ $*.F > $*.f90
.F.o:
	$(CPP) -P -I ./ $*.F > $*.f90
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
.f90.o:
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
src: $(OBJECTS)
	ar -rv srcs.a $(OBJECTS)
##################################################
# simulation  
##################################################
TARGET1 = generate_inpmat
$(TARGET1) : $(TARGET1).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)


##################################################
# Commands
##################################################
all : ${TARGET1}
clean: 
	$(RM) -f *.o *.f90 *.s core *~ *trace *.mod

