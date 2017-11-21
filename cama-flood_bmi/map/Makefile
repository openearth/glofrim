##################################################
include   ../adm/Mkinclude
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

TARGET2 = calc_outclm
$(TARGET2) : $(TARGET2).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)
	
TARGET3 = calc_rivwth
$(TARGET3) : $(TARGET3).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

TARGET4 = set_gwdlr
$(TARGET4) : $(TARGET4).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

TARGET5 = set_fldpth
$(TARGET5) : $(TARGET5).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

TARGET6 = set_rivhgt_pth
$(TARGET6) : $(TARGET6).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

##################################################
# Commands
##################################################
all : ${TARGET1} ${TARGET2} ${TARGET3} ${TARGET4} ${TARGET5} ${TARGET6}
clean: 
	$(RM) -f *.o *.f90 *.s core *~ *trace *.mod

