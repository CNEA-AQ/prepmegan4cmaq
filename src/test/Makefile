.SUFFIXES: .o .f90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
INC  = -I/usr/include
FCFLAGS = -O2 -ffree-line-length-none #-Wunused 
LDFLAGS = -O2
EXE  = test.exe
OBJS = module_utils.o module_proj.o module_griddesc_reader.o module_nc_handler.o test.o  #module_interp.o 

prepmegan4cmaq: $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@

clean:
		rm -f *.o *.mod
