.SUFFIXES: .o .f90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
INC  = -I/usr/include
FFLAGS = -O2 -ffree-line-length-none #-Wunused 
OBJS = prepmegan4cmaq.o
EXE  = prepmegan4cmaq.exe

.f90.o:
		${FC} ${FFLAGS} -c ${INC} $<

prepmegan4cmaq: ${OBJS}
		${FC} -o ${EXE} ${OBJS} ${LIBS} 

clean:
		rm -f *.o *.mod
