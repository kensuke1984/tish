# makefile for tish
PROGS = mpi-tish
FC	= mpif90
FC2=gfortran
FFLAGS	= -O

SRC = calmat.f90 trialf.f90 others.f90 dclisb.f90 \
	formpi.f90 mpi-tish.f90 parameters.f90
OBJS	= $(SRC:.f90=.o)
.SUFFIXES: .f90 .o

all:$(PROGS)

mpi-tish: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

mpi-tish.o: mpi-tish.f90
	$(FC) $(FFLAGS) -c mpi-tish.f90 -o $@

.f90.o:
	$(FC2) $(FFLAGS) -c $< 

%.o: %.mod

clean:
	rm -f $(OBJS) $(PROGS) $(OBJS_SINGLE) parameters.mod
