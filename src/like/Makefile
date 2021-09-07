FC = gfortran
FCFLAGS = -O3 -fopenmp -g -ffree-line-length-none -fbounds-check
OBJFILES = params.o neldermead.o util.o LabFuncs.o like.o
PROGRAMS = runDL runDL_2D

all: clean $(OBJFILES) $(PROGRAMS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

neldermead.o: neldermead.f
		$(FC) $(FCFLAGS) -c neldermead.f

$(PROGRAMS): $(OBJFILES)
	$(FC) $(FCFLAGS) $@.f95  -o $@ $^

clean:
	rm -f *.o *.mod
