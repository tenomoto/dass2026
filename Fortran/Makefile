.SUFFIXES:
.SUFFIXES: .o .f90
FC = gfortran
LDLIBS = 
LDFLAGS =
FFLAGS = -O2
SRCS = random_module.f90 settings_module.f90 eqocean_module.f90
OBJS = $(SRCS:.f90=.o)

check_adjoint : $(OBJS) check_adjoint.o
	$(FC) $(FFLAGS) $(LDFLAGS) $? $(LDLIBS) -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.f90.o:
	$(FC) $(FFLAGS) -c $<

check_adjoint.o : $(OBJS)
var.o : $(OBJS)
enkf.o : $(OBJS)
