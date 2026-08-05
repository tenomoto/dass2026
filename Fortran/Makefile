.SUFFIXES:
.SUFFIXES: .o .f90

.ALLOBJS ?= $^

FC = gfortran
LDLIBS = 
LDFLAGS =
FFLAGS = -O2
SRCS = random_module.f90 settings_module.f90 eqocean_module.f90 run_module.f90
OBJS = $(SRCS:.f90=.o)
MODS = $(SRCS:.f90=.mod)
#LIBS = -lopenblas
LIBS = -framework accelerate

TARGETS = check_adjoint var enkf

all: $(TARGETS)

check_adjoint : $(OBJS) check_adjoint.o
	$(FC) $(FFLAGS) $(LDFLAGS) $(.ALLOBJS) $(LDLIBS) -o $@

var : $(OBJS) var.o
	$(FC) $(FFLAGS) $(LDFLAGS) $(.ALLOBJS) $(LDLIBS) -o $@

enkf : $(OBJS) matrix_module.o enkf.o
	$(FC) $(FFLAGS) $(LDFLAGS) $(.ALLOBJS) $(LDLIBS) $(LIBS) -o $@

clean:
	rm -f $(OBJS) $(MODS) $(TARGETS) matrix_module.o matrix_module.mod

.f90.o:
	$(FC) $(FFLAGS) -c $<

check_adjoint.o : $(OBJS)
var.o : $(OBJS)
enkf.o : $(OBJS)
