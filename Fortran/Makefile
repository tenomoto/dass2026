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
# Windows
#LIBS = -L" C:\rtools45\x86_64-w64-mingw32.static.posix\lib" -llapack -lblas
# Mac
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
run_module.o : settings_module.o
var.o : $(OBJS)
enkf.o : $(OBJS)
