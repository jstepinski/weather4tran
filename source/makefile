EXEC = weather

FFLAGS = -g -fdefault-real-8 -fdefault-integer-8 \
        -Wall -Wextra -fbounds-check \
        -ffree-line-length-0

SRC = uvars.f90 pvars.f90 funcs.f90 xfun.f90 \
      readTableB1.f90 \
      interp.f90 find2els.f90 \
      extendMat.f90 \
      initializePsurf.f90 \
      calc_flux.f90 calc_pia.f90 calc_sigmaDot.f90 \
      calc_theta.f90 calc_phi.f90 calc_PP.f90 \
      calc_uwind.f90 calc_vwind.f90 \
      printToFile.f90 \
      weather.f90
OBJ = $(SRC:.f90=.o)

LIBS = -lm
FC = gfortran

.SUFFIXES: .f90 

all: $(SRC) $(EXEC)

$(EXEC):  $(OBJ)
	$(FC) $^ $(LIBS) -o $@

%.o: %.mod

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm *.o
	rm *.mod
