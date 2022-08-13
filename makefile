CXX = gfortran
OBJ = subroutines_matrix.o constants.o laguerrepol.o har_oscillator.o eigenvalues.o integration.o hamiltonian_elements.o variational.o

all: $(OBJ) main.exe

main.exe: $(OBJ) main.f90 
	$(CXX) -o $@ $^

%.o: %.f90
	$(CXX) -c $<


.PHONY: clean
clean:
	del *.o *~ *.exe *.mod *.a *#


