FC=gfortran -c
LD=gfortran
SRC=params.f90 allocatearrays.f90 deallocatearrays.f90 scheightorder.f90
OBJ=params.o allocatearrays.o deallocatearrays.o scheightorder.o
schrodinger:
	$(FC) $(SRC)
	$(LD) $(OBJ) -o scheightorder.x -llapack -lblas
	rm -rf *.o
clean:
	rm -rf *.o *.mod *.x
