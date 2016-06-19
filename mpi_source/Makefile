mpi_jayci.x: progheader.o mpi_jayci.o
	mpicc -I./include -o $@ progheader.o mpi_jayci.o

progheader.o:progheader.c
	cc -c -I./include progheader.c

mpi_jayci.o:mpi_jayci.c
	mpicc -c -I./include mpi_jayci.c

clean:
	rm -rf ./*.o
