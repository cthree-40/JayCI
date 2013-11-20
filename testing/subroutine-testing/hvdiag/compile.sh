#!/usr/bin/bash

ifort -c -cpp detci1.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c -cpp detci2.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c -cpp detci5.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c -cpp truncation.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c -cpp citrunc.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c -cpp hvdiag.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c iwfmt.f -g -traceback -fp-stack-check -check bounds -DHVDIAG
ifort -c -cpp test.f90 -g -traceback -fp-stack-check -check bounds -DHVDIAG

ifort -o test.x detci1.o detci2.o detci5.o truncation.o citrunc.o hvdiag.o iwfmt.o test.o colib.a blaswrapper.a

rm ./*.o ./*.mod
rm ./*.dets ./*.list
rm ./*.ref
