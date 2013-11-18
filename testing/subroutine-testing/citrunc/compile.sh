#!/usr/bin/bash

ifort -c detci1.f90 -g -traceback -fp-stack-check -check bounds
ifort -c detci2.f90 -g -traceback -fp-stack-check -check bounds
ifort -c truncation.f90 -g -traceback -fp-stack-check -check bounds
ifort -c citrunc.f90 -g -traceback -fp-stack-check -check bounds
ifort -c test.f90 -g -traceback -fp-stack-check -check bounds

ifort -o test.x detci1.o detci2.o truncation.o citrunc.o test.o

rm ./*.o ./*.mod
rm ./*.dets ./*.list
