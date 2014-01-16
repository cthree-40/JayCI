#


OBJS	:= orthogroutines.o test.o
LIBS	:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm
EXEC	:= test.x

test:$(OBJS)
	sh -c "ifort -i8 -o  $(EXEC) $(OBJS) $(LIBS) blaswrapper.a"

%.o:%.f90
	sh -c "ifort -cpp -o $@ $< "

clean:
	rm -r ./*.o ./*.mod
