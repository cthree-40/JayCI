## Makefile



OBJS			:= detci1.o detci2.o truncation.o citrunc.o test.o

PREFIX		:= sh -c
COMPILER	:= ifort -i8
EXEC			:= test.x

DEBUG			:= -g -traceback -fp-stack-check -check bounds

RM				:= rm -rf

test:$(OBJS)
	@echo "BUILDING"
	$(PREFIX) "$(COMPILER) -o $(EXEC) $(OBJS)"

%.o:./src/%.f90
	$(PREFIX) "$(COMPILER) -cpp -c $< -O2"

clean:
	$(RM) $(OBJS) ./*.list ./*.dets ./*.mod
