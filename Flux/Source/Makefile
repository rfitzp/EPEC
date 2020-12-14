F_OBJECT = NameListRead.o gFileRead.o gFileInterpolate.o
F_SOURCE = NameListRead.f90 gFileRead.f90 gFileInterpolate.f90

C_HEADER = Flux.h
C_OBJECT = main.o Flux.o Stage1.o Stage2.o Integrate.o\
	Interpolate.o Green.o ToroidalP.o gFileInterp.o 
C_SOURCE = main.cpp Flux.cpp Stage1.cpp Stage2.cpp Integrate.cpp\
	Interpolate.cpp Green.cpp ToroidalP.cpp gFileInterp.cpp

GIT_HASH             = `git rev-parse HEAD`
COMPILE_TIME         = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH           = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

include ../../Machines.in

EXECUTABLE = ../flux

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)

