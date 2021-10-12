# Makefile for Program NEOCLASSICAL

F_OBJECT = NamelistRead.o
F_SOURCE = NamelistRead.f90

C_HEADER = Neoclassical.h
C_OBJECT = main.o Neoclassical.o Field.o pFileRead.o Interpolate.o pFileInterpolate.o fFileInterpolate.o cFileRead.o cFileInterpolate.o Smoothing.o
C_SOURCE = main.cpp Neoclassical.cpp Field.cpp pFileRead.cpp Interpolate.cpp pFileInterpolate.cpp fFileInterpolate.cpp cFileRead.cpp cFileInterpolate.cpp Smoothing.cpp

GIT_HASH             = `git rev-parse HEAD`
COMPILE_TIME         = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH           = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

include ../../Machines.in

EXECUTABLE = ../neoclassical

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE)  $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)
