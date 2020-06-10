F_OBJECT = NamelistRead.o
F_SOURCE = NamelistRead.f90

C_HEADER = Phase.h
C_OBJECT = main.o Phase.o uFileInterpolate.o lFileInterpolate.o fFileInterpolate.o nFileInterpolate.o
C_SOURCE = main.cpp Phase.cpp uFileInterpolate.cpp lFileInterpolate.cpp fFileInterpolate.cpp nFileInterpolate.cpp

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	F_COMPILER = gfortran
	F_FLAGS    =
	F_INCLUDE  =

	C_COMPILER = clang++ -O3
	C_FLAGS    = -w
	C_INCLUDE  = -I/opt/local/include
	C_LIBS	   = -L/opt/local/lib -L/opt/local/lib/libgcc -lgfortran.5 -lgsl -lgslcblas
else ifeq ($(UNAME), Linux)
	F_COMPILER = gfortran
	F_FLAGS    =
	F_INCLUDE  = 

	C_COMPILER = clang++ -O3
	C_FLAGS    = -w
	C_INCLUDE  = -lgsl -lgslblas
	C_LIBS	   = -lgfortran -lgsl -lgslcblas
endif

EXECUTABLE = ../phase

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE)  $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)