# Makefile for program Integrals

COMPILER 	= clang++ -O3
CFLAGS 		= -w
LIBS		= -lblitz -L/opt/local/lib -lgsl -lgslcblas
INCLUDE	 	= -I/opt/local/include
EXECUTABLE	= integrals
OBJECT 		= main.o Integrals.o 

# Linking everything:

$(EXECUTABLE): $(OBJECT)
	       $(COMPILER) -o $(EXECUTABLE) $(OBJECT) $(LIBS)

%.o:	%.cpp %.h
	$(COMPILER) $(CFLAGS) -o $*.o -c $(INCLUDE) $*.cpp

# Cleanup

clean: 
	  rm -f $(EXECUTABLE) *.o

clear:
	  rm -f  *.o

