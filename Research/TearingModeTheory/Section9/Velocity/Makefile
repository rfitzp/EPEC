# Makefile for program Velociy

COMPILER 	= clang++ -O3
CFLAGS 		= -w
LIBS		= -lblitz -L/opt/local/lib -lgsl -lgslcblas
INCLUDE	 	= -I/opt/local/include
EXECUTABLE	= velocity
OBJECT 		= Velocity.o

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

