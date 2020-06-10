COMPILER 	= clang++ -O3
CFLAGS 		= 
LIBS 		= -lgsl -lgslcblas
INCLUDE 	= 
EXECUTABLE  	= test
OBJECT  	= main.o  ToroidalP.o

# Linking everything

all:		$(EXECUTABLE)

$(EXECUTABLE):	$(OBJECT)
		$(COMPILER) -o $(EXECUTABLE) $(OBJECT) $(LIBS)

%.o:		%.cpp 
		$(COMPILER) $(CFLAGS) -o $*.o -c $(INCLUDE) $*.cpp

# Cleanup

clean: 
		rm -f $(OBJECT)

clear:
		rm -f *.out *.eps 

