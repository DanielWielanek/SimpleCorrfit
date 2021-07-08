CC    = g++
LD    = gfortran
FF    = gfortran
#COPTS    = `root-config --cflags` -I/usr/local/root/include -g
#LDOPTS    = `root-config --libs` -g
COPTS    = -g -O  `root-config --cflags` -Wall -g
LDOPTS    = -g -O  `root-config --libs` -lstdc++ -Wall -g
#C++ Files
SOURCES =  Main.cpp Pair.cpp StandAloneFsiLednicky.cpp StandAloneSimpleFsi.cpp
OBJECTS = $(SOURCES:.cpp=.o)
#FORTRAN Files
FSOURCES= FsiTools.F FsiWeightLednicky4.F
FOBJECTS=$(FSOURCES:.F=.o)
EXECUTABLE=main.exe
all: $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS) $(FOBJECTS)
	$(LD) -o $@ $^ $(LDOPTS)
#C++ files
.cpp.o:
	$(CC) -o $@ $^ -c $(COPTS)
#FORTRAN files
.F.o:
	$(FF) -O3 -lgfortran -o $@ $^ -c $(COPTS)
clean:;         @rm -f $(OBJECTS)  $(EXECUTABLE) *.o  *.d