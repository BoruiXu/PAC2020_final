TARGET = Demo

SOURCES = Demo.cpp NUFFT3D.cpp
HEADERS =
OBJECTS = $(SOURCES:.cpp=.o)

CXX = icpc
FLAGS = -std=c++11 -O3 -xHost -Wall -fopenmp -fp-model precise -DTIMING
LIB = -fopenmp -L/path/to/fftw/lib -lfftw3f -lfftw3f_threads
INC = -I/path/to/fftw/include

################################################################################
################################################################################
################################################################################

%.o : %.cpp
	$(CXX) -c $(SOURCES) $(INC) $(FLAGS)

$(TARGET) : $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(LIB)

.PHONY: clean
clean:	
	rm -rf $(OBJECTS) $(TARGET)

