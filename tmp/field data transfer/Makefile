CXX := g++
CXXFLAGS:=  -std=c++11 -fopenmp -O2
LIBFLAGS:=   -lfftw3_omp -lfftw3 
object=datatrans.o


data: $(object)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

datatrans.o: 




clean:
	rm *.o 
