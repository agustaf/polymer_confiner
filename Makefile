
CXX = g++
CXX_FLAGS = -std=c++17 -O3
OUTPUTNAME = confiner.exe

all:
	$(CXX) $(CXX_FLAGS) -o $(OUTPUTNAME) polymer_confiner.cpp

clean:
	rm -f $(OUTPUTNAME)
