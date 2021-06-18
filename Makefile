CXX = dpcpp
CXXFLAGS = -O2 -g -std=c++17

SWM_EXE_NAME = swm_dpcpp
SWM_SOURCES = shallow_unroll.cpp wtime.cpp


all: build_swm

build_swm:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME) $(SWM_SOURCES)

run: 
	./$(SWM_EXE_NAME)

clean: 
	rm -f $(SWM_EXE_NAME)
