CXX = dpcpp
CXXFLAGS = -O2 -g -std=c++17

SWM_EXE_NAME_DPCPP = swm_dpcpp
SWM_SOURCES_DPCPP = shallow_unroll.dpcpp.cpp wtime.cpp

SWM_EXE_NAME_CPP = swm_cpp
SWM_SOURCES_CPP = shallow_unroll.cpp wtime.cpp

SWM_EXE_NAME_ACC_OMP = swm_acc_omp
SWM_SOURCES_ACC_OMP = shallow_unroll.acc.omp.cpp wtime.cpp

all: swm_dpcpp

swm_dpcpp:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_DPCPP) $(SWM_SOURCES_DPCPP)

swm_cpp:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_CPP) $(SWM_SOURCES_CPP)
    
swm_acc_omp:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_ACC_OMP) $(SWM_SOURCES_ACC_OMP)

clean: 
	rm -f $(SWM_EXE_NAME_DPCPP) $(SWM_EXE_NAME_CPP) $(SWM_EXE_NAME_ACC_OMP) 
