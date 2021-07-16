CXX = dpcpp
#Casper
#CXXFLAGS = -O2 -g --gcc-toolchain=/glade/u/apps/dav/opt/gnu/9.1.0 -std=c++17
#Intel DevCloud
CXXFLAGS = -O2 -g -std=c++17

SWM_EXE_NAME_DPCPP_USM = swm_dpcpp_usm
SWM_SOURCES_DPCPP_USM = shallow_unroll.dpcpp_usm.cpp wtime.cpp

SWM_EXE_NAME_DPCPP_USM_SHARED = swm_dpcpp_usm_shared
SWM_SOURCES_DPCPP_USM_SHARED = shallow_unroll.dpcpp_usm_shared.cpp wtime.cpp

SWM_EXE_NAME_DPCPP_BUF = swm_dpcpp_buf
SWM_SOURCES_DPCPP_BUF = shallow_unroll.dpcpp_buf.cpp wtime.cpp

SWM_EXE_NAME_CPP = swm_cpp
SWM_SOURCES_CPP = shallow_unroll.cpp wtime.cpp

SWM_EXE_NAME_ACC_OMP = swm_acc_omp
SWM_SOURCES_ACC_OMP = shallow_unroll.acc.omp.cpp wtime.cpp

all: swm_dpcpp_usm

swm_dpcpp_usm:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_DPCPP_USM) $(SWM_SOURCES_DPCPP_USM)
    
swm_dpcpp_usm_shared:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_DPCPP_USM_SHARED) $(SWM_SOURCES_DPCPP_USM_SHARED)

swm_dpcpp_buf:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_DPCPP_BUF) $(SWM_SOURCES_DPCPP_BUF)
    
swm_cpp:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_CPP) $(SWM_SOURCES_CPP)
    
swm_acc_omp:
	$(CXX) $(CXXFLAGS) -o $(SWM_EXE_NAME_ACC_OMP) $(SWM_SOURCES_ACC_OMP)

clean: 
	rm -f $(SWM_EXE_NAME_DPCPP_USM) $(SWM_EXE_NAME_DPCPP_USM_SHARED) $(SWM_EXE_NAME_DPCPP_BUF) $(SWM_EXE_NAME_CPP) $(SWM_EXE_NAME_ACC_OMP) 
