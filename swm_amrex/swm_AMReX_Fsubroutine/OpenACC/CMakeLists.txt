################################################################################
## Hack to run old Makefile
################################################################################
#add_custom_target(swm_amrex_fsubroutine_acc ALL
#    COMMAND ${CMAKE_MAKE_PROGRAM} -f GNUmakefile
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    COMMENT "Building with legacy GNUmakefile"
#)

###############################################################################
# Using CMake to buidd the SWM AMReX Fsubroutine with OpenACC support
###############################################################################
find_package(AMReX REQUIRED FORTRAN FINTERFACES) 

# TODO: This is basicaly a copy of the block from the CMakeLists.txt for the C OpenACC build. Maybe combine in a function?
# Set OpenACC flags based on compiler and device
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_COMPILE_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_LINK_FLAGS -acc=multicore)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma)
        #set(OPENACC_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma -mcmodel=medium)
        set(OPENACC_LINK_FLAGS -acc=gpu)
    endif()
else()
    # Maybe add a warning here that we are not using OpenACC for whatever compiler you are using.
    message(FATAL_ERROR "We have not added support for your compiler (${CMAKE_Fortran_COMPILER_ID}) for OpenACC build. Exiting configuration.")
endif()


set(SWM_AMREX_FSUBROUTINE_ACC_SOURCES
    ${PROJECT_SOURCE_DIR}/swm_amrex/common/main.cpp
    ${PROJECT_SOURCE_DIR}/swm_amrex/common/swm_mini_app_utils.cpp
    ../swm_mini_app_utils_fsubroutine.cpp
    swm_fortran_kernelsACC.F90
)

add_executable(swm_amrex_fsubroutine_acc
    ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES}
)

#set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_OPTIONS "${OPENACC_COMPILE_FLAGS}")
set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_OPTIONS "${OPENACC_COMPILE_FLAGS};-Mnomain")
#set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_OPTIONS "${OPENACC_COMPILE_FLAGS};-Mnomain;-gpu=cc80,fastmath;-cuda;-gpu=lineinfo;-gpu=keepptx;-gpu=maxregcount:255;-Mdclchk;-gpu=rdc")
#set_source_files_properties(swm_fortran_kernelsACC.F90 PROPERTIES COMPILE_OPTIONS "-O2;-fast;-Mnomain;-gpu=cc80,fastmath;-cuda;-gpu=lineinfo;-gpu=keepptx;-gpu=maxregcount:255;-Mdclchk;-acc=gpu;-Minfo=accel;-mcmodel=medium;-gpu=cc80")

#target_link_options(swm_amrex_fsubroutine_acc PRIVATE ${OPENACC_LINK_FLAGS})
target_link_options(swm_amrex_fsubroutine_acc PRIVATE "${OPENACC_LINK_FLAGS};-Mnomain")
#target_link_options(swm_amrex_fsubroutine_acc PRIVATE "${OPENACC_LINK_FLAGS};-cuda;-Mnomain;-gpu=maxregcount:255;-Mdclchk;-gpu=cc80")
#target_link_options(swm_amrex_fsubroutine_acc PRIVATE "-O2;-fast;-Mnomain;-gpu=cc80,fastmath;-cuda;-gpu=lineinfo;-gpu=keepptx;-gpu=maxregcount:255;-Mdclchk;-acc=gpu;-Minfo=accel;-mcmodel=medium;-gpu=cc80")

#target_link_options(swm_amrex_fsubroutine_acc
#  PRIVATE
#    $<$<OR:$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,Fortran>,$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,CUDA>>:${OPENACC_LINK_FLAGS}>
#)
#if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
#  target_link_options(swm_amrex_fsubroutine_acc PRIVATE ${OPENACC_LINK_FLAGS})
#endif()

if (SWM_CUDA)
  include("AMReXTargetHelpers")
  setup_target_for_cuda_compilation(swm_amrex_fsubroutine_acc
    SOURCES ${SWM_AMREX_FSUBROUTINE_ACC_SOURCES}
  )
endif()

#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gpu=rdc")

target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${PROJECT_SOURCE_DIR}/swm_amrex/common)
target_include_directories(swm_amrex_fsubroutine_acc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC AMReX::amrex_2d)
target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC m)

#target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC cuda curand cusparse)
#target_link_libraries(swm_amrex_fsubroutine_acc PUBLIC atomic nvf rt)