set(SWM_KOKKOS_FORTRAN_OPENACC_TARGET shallow_mdrangepolicy_fortran_acc)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_FORTRAN_COMPILE_FLAGS -fopenacc -foffload=disable)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_FORTRAN_COMPILE_FLAGS -fopenacc)
    endif()
    set(OPENACC_FORTRAN_LINK_FLAGS -fopenacc)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    if(SWM_DEVICE STREQUAL "cpu")
        set(OPENACC_FORTRAN_COMPILE_FLAGS -acc=multicore -Minfo=accel -Mnofma)
        set(OPENACC_FORTRAN_LINK_FLAGS -acc=multicore)
    elseif(SWM_DEVICE STREQUAL "gpu")
        set(OPENACC_FORTRAN_COMPILE_FLAGS -acc=gpu -Minfo=accel -Mnofma)
        set(OPENACC_FORTRAN_LINK_FLAGS -acc=gpu)
    endif()
else()
    message(WARNING "OpenACC Fortran flags are not configured for compiler ${CMAKE_Fortran_COMPILER_ID}; directives may be ignored.")
    set(OPENACC_FORTRAN_COMPILE_FLAGS "")
    set(OPENACC_FORTRAN_LINK_FLAGS "")
endif()

add_executable(${SWM_KOKKOS_FORTRAN_OPENACC_TARGET}
    shallow_mdrangepolicy.cpp
    shallow_fortran_acc.F90
)

target_link_libraries(${SWM_KOKKOS_FORTRAN_OPENACC_TARGET} Kokkos::kokkos)

foreach(flag IN LISTS OPENACC_FORTRAN_COMPILE_FLAGS)
    target_compile_options(${SWM_KOKKOS_FORTRAN_OPENACC_TARGET}
        PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${flag}>"
    )
endforeach()

target_link_options(${SWM_KOKKOS_FORTRAN_OPENACC_TARGET}
    PRIVATE ${OPENACC_FORTRAN_LINK_FLAGS}
)
