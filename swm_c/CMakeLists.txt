
cmake_minimum_required(VERSION 3.0)


# Add the source files
set(SOURCES
    shallow_swap.c
    wtime.c
)

# Generate the executable
add_executable(swm.exe ${SOURCES})

target_link_libraries(swm.exe m)