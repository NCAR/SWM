set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Add optimization flag
#add_compile_options(-O2)
# Maybe move this to be set if CMAKE_BUILD_TYPE is Release
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")

# Flags that get used when building  debug flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g")

add_subdirectory(c)

if (SWM_OPENACC)
  add_subdirectory(c_OpenACC)
endif()