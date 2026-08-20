# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

enable_language(CUDA)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Original CUDA implementation
add_executable(swm_c_cuda
    shallow_swap.cuda.cu
    ../common/wtime.c
)

target_include_directories(swm_c_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../common)

set_target_properties(swm_c_cuda PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

target_compile_options(swm_c_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>)

target_link_libraries(swm_c_cuda PRIVATE m)

# Fused CUDA implementation (all kernels fused into one)
add_executable(swm_c_cuda_fused
    shallow_swap.cuda.fused.cu
    ../common/wtime.c
)

target_include_directories(swm_c_cuda_fused PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../common)

set_target_properties(swm_c_cuda_fused PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

target_compile_options(swm_c_cuda_fused PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>)

target_link_libraries(swm_c_cuda_fused PRIVATE m)
