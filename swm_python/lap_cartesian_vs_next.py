#!/usr/bin/env python
# coding: utf-8

# GT4Py - GridTools for Python
# 
# Copyright (c) 2014-2023, ETH Zurich
# All rights reserved.
# 
# This file is part the GT4Py project and the GridTools framework.
# GT4Py is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or any later
# version. See the LICENSE.txt file at the top-level directory of this
# distribution for a copy of the license or check <https://www.gnu.org/licenses/>.
# 
# SPDX-License-Identifier: GPL-3.0-or-later

# # Demonstrates gt4py.cartesian with gt4py.next compatibility

# Imports

# In[19]:


import numpy as np

nx = 32
ny = 32
nz = 1
dtype = np.float64


# Storages
# --
# 
# We create fields using the gt4py.next constructors. These fields are compatible with gt4py.cartesian when we use "I", "J", "K" as the dimension names.

# In[8]:


import gt4py.next as gtx

#allocator = gtx.itir_python # should match the executor
#allocator = gtx.gtfn_cpu
allocator = gtx.gtfn_gpu

# Note: for gt4py.next, names don't matter, for gt4py.cartesian they have to be "I", "J", "K"
I = gtx.Dimension("I")
J = gtx.Dimension("J")
K = gtx.Dimension("K", kind=gtx.DimensionKind.VERTICAL)

domain = gtx.domain({I: nx, J: ny, K: nz})

inp = gtx.as_field(domain, np.fromfunction(lambda x, y, z: x**2+y**2, shape=(nx, ny, nz)), dtype, allocator=allocator)
out_cartesian = gtx.zeros(domain, dtype, allocator=allocator)
out_next = gtx.zeros(domain, dtype, allocator=allocator)


# In[25]:


#get_ipython().system('module load gcc')
#get_ipython().system('export BOOST_HOME=/glade/derecho/scratch/haiyingx/boost_1_84_0/include/boost')


# gt4py.cartesian
# --

# In[26]:


import gt4py.cartesian.gtscript as gtscript

#cartesian_backend = "numpy"
#cartesian_backend = "gt:cpu_ifirst"
cartesian_backend = "gt:gpu"

@gtscript.stencil(backend=cartesian_backend)
def lap_cartesian(
    inp: gtscript.Field[dtype],
    out: gtscript.Field[dtype],
):
    with computation(PARALLEL), interval(...):
        out = -4.0 * inp[0, 0, 0] + inp[-1, 0, 0] + inp[1, 0, 0] + inp[0, -1, 0] + inp[0, 1, 0]

#lap_cartesian(inp=inp, out=out_cartesian, origin=(1, 1, 0), domain=(nx-2, ny-2, nz))


# In[ ]:


from gt4py.next import Field

#next_backend = gtx.itir_python
# next_backend = gtx.gtfn_cpu
next_backend = gtx.gtfn_gpu

Ioff = gtx.FieldOffset("I", source=I, target=(I,))
Joff = gtx.FieldOffset("J", source=J, target=(J,))

@gtx.field_operator
def lap_next(inp: Field[[I, J, K], dtype]) -> Field[[I, J, K], dtype]:
    return -4.0 * inp + inp(Ioff[-1]) + inp(Ioff[1]) + inp(Joff[-1]) + inp(Joff[1])

@gtx.program(backend=next_backend)
def lap_next_program(inp: Field[[I, J, K], dtype], out: Field[[I, J, K], dtype]):
    lap_next(inp, out=out[1:-1, 1:-1, :])

lap_next_program(inp, out_next, offset_provider={"Ioff": I, "Joff": J})


# In[ ]:


assert np.allclose(out_cartesian.asnumpy(), out_next.asnumpy())


# In[ ]:




