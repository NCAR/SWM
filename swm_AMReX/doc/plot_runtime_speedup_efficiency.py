import numpy as np
import matplotlib.pyplot as plt

label_font_size=18
figure_size=(8,8)

# Fake data to test the plot
#n_procs = np.array([1,2,4,8]);
#t = np.array([1,0.5,0.28,0.2]);
#legend_label='test data'
#line_type='-*'

## Laptop MPI... 4096x4096, 1024 max chunk size, 100 time steps 
#n_procs = np.array([1, 2, 4, 8, 16])
#t = np.array([17.62, 11.69, 9.27, 8.843, 9.47])
#legend_label='Local Machine MPI'
#line_type='-*'

# Deriechio MPI... 4096x4096, 1024 max chunk size, 100 time steps 
n_procs = np.array([1, 2, 4, 8, 16])
t = np.array([25.01, 16.52, 16.42, 16.39, 15.44])
legend_label='Derecho MPI'
line_type='-o'

###############################################################################
# Runtime Plot
###############################################################################
plt.figure(figsize=figure_size)

plt.plot(n_procs, t, line_type, label=legend_label)

plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
plt.legend()
plt.tight_layout()

###############################################################################
# Speedup Plot
###############################################################################
plt.figure(figsize=figure_size)

# Ideal speedup for reference
plt.plot(n_procs, n_procs, '--k', label='Ideal')

# Actual speedup
assert(n_procs[0] == 1) # Serial time must be the first element or speedup calculation y label will be misleading
t_serial = t[0]
plt.plot(n_procs, t_serial/t, line_type, label=legend_label)

plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Speedup $= \frac{t_{\mathrm{serial}}}{t_\mathrm{parallel}}$',fontsize=label_font_size)
plt.legend()
plt.gca().set_aspect('equal')
plt.tight_layout()

###############################################################################
# Efficiency Plot
###############################################################################
plt.figure(figsize=figure_size)

# Ideal efficiency for reference
plt.plot(n_procs, 1.0*np.ones(len(n_procs)), '--k', label='Ideal')

# Actual efficiency
assert(n_procs[0] == 1) # Serial time must be the first element or speedup calculation y label will be misleading
t_serial = t[0]
plt.plot(n_procs, t_serial/(n_procs*t), line_type, label=legend_label)

plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Efficiency $= \frac{\mathrm{Speedup}}{N_{\mathrm{proc}}} =  \frac{t_{\mathrm{serial}}}{N_{\mathrm{proc}}t_{\mathrm{parallel}}}$',fontsize=label_font_size)
plt.legend()
plt.tight_layout()

plt.show()