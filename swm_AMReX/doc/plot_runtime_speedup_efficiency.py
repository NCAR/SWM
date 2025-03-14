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
# No MPI optimizations
derecho_MPI_old = {
    'n_procs': np.array([1, 2, 4, 8, 16]),
    't': np.array([25.01, 16.52, 16.42, 16.39, 15.44]),
    'legend_label': 'Derecho MPI - Original',
    'line_type': '-o'
}

# Deriechio MPI... 4096x4096, 1024 max chunk size, 100 time steps 
# Using swap instead of copy. Made sure that partitioning and mapping look ok. Removed several uneeded fill boundary calls for multfabs.
derecho_MPI_new = {
    'n_procs': np.array([1, 2, 4, 8, 16]),
    't': np.array([21.94, 13.92, 13.65, 13.83, 13.22]),
    'legend_label': 'Derecho MPI - Swap + Less fill boundary',
    'line_type': '-o'
}

cases_to_plot = [derecho_MPI_old, derecho_MPI_new]

###############################################################################
# Runtime Plot
###############################################################################
plt.figure(figsize=figure_size)

for case in cases_to_plot:
  plt.plot(case['n_procs'], case['t'], case['line_type'], label=case['legend_label'])

plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
plt.legend()
plt.tight_layout()


###############################################################################
# Speedup Plot
###############################################################################
plt.figure(figsize=figure_size)

# Ideal speedup for reference
n_procs_max = max([np.max(case['n_procs']) for case in cases_to_plot])
plt.plot([1,n_procs_max], [1,n_procs_max], '--k', label='Ideal')

# Actual speedup
for case in cases_to_plot:
  assert(case['n_procs'][0] == 1) # Serial time must be the first element or speedup calculation y label will be misleading
  t_serial = case['t'][0]
  plt.plot(case['n_procs'], t_serial/case['t'], case['line_type'], label=case['legend_label'])

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
plt.plot([1,n_procs_max], [1,1], '--k', label='Ideal')

# Actual efficiency
for case in cases_to_plot:
  assert(case['n_procs'][0] == 1) # Serial time must be the first element or speedup calculation y label will be misleading
  t_serial = case['t'][0]
  plt.plot(case['n_procs'], t_serial/(case['n_procs']*case['t']), case['line_type'], label=case['legend_label'])

plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Efficiency $= \frac{\mathrm{Speedup}}{N_{\mathrm{proc}}} =  \frac{t_{\mathrm{serial}}}{N_{\mathrm{proc}}t_{\mathrm{parallel}}}$',fontsize=label_font_size)
plt.legend()
plt.tight_layout()

plt.show()