import numpy as np
import os
import matplotlib.pyplot as plt

# SWM_AMREX_ROOT environment variable must be set
swm_amrex_root = os.getenv('SWM_AMREX_ROOT')
if swm_amrex_root is None:
    raise EnvironmentError("SWM_AMREX_ROOT environment variable is not set")


###############################################################################
# User Input
###############################################################################

strong_scaling_results_dir = os.path.join(swm_amrex_root, 'strong_scaling_runs')
if not os.path.isdir(strong_scaling_results_dir):
    raise NotADirectoryError(f"Directory not found: {strong_scaling_results_dir}")

# Controls how the plot looks
label_font_size=18
figure_size=(8,8)


###############################################################################
# Plot total runtime vs N_proc
###############################################################################

file_path = os.path.join(strong_scaling_results_dir, 'total_runtime.txt')
if not os.path.isfile(file_path):
    raise FileNotFoundError(f"File not found: {file_path}")


# Read the CSV file into numpy arrays
data = np.genfromtxt(file_path, delimiter=',', dtype=None, encoding=None, names=True)

#for key in data.dtype.names:
#    print(key)
#    print(data[key])

# Extract individual columns into separate numpy arrays
n_proc_values = data['n_proc']
run_idx_values = data['run_idx']
runtime_min_values = data['runtime_min_s']
runtime_avg_values = data['runtime_avg_s']
runtime_max_values = data['runtime_max_s']

# Print the arrays to verify
#print("n_proc:", n_proc_values)
#print("run_idx:", run_idx_values)
#print("runtime_min:", runtime_min_values)
#print("runtime_avg:", runtime_avg_values)
#print("runtime_max:", runtime_max_values)

###############################################################################
# Runtime Plot
###############################################################################
plt.figure(figsize=figure_size)

np_to_plot = []
t_max_to_plot_best = []
t_max_to_plot_worst = []

np_to_best_run_idx  = {}

# Could use this to get the most unloadbalanced run... would find the largest spread between min and max for a given run per each N_proc
#t_min_to_plot = []

for n_proc in np.unique(n_proc_values):

    # Get the data for this specific value of N_proc
    np_mask = (n_proc_values == n_proc)
    run_idx_values_for_np = run_idx_values[np_mask]
    runtime_min_values_for_np = runtime_min_values[np_mask]
    runtime_avg_values_for_np = runtime_avg_values[np_mask]
    runtime_max_values_for_np = runtime_max_values[np_mask]

    # Will be used as x axis
    np_to_plot.append(n_proc)

    # Get the fastest and slowest runtimes for this N_proc
    # Will be used as y axis
    t_max_best = min(runtime_max_values_for_np)
    t_max_worst = max(runtime_max_values_for_np)

    t_max_to_plot_best.append(t_max_best)
    t_max_to_plot_worst.append(t_max_worst)

    # keep track of run_idx for the best run
    np_to_best_run_idx[n_proc] = np.argmin(runtime_max_values_for_np)


print(np_to_best_run_idx)

base_filename = os.path.splitext(os.path.basename(file_path))[0]
plt.plot(np_to_plot, t_max_to_plot_best, label = base_filename)

#plt.fill_between(np_to_plot, t_max_to_plot_best, t_max_to_plot_worst, alpha=0.2)

plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
#plt.legend()
plt.tight_layout()

###############################################################################
# Plot runtime vs N_proc for specific timers
###############################################################################


file_basenames = ['exclusive_runtime_UpdateIntermediateVariables',
                   'exclusive_runtime_UpdateNewVariables',
                   'exclusive_runtime_UpdateOldVariables',
                   'inclusive_runtime_FabArray::FillBoundary'
                   ]
files = []
for file_basename in file_basenames:
    files.append(os.path.join(strong_scaling_results_dir, file_basename+'.txt'))

for file_path in files:

  if not os.path.isfile(file_path):
      raise FileNotFoundError(f"File not found: {file_path}")

  data = np.genfromtxt(file_path, delimiter=',', dtype=None, encoding=None, names=True)
  
  # Extract individual columns into separate numpy arrays
  n_proc_values = data['n_proc']
  run_idx_values = data['run_idx']
  runtime_min_values = data['runtime_min_s']
  runtime_avg_values = data['runtime_avg_s']
  runtime_max_values = data['runtime_max_s']
  max_percent_values = data['max_percent']
  
  np_to_plot = []
  t_max_to_plot = []
  
  for n_proc in np.unique(n_proc_values):
  
    # Get the runtimes for this N_proc
    np_mask = (n_proc_values == n_proc)
    runtime_min_values_for_np = runtime_min_values[np_mask]
    runtime_avg_values_for_np = runtime_avg_values[np_mask]
    runtime_max_values_for_np = runtime_max_values[np_mask]
    max_percent_values_for_np = max_percent_values[np_mask]
  
    # Will be used as x axis
    np_to_plot.append(n_proc)
  
    # Will be used as y axis
    idx = np_to_best_run_idx[n_proc]
    t_max_to_plot.append(runtime_max_values_for_np[idx])
  
  base_filename = os.path.splitext(os.path.basename(file_path))[0]
  plt.plot(np_to_plot, t_max_to_plot, label = base_filename)
  
plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
plt.legend()
plt.tight_layout()

plt.show()