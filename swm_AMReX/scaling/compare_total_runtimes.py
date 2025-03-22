import numpy as np
import os
import matplotlib.pyplot as plt

## SWM_AMREX_ROOT environment variable must be set
#swm_amrex_root = os.getenv('SWM_AMREX_ROOT')
#if swm_amrex_root is None:
#    raise EnvironmentError("SWM_AMREX_ROOT environment variable is not set")

###############################################################################
# User Input
###############################################################################

#for dir_name in ['strong_scaling_runs', 'strong_scaling_runs_1']:
#  strong_scaling_results_dirs.append(os.path.join(swm_amrex_root, dir_name))

strong_scaling_results_dirs = ['/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_016_chunk_size',
                               '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_032_chunk_size',
                               '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_064_chunk_size',
                               '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_128_chunk_size',
                               '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_256_chunk_size']

strong_scaling_results_dir_2_label = {'/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_016_chunk_size' : 'chunk =  16',
                                      '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_032_chunk_size' : 'chunk =  32',
                                      '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_064_chunk_size' : 'chunk =  64',
                                      '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_128_chunk_size' : 'chunk = 128',
                                      '/home/lalo/SWM/swm_AMReX/from_derecho/strong_scaling_runs_4096_mesh_256_chunk_size' : 'chunk = 256'}

# Controls how the plot looks
label_font_size=18
figure_size=(8,8)

###############################################################################
# Plot total runtime vs N_proc
###############################################################################

runtime_figure = plt.figure(figsize=figure_size)
speedup_figure = plt.figure(figsize=figure_size)
efficiency_figure = plt.figure(figsize=figure_size)

for strong_scaling_results_dir in strong_scaling_results_dirs:

  file_path = os.path.join(strong_scaling_results_dir, 'total_runtime.txt')
  if not os.path.isfile(file_path):
      raise FileNotFoundError(f"File not found: {file_path}")
  
  # Read the CSV file into numpy arrays
  data = np.genfromtxt(file_path, delimiter=',', dtype=None, encoding=None, names=True)
  
  # Extract individual columns into separate numpy arrays
  n_proc_values = data['n_proc']
  run_idx_values = data['run_idx']
  runtime_min_values = data['runtime_min_s']
  runtime_avg_values = data['runtime_avg_s']
  runtime_max_values = data['runtime_max_s']
  
  ###############################################################################
  # Runtime Plot
  ###############################################################################
  
  np_to_plot = []
  t_max_to_plot_best = []
  t_max_to_plot_worst = []
  
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
  
  plt.figure(runtime_figure.number)
  #base_filename = os.path.splitext(os.path.basename(file_path))[0]
  #plt.plot(np_to_plot, t_max_to_plot_best, label = os.path.basename(strong_scaling_results_dir))
  plt.plot(np_to_plot, t_max_to_plot_best, label = strong_scaling_results_dir_2_label[strong_scaling_results_dir])

  #plt.fill_between(np_to_plot, t_max_to_plot_best, t_max_to_plot_worst, alpha=0.2)

  plt.figure(speedup_figure.number)
  t_serial = t_max_to_plot_best[0]
  plt.plot(np_to_plot, t_serial/t_max_to_plot_best, label = strong_scaling_results_dir_2_label[strong_scaling_results_dir])

  plt.figure(efficiency_figure.number)
  t_serial = t_max_to_plot_best[0]
  plt.plot(np_to_plot, t_serial/(np_to_plot*np.array(t_max_to_plot_best)), label = strong_scaling_results_dir_2_label[strong_scaling_results_dir])


plt.figure(runtime_figure.number)
plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Runtime [s]',fontsize=label_font_size)
plt.legend()
plt.tight_layout()


plt.figure(speedup_figure.number)
plt.plot([1,np.max(np_to_plot)], [1,np.max(np_to_plot)], '--k', label='Ideal') # Ideal speedup for reference
plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Speedup $= \frac{t_{\mathrm{serial}}}{t_\mathrm{parallel}}$',fontsize=label_font_size)
plt.legend()
plt.gca().set_aspect('equal')
plt.tight_layout()


plt.figure(efficiency_figure.number)
plt.plot([1,np.max(np_to_plot)], [1,1], '--k', label='Ideal')
plt.xlabel(r'$N_{\mathrm{proc}}$',fontsize=label_font_size)
plt.ylabel(r'Efficiency $= \frac{\mathrm{Speedup}}{N_{\mathrm{proc}}} =  \frac{t_{\mathrm{serial}}}{N_{\mathrm{proc}}t_{\mathrm{parallel}}}$',fontsize=label_font_size)
plt.legend()
plt.tight_layout()

plt.show()