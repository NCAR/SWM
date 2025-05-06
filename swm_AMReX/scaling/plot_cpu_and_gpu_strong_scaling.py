import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

##############################################################################
# User Input
##############################################################################

# Controls which plots to generate
plot_cpu_core_single_node = False
plot_cpu_core = False
plot_cpu_node = False
plot_gpu = False
plot_compare_cpu_node_gpu = True

# Set default figure size in inches
plt.rcParams["figure.figsize"] = (8, 8) 

# Font size for x and y axis labels
mpl.rcParams['axes.labelsize'] = 24

# Font size for the x and y axis tick labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

# Font size for lebels in the legend
plt.rc('legend', fontsize=20) 

##############################################################################
# Strong Scaling Data to Plot
##############################################################################

# CPU - Total Runtime vs N CPU Core

# dictionary to hold the number of CPU cores and the corresponding runtimes for each mesh size
n_cpu_core = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]

t_cpu_core = {}
t_cpu_core[1024] = [1.341, 0.6346, 0.2891, 0.1434, 0.08059, 0.05177, 0.04403, 0.045, 0.03228, 0.02579, 0.02533, 0.02398, 0.02411]
t_cpu_core[2048] = [4.642, 2.285, 1.114, 0.596, 0.2524, 0.1446, 0.0987, 0.09266, 0.06734, 0.05125, 0.03412, 0.02748, 0.2614]
t_cpu_core[4096] = [17.38, 8.72, 4.347, 2.15, 1.41, 1.286, 1.307, 1.352, 0.4695, 0.1029, 0.07226, 0.05467, 0.03503]
t_cpu_core[8192] = [79.53, 39.91, 20.09, 10.12, 6.709, 6.196, 6.336, 6.432, 3.071, 1.353, 0.4719, 0.1066, 0.0781]
t_cpu_core[16384] = [335.6, 169.5, 89.94, 54.56, 27.09, 24.99, 25.62, 26.11, 13.05, 6.462, 3.083, 1.358, 0.4745]

# CPU - Total Runtime vs N CPU Node (Ran on Derecho so each node has 128 cores)... the times are the same as cpu_core data above but just from 128 cores to the end.
# dictionary to hold the number of CPU cores and the corresponding runtimes for each mesh size
n_cpu_node = [1, 2, 4, 8, 16, 32]

t_cpu_node = {}
for mesh_size in t_cpu_core:
    t_cpu_node[mesh_size] = t_cpu_core[mesh_size][7:]

# GPU - Total Runtime vs N GPU
n_gpu = [1, 2, 4, 8, 16, 32, 64]

t_gpu = {}
t_gpu[4096] =  [0.3004, 0.1902, 0.1207, 0.09364, 0.06838, 0.06885, 0.05205]
t_gpu[8192] =  [1.166, 0.6295, 0.35, 0.2103, 0.1281, 0.09877, 0.06967]
t_gpu[16384] = [4.633, 2.381, 1.254, 0.6688, 0.3576, 0.2139, 0.1288]


##############################################################################
# Helper Functions
##############################################################################
def speedup(t_serial, t_parallel):
    return t_serial/np.array(t_parallel)

def efficiency(t_serial, t_parallel, n_proc):
    return speedup(t_serial, t_parallel)/np.array(n_proc)

def create_runtime_plot(x_ticks, x_label_text):
    figure = plt.figure()

    # Set the x and y labels
    plt.xlabel(x_label_text)
    plt.ylabel(r'Runtime [s]')

    # Set x-axis ticks and labels to match the plotted points
    ax = figure.gca()  
    ax.set_xticks(x_ticks)  
    ax.set_xticklabels(x_ticks)  

    plt.tight_layout()

    return figure

def create_speedup_plot(x_ticks, x_label_text):
    figure = plt.figure()

    # Plot ideal speedup for reference. Should be a 45 degree line up and to the right.
    plt.plot([1,np.max(x_ticks)], [1,np.max(x_ticks)], '--k', label='Ideal') 

    # Set the x and y labels
    plt.xlabel(x_label_text)
    plt.ylabel(r'Speedup $= \frac{t_{\mathrm{serial}}}{t_\mathrm{parallel}}$')

    # Set x-axis ticks and labels to match the plotted points
    ax = figure.gca()  
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks)  

    # Set aspect ratio to be equal to make the plot square. Makes it so the ideal speedup line is a 45 degree line.
    ax.set_aspect('equal') 

    plt.legend()
    plt.tight_layout()

    return figure

def create_efficiency_plot(x_ticks, x_label_text):
    figure = plt.figure()

    # Plot ideal efficiency for reference. Should be a horizontal line at y=1.
    plt.plot([1,np.max(x_ticks)], [1,1], '--k', label='Ideal')

    # Set the x and y labels
    plt.xlabel(x_label_text)
    plt.ylabel(r'Efficiency $= \frac{\mathrm{Speedup}}{N_{\mathrm{proc}}} =  \frac{t_{\mathrm{serial}}}{N_{\mathrm{proc}}t_{\mathrm{parallel}}}$')

    # Set x-axis ticks and labels to match the plotted points
    ax1 = figure.gca()  
    ax1.set_xticks(x_ticks)  
    ax1.set_xticklabels(x_ticks)  

    plt.legend()
    plt.tight_layout()

    return figure


# Useful labels for legends
legend_label_cpu_core = 'CPU Cores'
legend_label_cpu_node = 'CPU Nodes'
legend_label_gpu = 'GPU'

legend_label_mesh = dict()
legend_label_mesh[1024] =   '1024 x 1024'
legend_label_mesh[2048] =   '2048 x 2048'
legend_label_mesh[4096] =   '4096 x 4096'
legend_label_mesh[8192] =   '8192 x 8192'
legend_label_mesh[16384] = '16384 x 16384'

# Useful labels for the x-axis
x_label_cpu_core = r'$N_{\mathrm{CPU \ Core}}$'
x_label_cpu_node = r'$N_{\mathrm{CPU \ NODE}}$'
x_label_gpu = r'$N_{\mathrm{GPU}}$'

# Useful for line type, color, and marker types
line_type = dict()
line_type[16384] = '-ro'
line_type[8192] = '-.bs'
line_type[4096] = ':gd'
line_type[2048] = '-kx'
line_type[1024] = '-.c^'

line_type['cpu'] = '-.ro'
line_type['gpu'] = '-gs'

##############################################################################
# Runtime Figures
##############################################################################
if plot_cpu_core_single_node:
  runtime_figure_cpu_core_single_node = create_runtime_plot(n_cpu_core[0:8], x_label_cpu_core)
  for mesh_size in reversed(t_cpu_core):
    plt.plot(n_cpu_core[0:8], t_cpu_core[mesh_size][0:8], line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_core_single_node.png", dpi=300)

if plot_cpu_core:
  runtime_figure_cpu_core = create_runtime_plot(n_cpu_core, x_label_cpu_core)
  for mesh_size in reversed(t_cpu_core):
     plt.plot(n_cpu_core, t_cpu_core[mesh_size], line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_core.png", dpi=300)

if plot_cpu_node:
  runtime_figure_cpu_node = create_runtime_plot(n_cpu_node, x_label_cpu_node)
  for mesh_size in reversed(t_cpu_node):
     plt.plot(n_cpu_node, t_cpu_node[mesh_size], line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_node.png", dpi=300)

if plot_gpu:
  runtime_figure_gpu = create_runtime_plot(n_gpu, x_label_gpu)
  for mesh_size in reversed(t_gpu):
     plt.plot(n_gpu, t_gpu[mesh_size], line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_runtime_gpu.png", dpi=300)

if plot_compare_cpu_node_gpu:
  # CPU NODE and GPU Runtime Comparison
  for mesh_size in reversed(t_gpu):
    runtime_figure_cpu_node = create_runtime_plot(n_cpu_node, x_label_cpu_node)
    #plt.plot(n_cpu_node, t_cpu_node[mesh_size], line_type['cpu'], label='CPU Node - '+legend_label_mesh[mesh_size])
    #plt.plot(n_gpu[:-1], t_gpu[mesh_size][:-1], line_type['gpu'], label='GPU Node - '+legend_label_mesh[mesh_size])
    plt.semilogy(n_cpu_node, t_cpu_node[mesh_size], line_type['cpu'], label='CPU Node - '+legend_label_mesh[mesh_size])
    plt.semilogy(n_gpu[:-1], t_gpu[mesh_size][:-1], line_type['gpu'], label='GPU Node - '+legend_label_mesh[mesh_size])
    plt.xlabel('N')
    plt.legend()
    plt.savefig(f"strong_scaling_runtime_cpu_node_and_gpu_{mesh_size}.png", dpi=300)
  

##############################################################################
# Speedup Figures
##############################################################################
if plot_cpu_core_single_node:
  speedup_figure_cpu_core_single_node  = create_speedup_plot(n_cpu_core[0:8], x_label_cpu_core)
  for mesh_size in reversed(t_cpu_core):
     plt.plot(n_cpu_core[0:8], speedup(t_cpu_core[mesh_size][0], t_cpu_core[mesh_size][0:8]), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_core_single_node.png", dpi=300)

if plot_cpu_core:
  speedup_figure_cpu_core  = create_speedup_plot(n_cpu_core, x_label_cpu_core)
  for mesh_size in reversed(t_cpu_core):
     plt.plot(n_cpu_core, speedup(t_cpu_core[mesh_size][0], t_cpu_core[mesh_size]), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_core.png", dpi=300)

if plot_cpu_node:
  speedup_figure_cpu_node = create_speedup_plot(n_cpu_node, x_label_cpu_node)
  for mesh_size in reversed(t_cpu_node):
     plt.plot(n_cpu_node, speedup(t_cpu_node[mesh_size][0], t_cpu_node[mesh_size]), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_node.png", dpi=300)

if plot_gpu:
  speedup_figure_gpu = create_speedup_plot(n_gpu, x_label_gpu)
  for mesh_size in reversed(t_gpu):
     plt.plot(n_gpu, speedup(t_gpu[mesh_size][0], t_gpu[mesh_size]), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_speedup_gpu.png", dpi=300)

#if plot_compare_cpu_node_gpu:
#  # CPU NODE and GPU Speedup Comparison
#  for mesh_size in reversed(t_gpu):
#    speedup_figure_cpu_node_and_gpu = create_speedup_plot(n_cpu_node, x_label_cpu_node)
#    plt.plot(n_cpu_node, speedup(t_cpu_node[mesh_size][0], t_cpu_node[mesh_size]), line_type['cpu'], label='CPU Node - '+legend_label_mesh[mesh_size])
#    plt.plot(n_gpu[:-1], speedup(t_gpu[mesh_size][0], t_gpu[mesh_size][:-1]), line_type['gpu'], label='GPU - '+legend_label_mesh[mesh_size])
#    plt.xlabel('N')
#    plt.legend()
#    plt.savefig(f"strong_scaling_speedup_cpu_node_and_gpu_{mesh_size}.png", dpi=300)

##############################################################################
# Efficiency Figures
##############################################################################
if plot_cpu_core_single_node:
  efficiency_figure_cpu_core_single_node = create_efficiency_plot(n_cpu_core[0:8], x_label_cpu_core)
  for mesh_size in reversed(t_cpu_core):
      plt.plot(n_cpu_core[0:8], efficiency(t_cpu_core[mesh_size][0], t_cpu_core[mesh_size][0:8], n_cpu_core[0:8]), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_core_single_node.png", dpi=300)

if plot_cpu_core:
  efficiency_figure_cpu_core = create_efficiency_plot(n_cpu_core, x_label_cpu_core)
  for mesh_size in reversed(t_cpu_core):
      plt.plot(n_cpu_core, efficiency(t_cpu_core[mesh_size][0], t_cpu_core[mesh_size], n_cpu_core), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_core.png", dpi=300)

if plot_cpu_node:
  efficiency_figure_cpu_node = create_efficiency_plot(n_cpu_node, x_label_cpu_node)
  for mesh_size in reversed(t_cpu_node):
      plt.plot(n_cpu_node, efficiency(t_cpu_node[mesh_size][0], t_cpu_node[mesh_size], n_cpu_node), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_node.png", dpi=300)

if plot_gpu:
  efficiency_figure_gpu = create_efficiency_plot(n_gpu, x_label_gpu)
  for mesh_size in reversed(t_gpu):
      plt.plot(n_gpu, efficiency(t_gpu[mesh_size][0], t_gpu[mesh_size], n_gpu), line_type[mesh_size], label=legend_label_mesh[mesh_size])
  plt.legend()
  plt.savefig("strong_scaling_efficiency_gpu.png", dpi=300)

#if plot_compare_cpu_node_gpu:
#  # CPU NODE and GPU Efficiency Comparison
#  for mesh_size in reversed(t_gpu):
#    efficiency_figure_cpu_node_and_gpu = create_efficiency_plot(n_cpu_node, x_label_cpu_node)
#    plt.plot(n_cpu_node, efficiency(t_cpu_node[mesh_size][0], t_cpu_node[mesh_size], n_cpu_node), line_type['cpu'], label='CPU Node - '+legend_label_mesh[mesh_size])
#    plt.plot(n_gpu[:-1], efficiency(t_gpu[mesh_size][0], t_gpu[mesh_size][:-1], n_gpu[:-1]), line_type['gpu'], label='GPU - '+legend_label_mesh[mesh_size])
#    plt.xlabel('N')
#    plt.legend()
#    plt.savefig(f"strong_scaling_efficiency_cpu_node_and_gpu_{mesh_size}.png", dpi=300)

plt.show()