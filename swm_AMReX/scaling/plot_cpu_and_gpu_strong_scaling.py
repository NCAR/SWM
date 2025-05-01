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
plot_gpu = True
plot_compare_cpu_node_gpu = False

# Controls how the plots looks

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

# Note: 
#     small  mesh is   4096 x  4096
#     medium mesh is   8192 x  8192
#     large  mesh is  16384 x 16384

# CPU - Total Runtime vs N CPU Core
n_cpu_core_small_mesh = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
t_cpu_core_small_mesh = np.array([17.38, 8.72, 4.347, 2.15, 1.41, 1.286, 1.307, 1.352, 0.4695, 0.1029, 0.07226, 0.05467, 0.03503])

n_cpu_core_medium_mesh = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
t_cpu_core_medium_mesh = np.array([79.53, 39.91, 20.09, 10.12, 6.709, 6.196, 6.336, 6.432, 3.071, 1.353, 0.4719, 0.1066, 0.0781])

n_cpu_core_large_mesh = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
t_cpu_core_large_mesh = np.array([335.6, 169.5, 89.94, 54.56, 27.09, 24.99, 25.62, 26.11, 13.05, 6.462, 3.083, 1.358, 0.4745])


# CPU - Total Runtime vs N CPU Node (Ran on Derecho so each node has 128 cores)... the times are the same as cpu_core data above but just from 128 cores to the end.
n_cpu_node_small_mesh = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node_small_mesh = np.array([1.352, 0.4695, 0.1029, 0.07226, 0.05467, 0.03503])

n_cpu_node_medium_mesh = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node_medium_mesh = np.array([6.558, 3.102, 1.382, 0.4829, 0.1214, 0.08429])

n_cpu_node_large_mesh = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node_large_mesh = np.array([26.11, 13.05, 6.462, 3.083, 1.358, 0.4745])


# GPU - Total Runtime vs N GPU
n_gpu_small_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_small_mesh = np.array([0.3004, 0.1902, 0.1207, 0.09364, 0.06838, 0.06885, 0.05205])

n_gpu_medium_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_medium_mesh = np.array([1.166, 0.6295, 0.35, 0.2103, 0.1281, 0.09877, 0.06967])

n_gpu_large_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_large_mesh = np.array([4.633, 2.381, 1.254, 0.6688, 0.3576, 0.2139, 0.1288])


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
legend_label_large_mesh  = '16384 x 16384'
legend_label_medium_mesh =   '8192 x 8192'
legend_label_small_mesh  =   '4096 x 4096'

# Useful labels for the x-axis
x_label_cpu_core = r'$N_{\mathrm{CPU \ Core}}$'
x_label_cpu_node = r'$N_{\mathrm{CPU \ NODE}}$'
x_label_gpu = r'$N_{\mathrm{GPU}}$'

# Useful for line type, color, and marker types
small_mesh_line_type = ':'
medium_mesh_line_type = '-.'
large_mesh_line_type = '-'

color_1 = 'r'
color_2 = 'b'
color_3 = 'g'

marker_type_1 = 'o'
marker_type_2 = 's'
marker_type_3 = 'd'

##############################################################################
# Runtime Figures
##############################################################################
if plot_cpu_core_single_node:
  runtime_figure_cpu_core_single_node = create_runtime_plot(n_cpu_core_large_mesh[0:8], x_label_cpu_core)
  plt.plot(n_cpu_core_large_mesh[0:8],  t_cpu_core_large_mesh[0:8],  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_core_medium_mesh[0:8], t_cpu_core_medium_mesh[0:8], medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_cpu_core_small_mesh[0:8],  t_cpu_core_small_mesh[0:8],  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_core_single_node.png", dpi=300)

if plot_cpu_core:
  runtime_figure_cpu_core = create_runtime_plot(n_cpu_core_large_mesh, x_label_cpu_core)
  plt.plot(n_cpu_core_large_mesh,  t_cpu_core_large_mesh,  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_core_medium_mesh, t_cpu_core_medium_mesh, medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_cpu_core_small_mesh,  t_cpu_core_small_mesh,  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_core.png", dpi=300)

if plot_cpu_node:
  runtime_figure_cpu_node = create_runtime_plot(n_cpu_node_large_mesh, x_label_cpu_node)
  plt.plot(n_cpu_node_large_mesh,  t_cpu_node_large_mesh,  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_node_medium_mesh, t_cpu_node_medium_mesh, medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)   
  plt.plot(n_cpu_node_small_mesh,  t_cpu_node_small_mesh,  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_node.png", dpi=300)

if plot_gpu:
  runtime_figure_gpu = create_runtime_plot(n_gpu_large_mesh, x_label_gpu)
  plt.plot(n_gpu_large_mesh,  t_gpu_large_mesh,  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_gpu_medium_mesh, t_gpu_medium_mesh, medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_gpu_small_mesh,  t_gpu_small_mesh,  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_runtime_gpu.png", dpi=300)

if plot_compare_cpu_node_gpu:
  # CPU NODE and GPU Runtime Comparison
  runtime_figure_cpu_node = create_runtime_plot(n_cpu_node_large_mesh, x_label_cpu_node)
  
  #plt.plot(n_cpu_node_large_mesh,  t_cpu_node_large_mesh,  large_mesh_line_type+color_1+marker_type_1,  label='CPU Node - '+legend_label_large_mesh)
  #plt.plot(n_gpu_large_mesh[:-1], t_gpu_large_mesh[:-1], large_mesh_line_type+color_2+marker_type_2, label='GPU - '+legend_label_large_mesh)
  
  plt.plot(n_cpu_node_medium_mesh, t_cpu_node_medium_mesh, medium_mesh_line_type+color_3+marker_type_3, label='CPU Node - '+legend_label_medium_mesh)
  plt.plot(n_gpu_medium_mesh[:-1], t_gpu_medium_mesh[:-1], medium_mesh_line_type+color_2+marker_type_2, label='GPU - '+legend_label_medium_mesh)
  
  #plt.plot(n_cpu_node_small_mesh,  t_cpu_node_small_mesh,  small_mesh_line_type+color_3+marker_type_3,  label='CPU Node - '+legend_label_small_mesh)
  #plt.plot(n_gpu_small_mesh[:-1], t_gpu_small_mesh[:-1], small_mesh_line_type+color_2+marker_type_2, label='GPU - '+legend_label_small_mesh)
  
  plt.xlabel('N')
  plt.legend()
  plt.savefig("strong_scaling_runtime_cpu_node_and_gpu.png", dpi=300)

##############################################################################
# Speedup Figures
##############################################################################
if plot_cpu_core_single_node:
  speedup_figure_cpu_core_single_node  = create_speedup_plot(n_cpu_core_large_mesh[0:8], x_label_cpu_core)
  plt.plot(n_cpu_core_large_mesh[0:8],  speedup(t_cpu_core_large_mesh[0],  t_cpu_core_large_mesh[0:8]),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_core_medium_mesh[0:8], speedup(t_cpu_core_medium_mesh[0], t_cpu_core_medium_mesh[0:8]), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)  
  plt.plot(n_cpu_core_small_mesh[0:8],  speedup(t_cpu_core_small_mesh[0],  t_cpu_core_small_mesh[0:8]),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_core_single_node.png", dpi=300)

if plot_cpu_core:
  speedup_figure_cpu_core  = create_speedup_plot(n_cpu_core_large_mesh, x_label_cpu_core)
  plt.plot(n_cpu_core_large_mesh,  speedup(t_cpu_core_large_mesh[0],  t_cpu_core_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_core_medium_mesh, speedup(t_cpu_core_medium_mesh[0], t_cpu_core_medium_mesh), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)  
  plt.plot(n_cpu_core_small_mesh,  speedup(t_cpu_core_small_mesh[0],  t_cpu_core_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_core.png", dpi=300)

if plot_cpu_node:
  speedup_figure_cpu_node = create_speedup_plot(n_cpu_node_large_mesh, x_label_cpu_node)
  plt.plot(n_cpu_node_large_mesh,  speedup(t_cpu_node_large_mesh[0],  t_cpu_node_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_node_medium_mesh, speedup(t_cpu_node_medium_mesh[0], t_cpu_node_medium_mesh), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_cpu_node_small_mesh,  speedup(t_cpu_node_small_mesh[0],  t_cpu_node_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_node.png", dpi=300)

if plot_gpu:
  speedup_figure_gpu = create_speedup_plot(n_gpu_large_mesh, x_label_gpu)
  plt.plot(n_gpu_large_mesh,  speedup(t_gpu_large_mesh[0],  t_gpu_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_gpu_medium_mesh, speedup(t_gpu_medium_mesh[0], t_gpu_medium_mesh), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_gpu_small_mesh,  speedup(t_gpu_small_mesh[0],  t_gpu_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_speedup_gpu.png", dpi=300)

if plot_compare_cpu_node_gpu:
  # CPU NODE and GPU Speedup Comparison
  speedup_figure_cpu_node_and_gpu = create_speedup_plot(n_cpu_node_large_mesh, x_label_cpu_node)
  
  #plt.plot(n_cpu_node_large_mesh,  speedup(t_cpu_node_large_mesh[0],  t_cpu_node_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label='CPU Node - '+legend_label_large_mesh)
  #plt.plot(n_gpu_large_mesh[:-1],  speedup(t_gpu_large_mesh[0],  t_gpu_large_mesh[:-1]),  large_mesh_line_type+color_2+marker_type_2,  label='GPU - '+legend_label_large_mesh)
  
  plt.plot(n_cpu_node_medium_mesh, speedup(t_cpu_node_medium_mesh[0], t_cpu_node_medium_mesh), medium_mesh_line_type+color_3+marker_type_3, label='CPU Node - '+legend_label_medium_mesh)
  plt.plot(n_gpu_medium_mesh[:-1], speedup(t_gpu_medium_mesh[0], t_gpu_medium_mesh[:-1]), medium_mesh_line_type+color_2+marker_type_2, label='GPU - '+legend_label_medium_mesh)
  
  #plt.plot(n_cpu_node_small_mesh,  speedup(t_cpu_node_small_mesh[0],  t_cpu_node_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label='CPU Node - '+legend_label_small_mesh)
  #plt.plot(n_gpu_small_mesh[:-1],  speedup(t_gpu_small_mesh[0],  t_gpu_small_mesh[:-1]),  small_mesh_line_type+color_2+marker_type_2,  label='GPU - '+legend_label_small_mesh)    
  
  plt.xlabel('N')
  plt.legend()
  plt.savefig("strong_scaling_speedup_cpu_node_vs_gpu.png", dpi=300)

##############################################################################
# Efficiency Figures
##############################################################################
if plot_cpu_core_single_node:
  efficiency_figure_cpu_core_single_node = create_efficiency_plot(n_cpu_core_large_mesh[0:8], x_label_cpu_core)
  plt.plot(n_cpu_core_large_mesh[0:8],  efficiency(t_cpu_core_large_mesh[0],  t_cpu_core_large_mesh[0:8], n_cpu_core_large_mesh[0:8]),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_core_medium_mesh[0:8], efficiency(t_cpu_core_medium_mesh[0], t_cpu_core_medium_mesh[0:8], n_cpu_core_medium_mesh[0:8]), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_cpu_core_small_mesh[0:8],  efficiency(t_cpu_core_small_mesh[0],  t_cpu_core_small_mesh[0:8], n_cpu_core_small_mesh[0:8]),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_core_single_node.png", dpi=300)

if plot_cpu_core:
  efficiency_figure_cpu_core = create_efficiency_plot(n_cpu_core_large_mesh, x_label_cpu_core)
  plt.plot(n_cpu_core_large_mesh,  efficiency(t_cpu_core_large_mesh[0],  t_cpu_core_large_mesh, n_cpu_core_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_core_medium_mesh, efficiency(t_cpu_core_medium_mesh[0], t_cpu_core_medium_mesh, n_cpu_core_medium_mesh), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_cpu_core_small_mesh,  efficiency(t_cpu_core_small_mesh[0],  t_cpu_core_small_mesh, n_cpu_core_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_core.png", dpi=300)

if plot_cpu_node:
  efficiency_figure_cpu_node = create_efficiency_plot(n_cpu_node_large_mesh, x_label_cpu_node)
  plt.plot(n_cpu_node_large_mesh,  efficiency(t_cpu_node_large_mesh[0],  t_cpu_node_large_mesh, n_cpu_node_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_cpu_node_medium_mesh, efficiency(t_cpu_node_medium_mesh[0], t_cpu_node_medium_mesh, n_cpu_node_medium_mesh), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_cpu_node_small_mesh,  efficiency(t_cpu_node_small_mesh[0],  t_cpu_node_small_mesh, n_cpu_node_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_node.png", dpi=300)

if plot_gpu:
  efficiency_figure_gpu = create_efficiency_plot(n_gpu_large_mesh, x_label_gpu)
  plt.plot(n_gpu_large_mesh,  efficiency(t_gpu_large_mesh[0],  t_gpu_large_mesh, n_gpu_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label=legend_label_large_mesh)
  plt.plot(n_gpu_medium_mesh, efficiency(t_gpu_medium_mesh[0], t_gpu_medium_mesh, n_gpu_medium_mesh), medium_mesh_line_type+color_2+marker_type_2, label=legend_label_medium_mesh)
  plt.plot(n_gpu_small_mesh,  efficiency(t_gpu_small_mesh[0],  t_gpu_small_mesh, n_gpu_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label=legend_label_small_mesh)
  plt.legend()
  plt.savefig("strong_scaling_efficiency_gpu.png", dpi=300)

if plot_compare_cpu_node_gpu:
  # CPU NODE and GPU Efficiency Comparison
  efficiency_figure_cpu_node_and_gpu = create_efficiency_plot(n_cpu_node_large_mesh, x_label_cpu_node)
  
  #plt.plot(n_cpu_node_large_mesh,  efficiency(t_cpu_node_large_mesh[0],  t_cpu_node_large_mesh, n_cpu_node_large_mesh),  large_mesh_line_type+color_1+marker_type_1,  label='CPU Node - '+legend_label_large_mesh)
  #plt.plot(n_gpu_large_mesh[:-1],  efficiency(t_gpu_large_mesh[0],  t_gpu_large_mesh[:-1], n_gpu_large_mesh[:-1]),  large_mesh_line_type+color_2+marker_type_2,  label='GPU - '+legend_label_large_mesh)
  
  plt.plot(n_cpu_node_medium_mesh, efficiency(t_cpu_node_medium_mesh[0], t_cpu_node_medium_mesh, n_cpu_node_medium_mesh), medium_mesh_line_type+color_3+marker_type_3, label='CPU Node - '+legend_label_medium_mesh)
  plt.plot(n_gpu_medium_mesh[:-1], efficiency(t_gpu_medium_mesh[0], t_gpu_medium_mesh[:-1], n_gpu_medium_mesh[:-1]), medium_mesh_line_type+color_2+marker_type_2, label='GPU - '+legend_label_medium_mesh)
  
  #plt.plot(n_cpu_node_small_mesh,  efficiency(t_cpu_node_small_mesh[0],  t_cpu_node_small_mesh, n_cpu_node_small_mesh),  small_mesh_line_type+color_3+marker_type_3,  label='CPU Node - '+legend_label_small_mesh)
  #plt.plot(n_gpu_small_mesh[:-1],  efficiency(t_gpu_small_mesh[0],  t_gpu_small_mesh[:-1], n_gpu_small_mesh[:-1]),  small_mesh_line_type+color_2+marker_type_2,  label='GPU - '+legend_label_small_mesh)
  
  plt.xlabel('N')
  plt.legend()
  plt.savefig("strong_scaling_efficiency_cpu_node_vs_gpu.png", dpi=300)

plt.show()