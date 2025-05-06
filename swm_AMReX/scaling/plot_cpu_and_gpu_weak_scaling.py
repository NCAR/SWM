import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

##############################################################################
# User Input
##############################################################################
# Controls how the plots looks

# Set default figure size in inches
plt.rcParams["figure.figsize"] = (8, 8) 

# Font size for x and y axis labels
mpl.rcParams['axes.labelsize'] = 20

# Font size for the x and y axis tick labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

# Font size for lebels in the legend
plt.rc('legend', fontsize=18) 

##############################################################################
# Weak Scaling Data 
##############################################################################

# Note: 
#     small  mesh is   4096 x  4096
#     medium mesh is   8192 x  8192
#     large  mesh is  16384 x 16384

# CPU - Total Runtime vs N CPU Node (Ran on Derecho so each node has 128 cores)... the times are the same as cpu_core data above but just from 128 cores to the end.
n_cpu_node_small_mesh = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node_small_mesh = np.array([1.349, 1.349, 1.354, 1.356, 1.358, 1.356])

n_cpu_node_medium_mesh = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node_medium_mesh = np.array([6.487, 6.466, 6.474, 6.483, 6.476, 6.489])

n_cpu_node_large_mesh = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node_large_mesh = np.array([26.13, 26.12, 26.13, 26.2, 26.18, 26.16])

# GPU - Total Runtime vs N GPU
n_gpu_small_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_small_mesh = np.array([0.3004, 0.3315, 0.35, 0.3569, 0.3573, 0.3604, 0.3598])

n_gpu_medium_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_medium_mesh = np.array([1.16, 1.201, 1.246, 1.254, 1.252, 1.252, 1.254])

n_gpu_large_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_large_mesh = np.array([4.598, 4.661, 4.763, 4.774, 4.772, 4.773, 4.78])

##############################################################################
# Helper Functions
##############################################################################

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

legend_label_cpu_node = 'CPU Nodes'
legend_label_gpu = 'GPU'
legend_label_large_mesh  = '16384 x 16384'
legend_label_medium_mesh =   '8192 x 8192'
legend_label_small_mesh  =   '4096 x 4096'

x_label_cpu_node=r'$N_{\mathrm{CPU \ Node}}$'
x_label_gpu=r'$N_{\mathrm{GPU}}$'

# Useful for line type, color, and marker types
small_mesh_line_type = ':ro'
medium_mesh_line_type = '-.bs'
large_mesh_line_type = '-gd'

gpu_line_type = '-go'
cpu_line_type = '-.rs'

###############################################################################
# Runtime - Mesh Size
###############################################################################

runtime_figure_cpu_node = create_runtime_plot(n_cpu_node_small_mesh, x_label_cpu_node)
plt.plot(n_cpu_node_large_mesh, t_cpu_node_large_mesh, large_mesh_line_type, label=legend_label_large_mesh)
plt.plot(n_cpu_node_medium_mesh, t_cpu_node_medium_mesh, medium_mesh_line_type, label=legend_label_medium_mesh)
plt.plot(n_cpu_node_small_mesh, t_cpu_node_small_mesh, small_mesh_line_type, label=legend_label_small_mesh)
plt.legend()
plt.savefig("weak_scaling_runtime_cpu_node.png", dpi=300)

runtime_figure_gpu = create_runtime_plot(n_cpu_node_large_mesh, x_label_gpu)
plt.plot(n_gpu_large_mesh, t_gpu_large_mesh, large_mesh_line_type, label=legend_label_large_mesh)
plt.plot(n_gpu_medium_mesh, t_gpu_medium_mesh, medium_mesh_line_type, label=legend_label_medium_mesh)
plt.plot(n_gpu_small_mesh, t_gpu_small_mesh, small_mesh_line_type, label=legend_label_small_mesh)
plt.legend()
plt.savefig("weak_scaling_runtime_gpu.png", dpi=300)

###############################################################################
# Runtime - CPU vs GPU
###############################################################################

runtime_figure_compare_cpu_gpu_4096 = create_runtime_plot(n_cpu_node_small_mesh, 'N')
plt.plot(n_cpu_node_small_mesh, t_cpu_node_small_mesh, cpu_line_type, label='CPU Nodes -'+legend_label_small_mesh)
plt.plot(n_gpu_small_mesh[:-1], t_gpu_small_mesh[:-1], gpu_line_type, label='GPU -'+legend_label_small_mesh)
plt.legend()
plt.savefig("weak_scaling_runtime_compare_cpu_gpu_4096.png", dpi=300)

runtime_figure_compare_cpu_gpu_8192 = create_runtime_plot(n_cpu_node_small_mesh, 'N')
plt.plot(n_cpu_node_medium_mesh, t_cpu_node_medium_mesh, cpu_line_type, label='CPU Nodes -'+legend_label_medium_mesh)
plt.plot(n_gpu_medium_mesh[:-1], t_gpu_medium_mesh[:-1], gpu_line_type, label='GPU -'+legend_label_medium_mesh)
plt.legend()
plt.savefig("weak_scaling_runtime_compare_cpu_gpu_8192.png", dpi=300)

runtime_figure_compare_cpu_gpu_16384 = create_runtime_plot(n_cpu_node_small_mesh, 'N')
plt.plot(n_cpu_node_large_mesh, t_cpu_node_large_mesh, cpu_line_type, label='CPU Nodes -'+legend_label_large_mesh)
plt.plot(n_gpu_large_mesh[:-1], t_gpu_large_mesh[:-1], gpu_line_type, label='GPU -'+legend_label_large_mesh)
plt.legend()
plt.savefig("weak_scaling_runtime_compare_cpu_gpu_16384.png", dpi=300)

plt.show()