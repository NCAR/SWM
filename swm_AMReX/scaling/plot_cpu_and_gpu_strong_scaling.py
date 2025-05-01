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
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14

# Font size for lebels in the legend
plt.rc('legend', fontsize=18) 

##############################################################################
# Strong Scaling Data to Plot
##############################################################################

# CPU - Total Runtime vs N CPU Core
n_cpu_core = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
t_cpu_core = np.array([79.71, 40.03, 20.26, 13.35, 9.476, 6.212, 6.373, 6.55, 3.102, 1.382, 0.4829, 0.1214, 0.08429])
x_label_cpu_core = r'$N_{\mathrm{CPU \ Core}}$'
legend_label_cpu_core = 'CPU Cores'
marker_type_cpu_core = '-.ro'

# CPU - Total Runtime vs N CPU Node (Ran on Derecho so each node has 128 cores)
n_cpu_node = np.array([1, 2, 4, 8, 16, 32])
t_cpu_node = np.array([6.558, 3.102, 1.382, 0.4829, 0.1214, 0.08429])
x_label_cpu_node = r'$N_{\mathrm{CPU \ NODE}}$'
legend_label_cpu_node = 'CPU Nodes'
marker_type_cpu_node = ':bD'

# GPU - Total Runtime vs N GPU
n_gpu = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu = np.array([1.166, 0.6295, 0.35, 0.2103, 0.1281, 0.09877, 0.06967])
x_label_gpu = r'$N_{\mathrm{GPU}}$'
legend_label_gpu = 'GPU'
marker_type_gpu = '-go'

##############################################################################
# Helper Functions
##############################################################################
def speedup(t_serial, t_parallel):
    return t_serial/np.array(t_parallel)

def efficiency(t_serial, t_parallel, n_proc):
    return speedup(t_serial, t_parallel)/np.array(n_proc)

def plot_runtime(n, t, marker_type,x_label_text, legend_label):
    figure = plt.figure()

    plt.plot(n,  t, marker_type, label=legend_label)

    # Set the x and y labels
    plt.xlabel(x_label_text)
    plt.ylabel(r'Runtime [s]')

    # Set x-axis ticks and labels to match the plotted points
    ax = figure.gca()  
    ax.set_xticks(n)  
    ax.set_xticklabels(n)  

    plt.legend()
    plt.tight_layout()

    return figure

def plot_speedup(n, t, marker_type, x_label_text, legend_label):
    figure = plt.figure()

    # Plot ideal speedup for reference. Should be a 45 degree line up and to the right.
    plt.plot([1,np.max(n)], [1,np.max(n)], '--k', label='Ideal') 

    # Plot the speedup. Assuming the first entry in t, t[0], is the reference time t_serial used to calculate speedup.
    plt.plot(n, speedup(t[0], t), marker_type, label=legend_label)

    # Set the x and y labels
    plt.xlabel(x_label_text)
    plt.ylabel(r'Speedup $= \frac{t_{\mathrm{serial}}}{t_\mathrm{parallel}}$')

    # Set x-axis ticks and labels to match the plotted points
    ax = figure.gca()  
    ax.set_xticks(n)  
    ax.set_xticklabels(n)  

    # Set aspect ratio to be equal to make the plot square. Makes it so the ideal speedup line is a 45 degree line.
    ax.set_aspect('equal') 

    plt.legend()
    plt.tight_layout()

    return figure

def plot_efficiency(n, t, marker_type,x_label_text, legend_label):
    figure = plt.figure()

    # Plot ideal efficiency for reference. Should be a horizontal line at y=1.
    plt.plot([1,np.max(n)], [1,1], '--k', label='Ideal')

    # Plot the efficiency. Assuming the first entry in t, t[0], is the reference time t_serial used to calculate efficiency.
    plt.plot(n, efficiency(t[0], t, n), marker_type, label=legend_label)

    # Set the x and y labels
    plt.xlabel(x_label_text)
    plt.ylabel(r'Efficiency $= \frac{\mathrm{Speedup}}{N_{\mathrm{proc}}} =  \frac{t_{\mathrm{serial}}}{N_{\mathrm{proc}}t_{\mathrm{parallel}}}$')

    # Set x-axis ticks and labels to match the plotted points
    ax1 = figure.gca()  
    ax1.set_xticks(n)  
    ax1.set_xticklabels(n)  

    plt.legend()
    plt.tight_layout()

    return figure


##############################################################################
# Runtime Figures
##############################################################################
#runtime_figure_cpu_core = plot_runtime(n_cpu_core, t_cpu_core, marker_type_cpu_core, x_label_cpu_core, legend_label_cpu_core)
#plt.savefig("strong_scaling_runtime_cpu_core.png", dpi=300)
#
#runtime_figure_cpu_node = plot_runtime(n_cpu_node, t_cpu_node, marker_type_cpu_node, x_label_cpu_node, legend_label_cpu_node)
#plt.savefig("strong_scaling_runtime_cpu_node.png", dpi=300)
#
#runtime_figure_gpu= plot_runtime(n_gpu, t_gpu, marker_type_gpu, x_label_gpu, legend_label_gpu)
#plt.savefig("strong_scaling_runtime_gpu.png", dpi=300)

runtime_figure_cpu_node_and_gpu = plot_runtime(n_cpu_node, t_cpu_node, marker_type_cpu_node, x_label_cpu_node, legend_label_cpu_node)
#plt.plot(n_gpu, t_gpu, marker_type_gpu, label=legend_label_gpu)
plt.plot(n_gpu[:-1], t_gpu[:-1], marker_type_gpu, label=legend_label_gpu)
plt.xlabel('N')
plt.legend()
plt.savefig("strong_scaling_runtime_cpu_node_and_gpu.png", dpi=300)

##############################################################################
# Speedup Figures
##############################################################################
#speedup_figure_cpu_core = plot_speedup(n_cpu_core, t_cpu_core, marker_type_cpu_core, x_label_cpu_core, legend_label_cpu_core)
#plt.savefig("strong_scaling_speedup_cpu_core.png", dpi=300)
#
#speedup_figure_cpu_node = plot_speedup(n_cpu_node, t_cpu_node, marker_type_cpu_node, x_label_cpu_node, legend_label_cpu_node)
#plt.savefig("strong_scaling_speedup_cpu_node.png", dpi=300)
#
#speedup_figure_gpu = plot_speedup(n_gpu, t_gpu, marker_type_gpu, x_label_gpu, legend_label_gpu)
#plt.savefig("strong_scaling_speedup_cpu_node.png", dpi=300)

speedup_figure_cpu_node_and_gpu = plot_speedup(n_cpu_node, t_cpu_node, marker_type_cpu_node, x_label_cpu_node, legend_label_cpu_node)
#plt.plot(n_gpu, speedup(t_gpu[0], t_gpu), marker_type_gpu, label=legend_label_gpu)
plt.plot(n_gpu[:-1], speedup(t_gpu[0], t_gpu[:-1]), marker_type_gpu, label=legend_label_gpu)
plt.xlabel('N')
plt.legend()
plt.savefig("strong_scaling_speedup_cpu_node_vs_gpu.png", dpi=300)

##############################################################################
# Efficiency Figures
##############################################################################

#efficiency_figure_cpu_core = plot_efficiency(n_cpu_core, t_cpu_core, marker_type_cpu_core, x_label_cpu_core, legend_label_cpu_core)
#plt.savefig("strong_scaling_efficiency_cpu_core.png", dpi=300)
#
#efficiency_figure_cpu_node = plot_efficiency(n_cpu_node, t_cpu_node, marker_type_cpu_node, x_label_cpu_node, legend_label_cpu_node)
#plt.savefig("strong_scaling_efficiency_cpu_node.png", dpi=300)
#
#efficiency_figure_gpu= plot_efficiency(n_gpu, t_gpu, marker_type_gpu, x_label_gpu, legend_label_gpu)
#plt.savefig("strong_scaling_efficiency_gpu.png", dpi=300)

efficiency_figure_cpu_node_and_gpu = plot_efficiency(n_cpu_node, t_cpu_node, marker_type_cpu_node, x_label_cpu_node, legend_label_cpu_node)
#plt.plot(n_gpu, efficiency(t_gpu[0], t_gpu, n_gpu), marker_type_gpu, label=legend_label_gpu)
plt.plot(n_gpu[:-1], efficiency(t_gpu[0], t_gpu[:-1], n_gpu[:-1]), marker_type_gpu, label=legend_label_gpu)
plt.xlabel('N')
plt.legend()
plt.savefig("strong_scaling_efficiency_cpu_node_and_gpu.png", dpi=300)

##############################################################################
# Extra Figures Comparing Different Mesh Sizes on GPUs 
##############################################################################

# GPU - Large Mesh - 16384
n_gpu_large_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_large_mesh = np.array([4.633, 2.381, 1.254, 0.6688, 0.3576, 0.2139, 0.1288])
marker_type_gpu_large_mesh = '-.bD'
legend_label_gpu_large_mesh = '16384 x 16384'

# GPU - Default Mesh - 8192
n_gpu_default_mesh = n_gpu
t_gpu_default_mesh = t_gpu
marker_type_gpu_default_mesh = '-go'
legend_label_gpu_default_mesh = '8192 x 8192'

# GPU - Small Mesh - 4096
n_gpu_small_mesh = np.array([1, 2, 4, 8, 16, 32, 64])
t_gpu_small_mesh = np.array([0.3004, 0.1902, 0.1207, 0.09364, 0.06838, 0.06885, 0.05205])
marker_type_gpu_small_mesh = ':rs'
legend_label_gpu_small_mesh = '4096 x 4096'

speedup_figure_gpu_multi_mesh_size = plot_speedup(n_gpu_large_mesh, t_gpu_large_mesh, marker_type_gpu_large_mesh, x_label_gpu, legend_label_gpu_large_mesh)
plt.plot(n_gpu, speedup(t_gpu[0], t_gpu), marker_type_gpu_default_mesh, label=legend_label_gpu_default_mesh)
plt.plot(n_gpu_small_mesh, speedup(t_gpu_small_mesh[0], t_gpu_small_mesh), marker_type_gpu_small_mesh, label=legend_label_gpu_small_mesh)
plt.legend()
plt.savefig("strong_scaling_speedup_gpu_mesh_compare.png", dpi=300)

speedup_figure_gpu_multi_mesh_size = plot_efficiency(n_gpu_large_mesh, t_gpu_large_mesh, marker_type_gpu_large_mesh, x_label_gpu, legend_label_gpu_large_mesh)
plt.plot(n_gpu, efficiency(t_gpu[0], t_gpu, n_gpu), marker_type_gpu_default_mesh, label=legend_label_gpu_default_mesh)
plt.plot(n_gpu_small_mesh, efficiency(t_gpu_small_mesh[0], t_gpu_small_mesh, n_gpu), marker_type_gpu_small_mesh, label=legend_label_gpu_small_mesh)
plt.legend()
plt.savefig("strong_scaling_speedup_gpu_mesh_compare.png", dpi=300)

plt.show()