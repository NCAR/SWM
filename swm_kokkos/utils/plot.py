import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

xt=["128x128","256x256","512x512","1024x1024","2048x2048"] # x axis labels
xx=[128*128,256*256,512*512,1024*1024,2048*2048] # x axis values
y0=[0.584025,2.319167,21.536601,87.819181,320.092663] # C, -O2

# Read the CSV file into a numpy array
# y1 = pd.read_csv('./data/O2/kokkos_rangepolicy.csv', header=None, sep=',').values
# y1 = pd.read_csv('./data/O2/kokkos_teampolicy.csv', header=None, sep=',').values
y1 = pd.read_csv('./data/O2/kokkos_mdrangepolicy.csv', header=None, sep=',').values
# y4 = pd.read_csv('./data/O3/kokkos_mdrangepolicy.csv', header=None, sep=',').values

# Statistics: compute element-wise ratio and the average ratio
ratios = np.array(y0) / np.array(y1[5,:])
avg_ratio = np.mean(ratios)
print("C/Kokkos ratio:", ratios)
print("Average C/Kokkos ratio:", avg_ratio)

# Figure 1: C vs. Kokkos RangePolicy (serial)
plt.figure(figsize=(8, 6))
plt.plot(xx, y0, marker='o', linewidth=2., color='k', label='C')
plt.plot(xx, y1[0,:], marker='s', linewidth=2., color='r', label='Kokkos (Serial)')
plt.xscale('log')
plt.yscale('log')
plt.xticks(xx, xt)
plt.ylim(0.1, 400)
plt.xlabel('Domain Size')
plt.ylabel('Execution Time (s)')
plt.title('Performance comparison on a single AMD Milan CPU core (log-log scale)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figure1.pdf")

# Figure 2: Kokkos RangePolicy (multi-threaded)
plt.figure(figsize=(8, 6))
plt.plot(xx, y0, marker='o', linewidth=2., color='k', label='C')
markers = ['s', '^', 'v', 'D', 'P', '*']
colors = ['r', 'g', 'c', 'm', 'y', 'b']
labels = ['Kokkos (Serial)', 'Kokkos (OpenMP, 1 thread)', 'Kokkos (OpenMP, 2 threads)',
          'Kokkos (OpenMP, 4 threads)', 'Kokkos (OpenMP, 8 threads)', 'Kokkos (OpenMP, 128 threads)']
for i in range(6):
    plt.plot(xx, y1[i,:], marker=markers[i], linewidth=2., color=colors[i], label=labels[i])
plt.xscale('log')
plt.yscale('log')
plt.xticks(xx, xt)
plt.ylim(0.1, 400)
plt.xlabel('Domain Size')
plt.ylabel('Execution Time (s)')
plt.title('Performance comparison on an AMD Milan CPU node (log-log scale)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figure2.pdf")

# Figure 3: Kokkos RangePolicy (CUDA)
plt.figure(figsize=(8, 6))
plt.plot(xx, y0, marker='o', linewidth=2., color='k', label='C')
plt.plot(xx, y1[0,:], marker='s', linewidth=2., color='r', label='Kokkos (Serial)')
plt.plot(xx, y1[5,:], marker='^', linewidth=2., color='g', label='Kokkos (OpenMP, 128 threads)')
plt.plot(xx, y1[6,:], marker='v', linewidth=2., color='c', label='Kokkos (CUDA, 1 GPU)')
plt.xscale('log')
plt.yscale('log')
plt.xticks(xx, xt)
plt.ylim(0.1, 400)
plt.xlabel('Domain Size')
plt.ylabel('Execution Time (s)')
plt.title('Performance comparison on an AMD Milan CPU node or one NVIDIA A100 GPU (log-log scale)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figure3.pdf")