import matplotlib.pyplot as plt

xt=["128x128","256x256","512x512","1024x1024","2048x2048"] # x axis labels
xx=[128*128,256*256,512*512,1024*1024,2048*2048] # x axis values
y1=[0.764340,2.776699,11.449694,49.772942,215.709304] # kokkos
y2=[0.584025,2.319167,21.536601,87.819181,320.092663] # C

plt.figure(figsize=(8, 6))
plt.plot(xx, y1, marker='o', label='Kokkos')
plt.plot(xx, y2, marker='s', label='C')
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
plt.savefig("performance_comparison.pdf")