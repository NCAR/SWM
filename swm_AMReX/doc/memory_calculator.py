import numpy as np

def memory_required_gb(n_double):
    # Each double precision number requires 8 bytes
    bytes_per_double = 8
    # Calculate the total number of bytes required
    total_bytes = n_double * bytes_per_double
    # Convert bytes to gigabytes
    bytes_per_gb = np.power(1024, 3)
    total_gb = total_bytes / float(bytes_per_gb)
    return total_gb

def print_memory_table(n_array):
    # Table header
    print(f"{'Nx':>10}{'Ny':>10}{'Memory (GB)':>15}")
    print("-" * 35)

    # Table Rows
    pow_low = 6
    pow_hi = 15
    for nx in [2**i for i in range(pow_low, pow_hi+1)]:
        ny = nx # Assume a 2D square domain
        n_cell = nx*ny
        memory_per_array_gb = memory_required_gb(n_cell)
        memory_gb = memory_per_array_gb * n_array
        print(f"{nx:>10} {ny:>10} {memory_gb:>15.6f}")

if __name__ == "__main__":
    for n_array in [1, 13]:
      print(f"Memory required for {n_array} arrays of double precision numbers")
      print_memory_table(n_array)
      print("")