import argparse
import jax.numpy as jnp


def main():
    parser = argparse.ArgumentParser(description="Shallow Water Model - JAX implementation")
    parser.add_argument('--M', type=int, default=64, help='Number of points in the x direction')
    parser.add_argument('--N', type=int, default=128, help='Number of points in the y direction')
    
    args = parser.parse_args()
    print(args.M, args.N)