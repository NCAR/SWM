import argparse

parser = argparse.ArgumentParser(description="Shallow Water Model")
parser.add_argument('--M', type=int, default=64, help='Number of points in the x direction')
parser.add_argument('--N', type=int, default=64, help='Number of points in the y direction')
parser.add_argument('--L_OUT', type=bool, default=True, help='a boolean for L_OUT')
parser.add_argument('--ITMAX', type=int, default=4000, help='a boolean for L_OUT')


args = parser.parse_args()

# Initialize model parameters
M = args.M
N = args.N
M_LEN = M + 1
N_LEN = N + 1
L_OUT = args.L_OUT
VAL=True
VAL_DEEP=True
VIS = False
VIS_DT = 10
ITMAX = args.ITMAX
dt = 90.
dt = dt
dx = 100000.
dy = 100000.
fsdx = 4. / (dx)
fsdy = 4. / (dy)
a = 1000000.
alpha = 0.001
