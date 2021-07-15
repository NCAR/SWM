import numpy as np
import argparse


if __name__ == "__main__":
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-r", "--ref",
        required=True,
        help="Path to reference .csv file")
    ap.add_argument("-t", "--test",
        required=True,
        help="Path to test .csv file")
    args = vars(ap.parse_args())

    ref_data = np.genfromtxt(args['ref'], delimiter=',')
    test_data = np.genfromtxt(args['test'], delimiter=',')

    # Perform L infinity norm on difference between results
    diff = np.subtract(ref_data,test_data)
    
    L_inf_norm = np.linalg.norm(diff, ord=np.inf)

    print("L_inf norm: {:.15f}".format(L_inf_norm))