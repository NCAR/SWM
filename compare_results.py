import numpy as np
import argparse

if __name__ == "__main__":
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-r", "--ref",
        required=True,
        help="Path to directory containing reference .csv files")
    ap.add_argument("-t", "--test",
        required=True,
        help="Path to directory containing test .csv files")
    ap.add_argument("-i","--refid",
        required=True,
        help="Unique identifier of reference files you want to compare")
    ap.add_argument("-d","--testid",
        required=True,
        help="Unique identifier of test files you want to compare")
    args = vars(ap.parse_args())

    dims = [48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072]

    ref_dir = args['ref']
    test_dir = args['test']
    ref_id = str(args['refid'])
    test_id = str(args['testid'])
    out_file_name = test_dir + "/L_inf_norms." + test_id + ".csv"
    out_file = open(out_file_name,'w')

    # Label columns
    out_file.write("Prb Dims,L_inf Norm\n")
    for dim in dims:
        # Write dimensions to file
        out_file.write("{:d}x{:d},".format(dim,dim))

        # Get reference result file for this dimension
        _, ref_tail = ref_dir.rsplit("/",1)
        ref_file_name = ref_dir + "/swm_end." + str(dim) + "." + str(dim) + "." + ref_tail + "." + ref_id + ".csv"
        ref_data = np.genfromtxt(ref_file_name, delimiter=',')

        # Get test result file for this dimension
        _, test_tail = test_dir.rsplit("/",1)
        test_file_name = test_dir + "/swm_end." + str(dim) + "." + str(dim) + "." + test_tail + "." + test_id + ".csv"
        test_data = np.genfromtxt(test_file_name, delimiter=',')

        # Perform L infinity norm on difference between results
        diff = np.subtract(ref_data,test_data)
        L_inf_norm = np.linalg.norm(diff, ord=np.inf)

        # Write L infinity norm to file
        out_file.write("{:.15f}\n".format(L_inf_norm))
        print("{:d}x{:d}: {:.15f}".format(dim,dim,L_inf_norm))
