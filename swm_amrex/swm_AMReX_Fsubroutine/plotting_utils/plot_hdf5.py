import argparse
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Read datasets from an HDF5 file. Save a plot of each dataset as a PNG file.')
    parser.add_argument('hdf5file', type=str, help='Path to the HDF5 file to read.')
    parser.add_argument('--output_dir', type=str, default=os.getcwd(), help='Directory to save the output PNG files to.')

    args = parser.parse_args()

    title_fontsize = 18
    label_fontsize = 12

    with h5py.File(args.hdf5file, 'r') as f:
        # Read the time attribute
        time_step = f.attrs['time_step']
        time = f.attrs['time']

        # Iterate over all datasets in the file
        for dataset_name in f:
            dataset = f[dataset_name]
            data = dataset[:]
            data_2d = np.array(data).reshape(dataset.shape)

            plt.figure()
            plt.imshow(data_2d, aspect='auto', cmap='viridis')
            plt.colorbar()
            plt.title(f"${dataset_name}$", fontsize=title_fontsize)
            plt.xlabel('$x$', fontsize=label_fontsize)
            plt.ylabel('$y$', fontsize=label_fontsize)

            # Print time step in the top left corner
            plt.gcf().text(0, 1.01,
                           f"Time Step: {time_step}",
                           transform=plt.gca().transAxes,
                           fontsize=label_fontsize,
                           verticalalignment='bottom',
                           horizontalalignment='left')

            # Print time in the top right corner
            plt.gcf().text(1.00, 1.01,
                           f"Time: {time}",
                           transform=plt.gca().transAxes,
                           fontsize=label_fontsize,
                           verticalalignment='bottom',
                           horizontalalignment='right')

            output_file = os.path.join(args.output_dir, f"{dataset_name}_{time_step:05d}.png")
            plt.savefig(output_file)

            print(f"Saved plot to {output_file}")