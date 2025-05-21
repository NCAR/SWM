import os
import re
import glob
import argparse

def extract_total_runtime(amrex_output_txt_file):
    # Read the output file
    with open(amrex_output_txt_file, "r") as f:
        lines = f.readlines()

    # Find the line containing the total runtime
    runtime_line = None
    for line in lines:
        if "TinyProfiler total time across processes" in line:
            runtime_line = line.strip()
            break

    if runtime_line is None:
        exit(f"Error: Could not find runtime information in {amrex_output_txt_file}")

    # Extract runtime_min, runtime_avg, runtime_max using regex... kind of a hack... maybe comeback to this later
    match = re.search(r"total time across processes .*: ([^ ]+) \.\.\. ([^ ]+) \.\.\. ([^ ]+)", runtime_line)

    if match:
        runtime_min = float(match.group(1))
        runtime_avg = float(match.group(2))
        runtime_max = float(match.group(3))
    else:
        exit(f"Error: Could not parse line {runtime_line} for runtime information from file {amrex_output_txt_file}")
    
    return runtime_min, runtime_avg, runtime_max


def extract_timer_data(amrex_output_txt_file, timer_names, timer_type):
    # Define the header pattern based on the timer type
    if timer_type == "exclusive":
        table_header_pattern = r"^Name\s+NCalls\s+Excl\. Min\s+Excl\. Avg\s+Excl\. Max\s+Max %$"
    elif timer_type == "inclusive":
        table_header_pattern = r"^Name\s+NCalls\s+Incl\. Min\s+Incl\. Avg\s+Incl\. Max\s+Max %$"
    else:
        raise ValueError(f"Invalid timer type: {timer_type}")

    with open(amrex_output_txt_file, "r") as f:
        lines = f.readlines()

    # Find the header line and extract the relevant table
    table_start = None
    for i, line in enumerate(lines):
        if re.match(table_header_pattern, line):
            table_start = i + 1  # Start after the header line
            break

    if table_start is None:
        exit(f"Error: Could not find {timer_type} timer table in {amrex_output_txt_file}")

    # Extract the table (assumes the table continues for 10 lines after the header)
    table_lines = lines[table_start:table_start + 10]

    # Parse the table and extract data for the specified timers
    timer_data = {}
    for timer_name in timer_names:
        for line in table_lines:
            if line.startswith(timer_name):
                parts = line.split()
                runtime_min = float(parts[2])
                runtime_avg = float(parts[3])
                runtime_max = float(parts[4])
                percent_max = float(parts[5].strip('%'))  # Remove the percent sign at the end
                timer_data[timer_name] = (runtime_min, runtime_avg, runtime_max, percent_max)
                break

    return timer_data

def main():

    ###########################################################################
    # Hard coded input
    ###########################################################################

    # Names of timers that will be extracted from the AMReX profiler output and saved to CSV files
    timer_names=["main", "UpdateIntermediateVariables", "UpdateNewVariables", "UpdateOldVariables", "FabArray::FillBoundary"]

    ###########################################################################
    # Command line input 
    ###########################################################################
    parser = argparse.ArgumentParser(description="Process AMReX scaling results.")

    parser.add_argument(
        "scaling_results_root_dir",
        type=str,
        help="Path to the root directory containing scaling results (e.g., it should have sub directories named n_proc_XXXX for various values of XXXX.)."
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default=os.getcwd(),
        help="Path to the output directory where resulting csv and text files will be saved. Defaults to the current working directory."
    )

    args = parser.parse_args()

    ###########################################################################
    # Setup
    ###########################################################################

    # Make sure the user specified input directory exists
    if not os.path.isdir(args.scaling_results_root_dir):
        print(f"Error: The directory you supplied via the command line '{args.scaling_results_root_dir}' does not exist.")
        exit(1)

    # Make sure the output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    
    # Print headers for each CSV file that are going to be create
    with open(os.path.join(args.outdir, "total_runtime.csv"), "w") as f:
        f.write("n_proc, run_idx, runtime_min [s], runtime_avg [s], runtime_max [s]\n")
    
    timer_types = ["inclusive", "exclusive"]
    for timer_type in timer_types:
        for timer_name in timer_names:
            file_path = os.path.join(args.outdir, f"{timer_type}_runtime_{timer_name}.csv")
            with open(file_path, "w") as f:
                f.write("n_proc, run_idx, runtime_min [s], runtime_avg [s], runtime_max [s], max_percent\n")
    
    ###########################################################################
    # Setup
    ###########################################################################

    # Loop over directories matching the pattern n_proc_XXXX
    np_dirs = sorted(glob.glob(os.path.abspath(os.path.join(args.scaling_results_root_dir, "n_proc_*"))))

    for np_dir in np_dirs:
        # Make sure this is a valid directory
        if not os.path.isdir(np_dir):
            os.exit(f"Error: The path '{np_dir}' is not a directory.")

        # Match the pattern n_proc_XXXX and extract the number at the end
        np_match = re.match(r"n_proc_(\d+)", os.path.basename(np_dir))

        # Make sure there was a pattern match for the number of processors
        if not np_match:
            os.exit(f"Error: The directory name '{np_dir}' does not match the expected pattern 'n_proc_XXXX'.")

        n_proc = int(np_match.group(1))

        print(f"Found directory: {np_dir}, n_proc: {n_proc}")

        run_dirs = sorted(glob.glob(os.path.abspath(os.path.join(np_dir, "run_number_*"))))

        for run_dir in run_dirs:
            if not os.path.isdir(run_dir):
                os.exit(f"Error: The path '{run_dir}' is not a directory.")
            # Match the pattern run_number_XXXX and extract the number at the end
            run_match = re.match(r"run_number_(\d+)", os.path.basename(run_dir))
            # Make sure there was a pattern match for the run number
            if not run_match:
                os.exit(f"Error: The directory name '{run_dir}' does not match the expected pattern 'run_number_XXXX'.")
            run_number = int(run_match.group(1))
            print(f"Found directory: {run_dir}, run_number: {run_number}")

            amrex_output_txt_file = os.path.join(run_dir, "output.txt")
            
            # Write the extracted runtime data to the total_runtime.csv file
            (runtime_min, runtime_avg, runtime_max) = extract_total_runtime(amrex_output_txt_file)
            with open(os.path.join(args.outdir, "total_runtime.csv"), "a") as f:
                f.write(f"{n_proc}, {run_number}, {runtime_min}, {runtime_avg}, {runtime_max}\n")


            # Write the specific timers we asked for to total_runtime.csv file
            for timer_type in timer_types:
                timer_data = extract_timer_data(amrex_output_txt_file, timer_names, timer_type)

                # Write the extracted data to CSV files
                for timer_name, (runtime_min, runtime_avg, runtime_max, percent_max) in timer_data.items():
                    csv_file = os.path.join(args.outdir, f"{timer_type}_runtime_{timer_name}.csv")
                    with open(csv_file, "a") as f:
                        f.write(f"{n_proc}, {run_number}, {runtime_min}, {runtime_avg}, {runtime_max}, {percent_max}\n")

if __name__ == "__main__":
    main()