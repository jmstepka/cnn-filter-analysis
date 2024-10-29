import argparse
import os
from src.create_pfm import create_pfm_matrices
from src.calculate_ic import calculate_ic
from src.ic_analysis import ic_plots
from src.calculate_kl import calculate_kl_distance
from src.kl_analysis import perform_kl_analysis

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def generate_run_id(base_dir="./data"):
    """
    Generate a new unique run_id in a simpler format: run_1, run_2, etc.

    Args:
        base_dir (str): The base directory where runs are stored.

    Returns:
        str: A new run_id in the format 'run_N'.
    """
    existing_runs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    run_numbers = []
    for run in existing_runs:
        if run.startswith('run_'):
            try:
                run_num = int(run.split('_')[1])
                run_numbers.append(run_num)
            except ValueError:
                continue
    if run_numbers:
        next_run_number = max(run_numbers) + 1
    else:
        next_run_number = 1
    return f"run_{next_run_number}"

def create_directories_for_run(run_id, base_dir="./data"):
    """Create the main directory for the run and its subdirectories."""
    run_dir = os.path.join(base_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)
    
    # Creating subdirectories in the run_dir
    pfm_dir = os.path.join(run_dir, "pfm_matrices")
    ic_dir = os.path.join(run_dir, "calculated_ic")
    kl_dir = os.path.join(run_dir, "kl_distance")
    ic_plots_dir = os.path.join(run_dir, "ic_plots")
    kl_plots_dir = os.path.join(run_dir, "kl_plots")
    
    os.makedirs(pfm_dir, exist_ok=True)
    os.makedirs(ic_dir, exist_ok=True)
    os.makedirs(kl_dir, exist_ok=True)
    os.makedirs(ic_plots_dir, exist_ok=True)
    os.makedirs(kl_plots_dir, exist_ok=True)
    
    return run_dir, pfm_dir, ic_dir, kl_dir, ic_plots_dir, kl_plots_dir

def compute_default_ic_threshold(k):
    """Compute the default IC threshold as (k + 1) / 2."""
    return (k + 1) / 2

def main():
    parser = argparse.ArgumentParser(description="Bioinformatics CLI for matrix operations and information theory analysis.")
    
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands for different operations")

    # Input and output options for each command
    def add_input_output_options(subparser):
        subparser.add_argument("--input_dir", required=True, help="Input directory for the subcommand")
        subparser.add_argument("--output_dir", required=True, help="Output directory for the subcommand")

    # Sub-command: create_pfm_matrices
    parser_create = subparsers.add_parser("create_pfm_matrices", help="Create PFM matrices from input")
    add_input_output_options(parser_create)
    parser_create.add_argument("--filter_mean_threshold", type=float, default=1e-5, help="Mean threshold for filters to be used")
    # Note: 'filter_mean_threshold' is optional with default value 1e-5

    # Sub-command: calculate_ic
    parser_ic = subparsers.add_parser("calculate_ic", help="Calculate Information Content from PFM matrices")
    add_input_output_options(parser_ic)
    parser_ic.add_argument("-k", "--kmer_length", type=int, default=7, help="Length of k-mers to use in IC calculation")

    # Sub-command: ic_analysis
    parser_analysis = subparsers.add_parser("ic_analysis", help="Perform IC analysis with a threshold")
    add_input_output_options(parser_analysis)
    parser_analysis.add_argument("-k", "--kmer_length", type=int, default=7, help="Length of k-mers to use in IC analysis")

    # Sub-command: calculate_kl
    parser_kl = subparsers.add_parser("calculate_kl", help="Calculate KL divergence between two PFM matrices")
    add_input_output_options(parser_kl)
    parser_kl.add_argument("-ic_thresh", "--ic_threshold", type=float, default=None, help="IC threshold for filtering motifs")
    parser_kl.add_argument("-pfm", "--pfm_folder", required=True, help="Folder containing the PFM matrices")
    parser_kl.add_argument("--hocomoco_path", required=True, help="Path to HOCOMOCO models in meme format")
    # Note: 'ic_threshold' is optional

    # Sub-command: kl_analysis
    parser_kl_analysis = subparsers.add_parser("kl_analysis", help="Perform analysis on KL divergence results")
    add_input_output_options(parser_kl_analysis)
    parser_kl_analysis.add_argument("--tf_family_file", required=True, help="Path to the TF family CSV file")

    # Updated mutually exclusive group to include --test_all
    group_kl_analysis = parser_kl_analysis.add_mutually_exclusive_group(required=True)
    group_kl_analysis.add_argument("-f", "--family_files", nargs="+", help="List of motif family files")
    group_kl_analysis.add_argument("-c", "--class_descriptors", nargs="+", help="List of class descriptors in 'x.y.z' format")
    group_kl_analysis.add_argument("--test_all", action='store_true', help="Test all superclasses or classes")

    parser_kl_analysis.add_argument("-n", "--family_names", nargs="+", help="List of family names corresponding to family files")
    parser_kl_analysis.add_argument("--test_level", choices=['superclass', 'class'], default='superclass', help="Level to test when testing all")
    parser_kl_analysis.add_argument("--correction_method", default='fdr_bh', help="Method for multiple testing correction (e.g., 'bonferroni', 'fdr_bh')")

    # New sub-command: create_new_run
    parser_new_run = subparsers.add_parser("create_new_run", help="Generate a new run and run all steps of the pipeline")
    parser_new_run.add_argument("--input_dir", required=True, help="Input directory for the analysis")
    parser_new_run.add_argument("-k", "--kmer_length", type=int, default=7, help="Length of k-mers to use for analysis")
    parser_new_run.add_argument("--tf_family_file", required=True, help="Path to the TF family CSV file")
    parser_new_run.add_argument("--hocomoco_path", required=True, help="Path to HOCOMOCO models in meme format")
    parser_new_run.add_argument("--filter_mean_threshold", type=float, default=1e-5, help="Mean threshold for filters to be used")
    parser_new_run.add_argument("--run_name", help="Name of the run. If not specified, a unique run_id will be generated.")
    parser_new_run.add_argument("-t", "--ic_threshold", default=None, help="IC threshold for analysis; can be 'manual' to input after IC plots are generated, or a numeric value. Default is (k + 1) / 2")
    # Note: 'ic_threshold' is optional

    # Updated mutually exclusive group to include --test_all
    group_new_run = parser_new_run.add_mutually_exclusive_group(required=True)
    group_new_run.add_argument("-f", "--family_files", nargs="+", help="List of motif family files")
    group_new_run.add_argument("-c", "--class_descriptors", nargs="+", help="List of class descriptors in 'x.y.z' format")
    group_new_run.add_argument("--test_all", action='store_true', help="Test all superclasses or classes")

    parser_new_run.add_argument("-n", "--family_names", nargs="+", help="List of family names corresponding to family files")
    parser_new_run.add_argument("--test_level", choices=['superclass', 'class'], default='superclass', help="Level to test when testing all")
    parser_new_run.add_argument("--correction_method", default='fdr_bh', help="Method for multiple testing correction (e.g., 'bonferroni', 'fdr_bh')")

    args = parser.parse_args()

    # Command handling
    if args.command == "create_pfm_matrices":
        create_pfm_matrices(args.input_dir, args.output_dir, threshold=args.filter_mean_threshold)

    elif args.command == "calculate_ic":
        calculate_ic(args.input_dir, args.output_dir, k=args.kmer_length)

    elif args.command == "ic_analysis":
        ic_plots(args.input_dir, args.output_dir, k=args.kmer_length)

    elif args.command == "calculate_kl":
        if args.ic_threshold is None:
            ic_threshold = compute_default_ic_threshold(args.kmer_length)
        else:
            ic_threshold = args.ic_threshold
        calculate_kl_distance(args.input_dir, args.pfm_folder, args.output_dir, hocomoco_models_path=args.hocomoco_path, ic_threshold=ic_threshold)

    elif args.command == "kl_analysis":
        # Validation for family_names when family_files are used
        if args.family_files and not args.family_names:
            parser.error("Argument --family_names is required when --family_files is specified")
        perform_kl_analysis(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            tf_family_file=args.tf_family_file,
            family_files=args.family_files,
            family_names=args.family_names,
            class_descriptors=args.class_descriptors,
            test_all=args.test_all,
            test_level=args.test_level,
            correction_method=args.correction_method
        )

    elif args.command == "create_new_run":
        # Step 1: Determine run_id
        if args.run_name:
            run_id = args.run_name
        else:
            run_id = generate_run_id()
        print(f"Using run_id: {run_id}")

        # Step 2: Create the main run directory and its subdirectories
        run_dir, pfm_dir, ic_dir, kl_dir, ic_plots_dir, kl_plots_dir = create_directories_for_run(run_id)
        print(f"Created run directory: {run_dir}")

        # Step 3: Determine ic_threshold
        if args.ic_threshold is None:
            ic_threshold = compute_default_ic_threshold(args.kmer_length)
            print(f"No IC threshold provided. Using default IC threshold: {ic_threshold}")
        elif args.ic_threshold.lower() == 'manual':
            ic_threshold = None  # Will prompt user later
        else:
            try:
                ic_threshold = float(args.ic_threshold)
                print(f"Using user-provided IC threshold: {ic_threshold}")
            except ValueError:
                parser.error("Invalid value for --ic_threshold. Must be a number or 'manual'.")

        # (a) Create PFM matrices and save to pfm_dir
        create_pfm_matrices(args.input_dir, pfm_dir, threshold=args.filter_mean_threshold)

        # (b) Calculate Information Content using the output of PFM as input, saving to ic_dir
        calculate_ic(pfm_dir, ic_dir, k=args.kmer_length)

        # (c) Perform IC analysis using the IC directory, saving to ic_plots_dir
        ic_plots(ic_dir, ic_plots_dir, k=args.kmer_length)

        # If ic_threshold is 'manual', prompt the user after IC plots are generated
        if ic_threshold is None:
            print(f"IC plots have been generated and saved to {ic_plots_dir}.")
            print("Please review the plots to determine an appropriate IC threshold.")
            user_input = input("Enter the IC threshold value: ")
            try:
                ic_threshold = float(user_input)
                print(f"Using user-provided IC threshold: {ic_threshold}")
            except ValueError:
                print("Invalid input. IC threshold must be a numeric value.")
                return  # Exit the script or handle the error as appropriate

        # (d) Calculate KL divergence using the IC CSVs and PFM folder, saving to kl_dir
        calculate_kl_distance(ic_dir, pfm_dir, kl_dir, hocomoco_models_path=args.hocomoco_path, ic_threshold=ic_threshold)

        # Validation for family_names when family_files are used
        if args.family_files and not args.family_names:
            parser.error("Argument --family_names is required when --family_files is specified")

        # (e) Perform KL analysis using the KL directory, saving to kl_plots_dir
        perform_kl_analysis(
            input_dir=kl_dir,
            output_dir=kl_plots_dir,
            tf_family_file=args.tf_family_file,
            family_files=args.family_files,
            family_names=args.family_names,
            class_descriptors=args.class_descriptors,
            test_all=args.test_all,
            test_level=args.test_level,
            correction_method=args.correction_method
        )

        print(f"Pipeline completed successfully for run_id: {run_id}")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
