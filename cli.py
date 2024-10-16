import argparse
import os
import uuid
from src.external.create_pfm import create_pfm_matrices
from src.external.calculate_ic import calculate_ic
from src.external.ic_analysis import ic_plots

def generate_run_id():
    """Generate a new unique run_id."""
    return str(uuid.uuid4())  # Generates a unique run_id using UUID

def create_directories_for_run(run_id, base_dir="/path/to/data"):
    """Create input and output directories for a new run."""
    input_dir = os.path.join(base_dir, run_id, "input")
    output_dir = os.path.join(base_dir, run_id, "output")
    
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    return input_dir, output_dir

def get_directories(run_id=None, input_dir=None, output_dir=None):
    """Helper function to resolve the input and output directories based on run_id or input_dir/output_dir."""
    if run_id:
        # Assuming run_id maps to some predefined input/output directories
        input_directory = f"/path/to/data/{run_id}/input"
        output_directory = f"/path/to/data/{run_id}/output"
    else:
        input_directory = input_dir
        output_directory = output_dir
    return input_directory, output_directory

def main():
    parser = argparse.ArgumentParser(description="Bioinformatics CLI for matrix operations and information theory analysis.")
    
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands for different operations")

    # Input and output options for each command
    def add_input_output_options(subparser):
        input_group = subparser.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--run_id", help="Run ID to point to predefined input and output directories")

        # Both input_dir and output_dir must be provided together
        input_dir_group = input_group.add_argument_group()
        input_dir_group.add_argument("--input_dir", help="Direct input directory")
        input_dir_group.add_argument("--output_dir", help="Direct output directory")

    # Sub-command: create_pfm_matrices
    parser_create = subparsers.add_parser("create_pfm_matrices", help="Create PFM matrices from input")
    add_input_output_options(parser_create)

    # Sub-command: calculate_ic
    parser_ic = subparsers.add_parser("calculate_ic", help="Calculate Information Content from PFM matrices")
    add_input_output_options(parser_ic)

    # Sub-command: ic_analysis
    parser_analysis = subparsers.add_parser("ic_analysis", help="Perform IC analysis with a threshold")
    add_input_output_options(parser_analysis)
    parser_analysis.add_argument("-t", "--threshold", type=float, required=True, help="Threshold for IC analysis")

    # Sub-command: calculate_kl
    parser_kl = subparsers.add_parser("calculate_kl", help="Calculate KL divergence between two PFM matrices")

    kl_input_group_1 = parser_kl.add_mutually_exclusive_group(required=True)
    kl_input_group_1.add_argument("--run_id_1", help="Run ID for the first PFM input")
    kl_input_group_1.add_argument("--input_dir_1", help="Direct input directory for the first PFM")
    kl_input_group_1.add_argument("--output_dir_1", help="Direct output directory for the first PFM")

    kl_input_group_2 = parser_kl.add_mutually_exclusive_group(required=True)
    kl_input_group_2.add_argument("--run_id_2", help="Run ID for the second PFM input")
    kl_input_group_2.add_argument("--input_dir_2", help="Direct input directory for the second PFM")
    kl_input_group_2.add_argument("--output_dir_2", help="Direct output directory for the second PFM")

    # Sub-command: kl_analysis
    parser_kl_analysis = subparsers.add_parser("kl_analysis", help="Perform analysis on KL divergence results")
    add_input_output_options(parser_kl_analysis)
    parser_kl_analysis.add_argument("-s", "--significance_level", type=float, required=True, help="Significance level for KL analysis")

    # New sub-command: create_new_run
    parser_new_run = subparsers.add_parser("create_new_run", help="Generate a new run_id and run all steps of the pipeline")
    parser_new_run.add_argument("-t", "--threshold", type=float, required=True, help="Threshold for IC analysis")
    parser_new_run.add_argument("-s", "--significance_level", type=float, required=True, help="Significance level for KL analysis")

    args = parser.parse_args()

    # Command handling
    if args.command == "create_pfm_matrices":
        input_dir, output_dir = get_directories(args.run_id, args.input_dir, args.output_dir)
        create_pfm_matrices(input_dir, output_dir)
    elif args.command == "calculate_ic":
        input_dir, output_dir = get_directories(args.run_id, args.input_dir, args.output_dir)
        calculate_ic(input_dir, output_dir)
    elif args.command == "ic_analysis":
        input_dir, output_dir = get_directories(args.run_id, args.input_dir, args.output_dir)
        ic_analysis(input_dir, args.threshold)
    elif args.command == "calculate_kl":
        input_dir_1, output_dir_1 = get_directories(args.run_id_1, args.input_dir_1, args.output_dir_1)
        input_dir_2, output_dir_2 = get_directories(args.run_id_2, args.input_dir_2, args.output_dir_2)
        calculate_kl(input_dir_1, input_dir_2, output_dir_1)
    elif args.command == "kl_analysis":
        input_dir, output_dir = get_directories(args.run_id, args.input_dir, args.output_dir)
        kl_analysis(input_dir, args.significance_level)
    
    elif args.command == "create_new_run":
        # Step 1: Generate a new run_id
        new_run_id = generate_run_id()
        print(f"Generated new run_id: {new_run_id}")

        # Step 2: Create input and output directories for this run_id
        input_dir, output_dir = create_directories_for_run(new_run_id)
        print(f"Created input directory: {input_dir}")
        print(f"Created output directory: {output_dir}")

        # Step 3: Run all steps of the pipeline
        print("Running the full pipeline...")

        # (a) Create PFM matrices
        create_pfm_matrices(input_dir, output_dir)
        # (b) Calculate Information Content
        calculate_ic(input_dir, output_dir)
        # (c) Perform IC analysis
        ic_analysis(input_dir, args.threshold)
        # (d) Calculate KL divergence (using the same input/output for simplicity)
        calculate_kl(input_dir, input_dir, output_dir)
        # (e) Perform KL analysis
        kl_analysis(input_dir, args.significance_level)
        
        print(f"Pipeline completed successfully for run_id: {new_run_id}")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
