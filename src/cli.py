import argparse
from bioinformatics_tools import create_pfm_matrices, calculate_ic, ic_analysis, calculate_kl, kl_analysis

def get_input_directory(run_id=None, input_dir=None):
    """Helper function to resolve the input directory based on run_id or input_dir."""
    if run_id:
        # Map the run_id to the input directory, customize this logic if needed
        input_directory = f"/path/to/data/{run_id}"
    else:
        input_directory = input_dir
    return input_directory

def main():
    parser = argparse.ArgumentParser(description="Bioinformatics CLI for matrix operations and information theory analysis.")

    subparsers = parser.add_subparsers(dest="command", help="Sub-commands for different operations")

    # Input directory option (shared across all commands)
    def add_input_options(subparser):
        input_group = subparser.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--run_id", help="Run ID to point to the input directory")
        input_group.add_argument("--input_dir", help="Direct input directory")

    # Sub-command: create_pfm_matrices
    parser_create = subparsers.add_parser("create_pfm_matrices", help="Create PFM matrices from input")
    add_input_options(parser_create)
    parser_create.add_argument("-o", "--output", required=True, help="Output file to save the PFM matrices")
    
    # Sub-command: calculate_ic
    parser_ic = subparsers.add_parser("calculate_ic", help="Calculate Information Content from PFM matrices")
    add_input_options(parser_ic)
    parser_ic.add_argument("-o", "--output", required=True, help="Output file to save IC results")
    
    # Sub-command: ic_analysis
    parser_analysis = subparsers.add_parser("ic_analysis", help="Perform IC analysis with a threshold")
    add_input_options(parser_analysis)
    parser_analysis.add_argument("-t", "--threshold", type=float, required=True, help="Threshold for IC analysis")
    
    # Sub-command: calculate_kl
    parser_kl = subparsers.add_parser("calculate_kl", help="Calculate KL divergence between two PFM matrices")
    
    kl_input_group_1 = parser_kl.add_mutually_exclusive_group(required=True)
    kl_input_group_1.add_argument("--run_id_1", help="Run ID for the first PFM input directory")
    kl_input_group_1.add_argument("--input_dir_1", help="Direct input directory for the first PFM")

    kl_input_group_2 = parser_kl.add_mutually_exclusive_group(required=True)
    kl_input_group_2.add_argument("--run_id_2", help="Run ID for the second PFM input directory")
    kl_input_group_2.add_argument("--input_dir_2", help="Direct input directory for the second PFM")

    parser_kl.add_argument("-o", "--output", required=True, help="Output file to save KL divergence results")
    
    # Sub-command: kl_analysis
    parser_kl_analysis = subparsers.add_parser("kl_analysis", help="Perform analysis on KL divergence results")
    add_input_options(parser_kl_analysis)
    parser_kl_analysis.add_argument("-s", "--significance_level", type=float, required=True, help="Significance level for KL analysis")

    args = parser.parse_args()

    # Command handling
    if args.command == "create_pfm_matrices":
        input_dir = get_input_directory(args.run_id, args.input_dir)
        create_pfm_matrices(input_dir, args.output)
    elif args.command == "calculate_ic":
        input_dir = get_input_directory(args.run_id, args.input_dir)
        calculate_ic(input_dir, args.output)
    elif args.command == "ic_analysis":
        input_dir = get_input_directory(args.run_id, args.input_dir)
        ic_analysis(input_dir, args.threshold)
    elif args.command == "calculate_kl":
        input_dir_1 = get_input_directory(args.run_id_1, args.input_dir_1)
        input_dir_2 = get_input_directory(args.run_id_2, args.input_dir_2)
        calculate_kl(input_dir_1, input_dir_2, args.output)
    elif args.command == "kl_analysis":
        input_dir = get_input_directory(args.run_id, args.input_dir)
        kl_analysis(input_dir, args.significance_level)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
