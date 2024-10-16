import argparse
import os
import uuid
from create_pfm import create_pfm_matrices
from calculate_ic import calculate_ic
from ic_analysis import ic_plots
from calculate_kl import calculate_kl_distance
from kl_analysis import perform_kl_analysis


def generate_run_id():
    """Generate a new unique run_id."""
    return str(uuid.uuid4())  # Generates a unique run_id using UUID

def create_directories_for_run(run_id, base_dir="./data"):
    """Create the main directory for the run, and input/output directories."""
    run_dir = os.path.join(base_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)
    
    input_dir = os.path.join(run_dir, "input")
    os.makedirs(input_dir, exist_ok=True)
    
    return run_dir, input_dir

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
    parser_create.add_argument("-k", "--kmer_length", type=int, default=7, help="Length of k-mers to extract (default: 7)")

    # Sub-command: calculate_ic
    parser_ic = subparsers.add_parser("calculate_ic", help="Calculate Information Content from PFM matrices")
    add_input_output_options(parser_ic)
    parser_ic.add_argument("-k", "--kmer_length", type=int, default=7, help="Length of k-mers to use in IC calculation")

    # Sub-command: ic_analysis
    parser_analysis = subparsers.add_parser("ic_analysis", help="Perform IC analysis with a threshold")
    add_input_output_options(parser_analysis)
    parser_analysis.add_argument("-t", "--threshold", type=float, required=True, help="Threshold for IC analysis")

    # Sub-command: calculate_kl
    parser_kl = subparsers.add_parser("calculate_kl", help="Calculate KL divergence between two PFM matrices")
    add_input_output_options(parser_kl)
    parser_kl.add_argument("-ic_thresh", "--ic_threshold", type=float, default=6.5, help="IC threshold for filtering motifs")
    parser_kl.add_argument("-pfm", "--pfm_folder", required=True, help="Folder containing the PFM matrices")
    parser_kl.add_argument("--hocomoco_path", required=True, help="Path to HOCOMOCO models in meme format")

    # Sub-command: kl_analysis
    parser_kl_analysis = subparsers.add_parser("kl_analysis", help="Perform analysis on KL divergence results")
    add_input_output_options(parser_kl_analysis)
    parser_kl_analysis.add_argument("-f", "--family_files", nargs="+", required=True, help="List of motif family files")
    parser_kl_analysis.add_argument("-n", "--family_names", nargs="+", required=True, help="List of family names corresponding to family files")

    # New sub-command: create_new_run
    parser_new_run = subparsers.add_parser("create_new_run", help="Generate a new run_id and run all steps of the pipeline")
    parser_new_run.add_argument("-t", "--ic_threshold", type=float, required=True, help="Threshold for IC analysis")
    parser_new_run.add_argument("-k", "--kmer_length", type=int, default=7, help="Length of k-mers to use for analysis")
    parser_new_run.add_argument("-f", "--family_files", nargs="+", required=True, help="List of motif family files")
    parser_new_run.add_argument("-n", "--family_names", nargs="+", required=True, help="List of family names corresponding to family files")
    parser_new_run.add_argument("--hocomoco_path", required=True, help="Path to HOCOMOCO models in meme format")

    args = parser.parse_args()

    # Command handling
    if args.command == "create_pfm_matrices":
        create_pfm_matrices(args.input_dir, args.output_dir)

    elif args.command == "calculate_ic":
        calculate_ic(args.input_dir, args.output_dir, k=args.kmer_length)

    elif args.command == "ic_analysis":
        ic_plots(args.input_dir, args.output_dir, threshold=args.threshold)

    elif args.command == "calculate_kl":
        calculate_kl_distance(args.input_dir, args.pfm_folder, args.output_dir, hocomoco_models_path=args.hocomoco_path, ic_threshold=args.ic_threshold)

    elif args.command == "kl_analysis":
        perform_kl_analysis(args.input_dir, args.family_files, args.family_names, args.output_dir)

    elif args.command == "create_new_run":
        # Step 1: Generate a new run_id
        new_run_id = generate_run_id()
        print(f"Generated new run_id: {new_run_id}")

        # Step 2: Create the main run directory and subdirectories for input
        run_dir, input_dir = create_directories_for_run(new_run_id)
        print(f"Created input directory: {input_dir}")

        # Create subdirectories for each step of the pipeline
        pfm_dir = os.path.join(run_dir, "pfm")
        ic_dir = os.path.join(run_dir, "ic")
        ic_analysis_dir = os.path.join(run_dir, "ic_analysis")
        kl_dir = os.path.join(run_dir, "kl")
        kl_analysis_dir = os.path.join(run_dir, "kl_analysis")

        # (a) Create PFM matrices and save to pfm_dir
        create_pfm_matrices(input_dir, pfm_dir)

        # (b) Calculate Information Content using the output of PFM as input, saving to ic_dir
        calculate_ic(pfm_dir, ic_dir, k=args.kmer_length)

        # (c) Perform IC analysis using the IC directory, saving to ic_analysis_dir
        ic_plots(ic_dir, ic_analysis_dir, threshold=args.ic_threshold)

        # (d) Calculate KL divergence using the IC CSVs and PFM folder, saving to kl_dir
        calculate_kl_distance(ic_dir, pfm_dir, kl_dir, hocomoco_models_path=args.hocomoco_path, ic_threshold=args.ic_threshold)

        # (e) Perform KL analysis using the KL directory, saving to kl_analysis_dir
        perform_kl_analysis(kl_dir, args.family_files, args.family_names, kl_analysis_dir)

        print(f"Pipeline completed successfully for run_id: {new_run_id}")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
