# CNN Filter Analysis Package

This tool provides a command-line interface (CLI) for performing various operations such as generating Position Frequency Matrices (PFM), calculating Information Content (IC), conducting Kullback-Leibler (KL) distance analysis, and visualizing results. The CLI includes several subcommands to automate these tasks.

## Table of Contents
- [Installation](#installation)
- [Commands](#commands)
  - [create_new_run](#create_new_run)
  - [create_pfm_matrices](#create_pfm_matrices)
  - [calculate_ic](#calculate_ic)
  - [ic_analysis](#ic_analysis)
  - [calculate_kl](#calculate_kl)
  - [kl_analysis](#kl_analysis)
- [Examples](#examples)
  - [Sample Usage](#sample-usage)
- [Contact](#contact)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/cnn-filter-analysis.git
   ```

2. Install the required libraries:
   ```bash
   pip install -r requirements.txt
   ```

## Commands

### `create_new_run`
Generates a new unique `run_id` and executes all data processing steps, including generating PFM matrices, IC calculations, KL calculations, and analyses.

#### Usage:
```bash
python cli.py create_new_run --input_dir <input_directory> -t <ic_threshold> -k <kmer_length> -f <family_files> -n <family_names> --hocomoco_path <hocomoco_model_path> --filter_mean_threshold <threshold>
```
- **`--input_dir`**: Input directory containing the data to process.
- **`-t` / `--ic_threshold`**: IC threshold for analysis (e.g., 6.5).
- **`-k` / `--kmer_length`**: Length of k-mers used for analysis (default is 7).
- **`-f` / `--family_files`**: List of files containing motif families.
- **`-n` / `--family_names`**: List of family names corresponding to the motif files.
- **`--hocomoco_path`**: Path to HOCOMOCO models in MEME format.
- **`--filter_mean_threshold`**: Mean threshold for filtering results (default is 1e-5).

### `create_pfm_matrices`
Generates PFM matrices based on the input data. The results will be saved in the specified output directory.

#### Usage:
```bash
python cli.py create_pfm_matrices --input_dir <input_directory> --output_dir <output_directory> --filter_mean_threshold <threshold>
```
- **`--input_dir`**: Input directory containing the data to process.
- **`--output_dir`**: Output directory where the PFM matrices will be saved.
- **`--filter_mean_threshold`**: Optional mean threshold for selecting filters (default is 1e-5).

### `calculate_ic`
Calculates Information Content (IC) from previously generated PFM matrices.

#### Usage:
```bash
python cli.py calculate_ic --input_dir <input_directory> --output_dir <output_directory> -k <kmer_length>
```
- **`--input_dir`**: Input directory containing PFM matrices.
- **`--output_dir`**: Output directory where the IC results will be saved.
- **`-k` / `--kmer_length`**: Length of k-mers used for calculations (default is 7).

### `ic_analysis`
Performs IC analysis based on the results saved in the input directory.

#### Usage:
```bash
python cli.py ic_analysis --input_dir <input_directory> --output_dir <output_directory>
```
- **`--input_dir`**: Input directory containing IC data.
- **`--output_dir`**: Output directory where the analysis results will be saved.

### `calculate_kl`
Calculates the KL (Kullback-Leibler) distance between two PFM matrices based on IC values.

#### Usage:
```bash
python cli.py calculate_kl --input_dir <input_directory> --output_dir <output_directory> -ic_thresh <ic_threshold> --pfm <pfm_folder> --hocomoco_path <hocomoco_model_path>
```
- **`--input_dir`**: Input directory containing IC results.
- **`--output_dir`**: Output directory where the KL results will be saved.
- **`-ic_thresh`**: IC threshold for filtering motifs (default is 6.5).
- **`--pfm`**: Directory containing PFM matrices for comparison.
- **`--hocomoco_path`**: Path to HOCOMOCO models in MEME format.

### `kl_analysis`
Performs analysis on the results of the KL distance calculations.

#### Usage:
```bash
python cli.py kl_analysis --input_dir <input_directory> --output_dir <output_directory> -f <family_files> -n <family_names>
```
- **`--input_dir`**: Input directory containing KL distance results.
- **`--output_dir`**: Output directory where the analysis results will be saved.
- **`-f` / `--family_files`**: List of files with motif families for comparison.
- **`-n` / `--family_names`**: List of family names corresponding to the motif files.

## Examples

### Sample Usage

   ```bash
   python cli.py create_new_run --input_dir ./data/input_example --ic_threshold 5 --kmer_length 7 -f data/motifs_families/HUMAN_mono_motifs_zinc.tsv data/motifs_families/HUMAN_mono_motifs_bhlh.tsv -n zinc bhlh --hocomoco_path data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme --filter_mean_threshold 1e-5
   ```