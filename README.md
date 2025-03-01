# CNN Filter Analysis Package

This tool provides a command-line interface (CLI) for performing various operations such as generating Position Frequency Matrices (PFM), calculating Information Content (IC), conducting Kullback-Leibler (KL) distance analysis, and visualizing results. The CLI includes several subcommands to automate these tasks.

## Table of Contents
- [Installation](#installation)
- [Commands](#commands)
  - [`create_new_run`](#create_new_run)
  - [`create_pfm_matrices`](#create_pfm_matrices)
  - [`calculate_ic`](#calculate_ic)
  - [`ic_analysis`](#ic_analysis)
  - [`calculate_kl`](#calculate_kl)
  - [`kl_analysis`](#kl_analysis)
- [Examples](#examples)
  - [Sample Usage](#sample-usage)
    - [Hypothesis Testing on All Classes](#hypothesis-testing-on-all-classes)
    - [Testing Specific Class Descriptors](#testing-specific-class-descriptors)
- [Contact](#contact)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/jmstepka/cnn-filter-analysis
   ```

2. Install the required libraries. The package was developed using Python version **3.10.14**:
   ```bash
   pip install -r requirements.txt
   ```

## Commands

### `create_new_run`

Generates a new run and executes all data processing steps, including generating PFM matrices, IC calculations, KL calculations, and analyses. You can optionally name your run using the `--run_name` option; otherwise, a unique `run_id` will be generated.

#### Usage:
```bash
python cli.py create_new_run \
  --input_dir <input_directory> \
  -t <ic_threshold> \
  -k <kmer_length> \
  --hocomoco_path <hocomoco_model_path> \
  --tf_family_file <tf_family_csv> \
  [ -c <class_descriptors> ] \
  [ --test_all ] \
  [ --test_level <level> ] \
  [ --correction_method <method> ] \
  --filter_mean_threshold <threshold> \
  [ --run_name <run_name> ]
```

- **`--input_dir`**: Input directory containing the data to process.
- **`-t` / `--ic_threshold`**: IC threshold for analysis. Can be a numeric value, or `'manual'` to input after IC plots are generated. Default is `(k + 1) / 2`.
- **`-k` / `--kmer_length`**: Length of k-mers used for analysis (default is 7).
- **`--hocomoco_path`**: Path to HOCOMOCO models in MEME format.
- **`--tf_family_file`**: Path to the TF family CSV file.
- **`-c` / `--class_descriptors`**: (Optional) List of class descriptors to test. Each descriptor can be in `'x'`, `'x.y'`, or `'x.y.z'` format.
- **`--test_all`**: (Optional) If provided, test all superclasses or classes.
- **`--test_level`**: (Optional) Level to test when `--test_all` is used (`superclass` or `class`). Default is `superclass`.
- **`--correction_method`**: (Optional) Method for multiple testing correction (`bonferroni`, `fdr_bh`, etc.). Default is `fdr_bh`.
- **`--filter_mean_threshold`**: Mean threshold for filtering results (default is 1e-5).
- **`--run_name`**: (Optional) Name of the run. If not specified, a unique `run_id` will be generated (e.g., `run_1`, `run_2`, etc.).

#### Notes:

- **Naming Your Run:**
  - Use the `--run_name` option to specify a custom name for your run. This name will be used to create directories and organize output files.
  - If `--run_name` is not provided, the script will automatically generate a unique `run_id`.

- **IC Threshold:**
  - If you set `--ic_threshold manual`, the script will generate IC plots and prompt you to input an appropriate IC threshold after reviewing the plots.

- **Testing Options:**
  - You must specify one of the following:
    - `-c` / `--class_descriptors` to test specific class descriptors.
    - `--test_all` to test all superclasses or classes.
    - `-f` / `--family_files` to specify motif family files (requires `-n` / `--family_names`).

### `create_pfm_matrices`

Generates PFM matrices based on the input data. The results will be saved in the specified output directory.

#### Usage:
```bash
python cli.py create_pfm_matrices \
  --input_dir <input_directory> \
  --output_dir <output_directory> \
  --filter_mean_threshold <threshold>
```

- **`--input_dir`**: Input directory containing the data to process.
- **`--output_dir`**: Output directory where the PFM matrices will be saved.
- **`--filter_mean_threshold`**: Optional mean threshold for selecting filters (default is 1e-5).

### `calculate_ic`

Calculates Information Content (IC) from previously generated PFM matrices.

#### Usage:
```bash
python cli.py calculate_ic \
  --input_dir <input_directory> \
  --output_dir <output_directory> \
  -k <kmer_length>
```

- **`--input_dir`**: Input directory containing PFM matrices.
- **`--output_dir`**: Output directory where the IC results will be saved.
- **`-k` / `--kmer_length`**: Length of k-mers used for calculations (default is 7).

### `ic_analysis`

Performs IC analysis based on the results saved in the input directory.

#### Usage:
```bash
python cli.py ic_analysis \
  --input_dir <input_directory> \
  --output_dir <output_directory> \
  -k <kmer_length>
```

- **`--input_dir`**: Input directory containing IC data.
- **`--output_dir`**: Output directory where the analysis results will be saved.
- **`-k` / `--kmer_length`**: Length of k-mers used for analysis (default is 7).

### `calculate_kl`

Calculates the KL (Kullback-Leibler) distance between two PFM matrices based on IC values.

#### Usage:
```bash
python cli.py calculate_kl \
  --input_dir <input_directory> \
  --output_dir <output_directory> \
  -t <ic_threshold> \
  --pfm <pfm_folder> \
  --hocomoco_path <hocomoco_model_path>
```

- **`--input_dir`**: Input directory containing IC results.
- **`--output_dir`**: Output directory where the KL results will be saved.
- **`-t` / `--ic_threshold`**: IC threshold for filtering motifs (default is 6.5).
- **`--pfm`**: Directory containing PFM matrices for comparison.
- **`--hocomoco_path`**: Path to HOCOMOCO models in MEME format.

### `kl_analysis`

Performs analysis on the results of the KL distance calculations, including plotting and family enrichment analysis.

#### Usage:
```bash
python cli.py kl_analysis \
  --input_dir <input_directory> \
  --output_dir <output_directory> \
  --tf_family_file <tf_family_csv> \
  [ -c <class_descriptors> ] \
  [ --test_all ] \
  [ --test_level <level> ] \
  [ --correction_method <method> ]
```

- **`--input_dir`**: Input directory containing KL distance results.
- **`--output_dir`**: Output directory where the analysis results will be saved.
- **`--tf_family_file`**: Path to the TF family CSV file.
- **`-c` / `--class_descriptors`**: (Optional) List of class descriptors to test. Each descriptor can be in `'x'`, `'x.y'`, or `'x.y.z'` format.
- **`--test_all`**: (Optional) If provided, test all superclasses or classes.
- **`--test_level`**: (Optional) Level to test when `--test_all` is used (`superclass` or `class`). Default is `superclass`.
- **`--correction_method`**: (Optional) Method for multiple testing correction (`bonferroni`, `fdr_bh`, etc.). Default is `fdr_bh`.


## Examples

### Sample Usage

#### Hypothesis Testing on All Classes

To perform a full run, including KL analysis and hypothesis testing on all classes, using the `create_new_run` command and naming the run `my_full_run`:

```bash
python cli.py create_new_run \
  --input_dir ./data/input_example \
  --kmer_length 7 \
  --hocomoco_path data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
  --tf_family_file data/tf_families.csv \
  --test_all \
  --test_level class \
  --filter_mean_threshold 1e-5 \
```

#### Testing Specific Class Descriptors

Alternatively, to test specific class descriptors (e.g., `2.2` and `2.3`), and naming the run `my_specific_run`, you can use:

```bash
python cli.py create_new_run \
  --input_dir ./data/input_example \
  --kmer_length 7 \
  --hocomoco_path data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
  --tf_family_file data/tf_families.csv \
  -c 2.2 2.3 \
  --filter_mean_threshold 1e-5 \
```
