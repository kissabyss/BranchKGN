# Branch-Specific Gene Discovery in Cell Differentiation Using Multi-Omics Graph Attention

This repository is suitable for the article "Branch-Specific Gene Discovery in Cell Differentiation Using Multi-Omics Graph Attention". The script contains all the code and raw data involved in the paper.

## Requirements

### OS Requirements

```{bash}
Python Version: 3.8.18
R Version:4.4.1
CUDA Version: 12.1
GPU: 2x NVIDIA RTX A4000
Operating System: Windows 10 Professional Workstation Edition
```

### Python Dependencies

All required Python packages are listed in `requirements.txt`. You can install them using:

```bash
pip install -r requirements.txt
```

### R Dependencies

All required R packages are listed in `gas_install.R`. You can install them using:

```R
source("gas_install.R")
```

### Manual installation (Optional)

If you encounter issues with the automatic installation, you can manually install the following packages:

```bash
torch                     2.1.0+cu121              
torch-cluster             1.6.2+pt21cu121          
torch-geometric           2.1.0.post1              
torch-scatter             2.1.2+pt21cu121          
torch-sparse              0.6.18+pt21cu121          
```

## Installation and Usage

### Download the Repository

1. Clone the repository from GitHub:

   ```bash
   git clone https://github.com/kissabyss/BranchKGN.git
   cd BranchKGN
   ```

2. Create a virtual environment using Anaconda:

   ```bash
   conda create -n branchK python=3.8.18
   conda activate branchK
   ```

3. Install the Python dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. Install the R dependencies:

   ```R
   source("gas_install.R")
   ```


### Set Up the Environment

Ensure that your CUDA version matches the requirements and that your GPU drivers are up to date. You may need to adjust the CUDA version in the `requirements.txt` file if necessary.

### Data sources

The data used in this repository is sourced from 10x Genomics. You can download the datasets from the following links:

- [PBMC from a Healthy Donor - Granulocytes Removed (3k)](https://www.10xgenomics.com/cn/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0)
- [PBMC from a Healthy Donor - No Cell Sorting (3k)](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0)
- [PBMC from a Healthy Donor - No Cell Sorting (10k)](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-2-0-0)

### File Description

This repository contains the following files and folders:

- **Data**: Contains the raw and processed data files.
- **pyHGT**: Python scripts for building a Graph Attention Transformer for bioinformatics tasks.
- **Result**: Contains the output files, including charts, reports, and analysis results.
- **arg.py**: Python script defining command-line arguments.
- **gas.R**: R script for data analysis and processing.
- **gas_gene.R**: R script focused on gene-related analysis.
- **gas_install.R**: R script to install required packages.
- **gas_simple.R**: Simplified R script for basic data analysis.

## Example Workflow

To run the analysis using the `pbmc_unsorted_3k` dataset, follow these steps:

1. Ensure the dataset is downloaded and placed in the `Data` folder.
2. Open `gas_simple.R` and run the script line by line.
3. Verify the dataset input location and format to ensure compatibility.

## License

This repository is licensed under the [MIT License](LICENSE).

