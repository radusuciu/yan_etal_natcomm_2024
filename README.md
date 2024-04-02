# Inflammation causes insulin resistance via interferon regulatory factor 3 (IRF3)-2 mediated reduction in FAHFA levels

This repository contains code used for post-processing of proteomics data associated with the publication submitted by Yan et al to Nature Communications in 2024. The proteomics experiments looked at the selectivity of an inhibitor ABD-110000 in a competitive format using a serine hydrolase selective probe (FP-biotin) at a range of concentrations in mouse brain and kidney membranes (profiled by reductive dimethylation labeling at two concentrations), and white adipose tissue (WAT; profiled by TMT18 labeling at 5 concentrations). PRM experiments were also performed to measure in vivo target engagement, however, these were analyzed with Skyline and no custom code was used in the processing of that data. Please refer to the paper for full methods, or feel free to open an issue if you have any questions or are seeking help reproducing these results.

## Usage

All scripts in this repository were run in a Ubuntu Linux 22.04 environment with Python 3.10.12 installed. The following commands were issued to setup and activate the environment, as well as to install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

To re-generate the data that was plotted in Supplementary Figure 7, please read the below sections that described the use of the two scripts included in this reposistory: `tmt.py` and `redime.py`.

### Reductive dimethylation datasets

Raw files were searched and quantified using ProteomeDiscoverer (PD). Supplementary Figures 7e and 7f were generated in Graphpad Prism using the tables generated in the script `redime.py`. This script requires that all `*_Proteins.txt` files output by PD are in the `data/input` directory, relative to `redime.py`. These data files must be obtained from PRIDE (see below section). The script extracts ratios corresponding to serine hydrolases from each input file, converts them into percent competition values, and outputs a CSV report for each tissue into the `data/output` directory. The script also normalizes the ratio for each file so that the median ratio for all proteins is equal to 1.

To execute this script, ensure that the environment has been setup as described below and then issue the command `python redime.py`.

### TMT datasets

Raw files were searched using PD, and quantification was carried out using the code included in the script `tmt.py` using the PSM level export from PD, and was plotted in Supplementary Figure 7g using Graphpad Prism. This script requires that all `*_PSMs.txt` files output by PD are in the `data/input` directory, relative to `tmt.py`. These data files must be obtained from PRIDE (see below section). The script sums PSM level signal-to-noise (SNR) values for each protein, converts them into percent competition values relative to the mean of the protein-level SNR, and outputs a CSV report into `data/output` directory. The script applies two kinds of normalization: 

1. Total channel SNR normalization (by condition): The SNR of each PSM is adjusted so that each individual channel within that condition (a replicate) has a total SNR equal to the mean across all channels belonging to that condition.
2. Median normalization: The competition ratios of each protein in an individual channel (replicate), are normalized so that the median ratio for all proteins is equal to 1.

To execute this script, ensure that the environment has been setup as described below and then issue the command `python tmt.py`.

## Methods

For all proteomics methods, including versions of all software used, please refer to the publication (URL to be added upon publication).

## Data Files

All data files have been uploaded to PRIDE at this URL: https://www.ebi.ac.uk/pride/archive/projects/PXD047864/
