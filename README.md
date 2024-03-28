# Inflammation causes insulin resistance via interferon regulatory factor 3 (IRF3)-2 mediated reduction in FAHFA levels

This repository contains code used for post-processing of proteomics data associated with the publication submitted by Yan et al to Nature Communications in 2024. The proteomics experiments looked at the selectivity of an inhibitor ABD-110000 in a competitive format using a serine hydrolase selective probe (FP-biotin) at a range of concentrations in mouse brain and kidney membranes (profiled by reductive dimethylation labeling at two concentrations), and white adipose tissue (WAT; profiled by TMT18 labeling at 5 concentrations). PRM experiments were also performed to measure in vivo target engagement, however, these were analyzed with Skyline and no custom code was used in the processing of that data. Please refer to the paper for full methods, or feel free to open an issue if you have any questions or are seeking help reproducing these results.

## Usage

All scripts in this repository were run in a Ubuntu Linux 22.04 environment with Python 3.10 installed. The following commands were issued to setup and activate the environment, as well as to install :

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Reductive dimethylation datasets

Raw files were searched and quantified using ProteomeDiscoverer. Supplementary Figures 7e and 7f were generated in Graphpad Prism using the tables generated in the script `redime.py`. This script requires that all `*_Proteins.txt` files output by ProteomeDiscoverer are in the `data/input` directory, relative to `redime.py`. These data files must be obtained from PRIDE (see below section). The script aggregates data across replicates for each protein at each tested concentration and converts ratios to percent inhibition values, which are plotted in the mentioned figures.

### TMT datasets


## Methods

For all proteomics methods, including versions of all software used, please refer to the publication (URL to be added upon publication).

## Data Files

All data files have been uploaded to PRIDE at this URL: https://www.ebi.ac.uk/pride/archive/projects/PXD047864/
