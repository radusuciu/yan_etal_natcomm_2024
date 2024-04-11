# Inflammation causes insulin resistance via interferon regulatory factor 3 (IRF3)-2 mediated reduction in FAHFA levels

This repository contains code used for the processecing of the proteomics data associated with the publication submitted by Yan et al. to Nature Communications in 2024. The proteomics experiments aimed to characterize the selectivity of the inhibitor ABD-110000 in a competitive format using a serine hydrolase selective probe (FP-biotin) at a range of concentrations in mouse brain and kidney membranes (profiled by reductive dimethylation (ReDiMe) labeling at two concentrations), and white adipose tissue (WAT; profiled by TMTPro18 labeling at 5 concentrations).

PRM experiments were also performed to measure in vivo target engagement, however, these were analyzed with Skyline and no custom code was used in the processing of that data. Please refer to the paper for full methods, or feel free to open an issue if you have any questions or are seeking help reproducing these results.

The code provided here includes the custom pipeline used to process both the ReDiMe and TMT data, as well as code used to generate tabular reports that served as the source for the proteomics plots in Supplementary Figure 7. The following sections will describe the usage and function of all code contained in this repository. Please read the paper and methods for more context and additional details.

## Usage

All scripts in this repository were run in a Ubuntu Linux 22.04 environment with Python 3.11.9 installed using [`pyenv`](https://github.com/pyenv/pyenv). Python 3.11 is required primarily for [`snakemake`](https://snakemake.github.io/) which is used for pipeline orchestration. Assuming Python is installed and that `python` points to a Python 3.11 binary, you can issue the following commands to setup and activate the environment, as well as to install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**NOTE**: Unless otherwise noted, all future commands presented in this document assume that the environment setup above is activated and that all dependencies are present.

The one other other requirement is having [`docker`](https://github.com/docker) which is used to distribute an image that contains all software used in the pipeline. This image is defined in the [`Dockerfile`](Dockerfile) and has been uploaded to the GitHub container registry at this URL: https://ghcr.io/radusuciu/yan_etal_natcomm_2024. Docker version 24.0.2 was used, but any recent version will likely work. The image will be automatically downloaded if you execute any of the pipelines, but you can also download it separately like so: `docker pull ghcr.io/radusuciu/docker-percolator`.

Finally, though the code here was run on Ubuntu 22.04, any Linux operating system where you can install docker should work. You should be able to use Windows as well, though the commands for activating the virtual environment need to be modified accordingly. I would recommend that you [`use WSL`](https://learn.microsoft.com/en-us/windows/wsl/install) if you're on Windows.

### Data Files

All data files have been uploaded to PRIDE at this URL: https://www.ebi.ac.uk/pride/archive/projects/PXD047864/. For convenience, a script, `download.py` is provided which will download all raw files necessary for re-running the pipeline and generating reports, using [`ppx`](https://github.com/wfondrie/ppx). To use this script, issue the command:

```bash
python download.py
```

This may take a while depending on your connection! After running this script, all of the files required for running both TMT and ReDiMe pipelines, and for generating the respective reports, will be placed in the `data/input` folder.

**NOTE**: The PRIDE dataset also includes data from a dose response TMT experiment in mouse brain, but these data were not referenced in the final final submission, and are thus not included in the analyses here. Uou may of course modify the code to take a look at that data as well -- please feel free to open an issue if you need assistance with this.

### ReDiMe

ReDiMe raw files were converted using [`ThermoRawFileParser`](https://github.com/compomics/ThermoRawFileParser), searched with [`comet`](https://uwpr.github.io/Comet/), processed with a number of tools from the [`OpenMS`](https://github.com/openms/openms) project, post-processed with [`percolator`](https://github.com/percolator/percolator), and quantified using [`Cimage`](https://github.com/wangchulab/CIMAGE).

#### Pipeline

Assuming all input files are in `data/input` (see [data files](#data-files) section if not), you can run the pipeline by issuing the following command:

```bash
# edit the number of --cores if you want to go faster!
snakemake --snakefile redime.smk --cores 1
```

Please inspect [`redime.smk`](redime.smk) to see exact parameters used for each step of the workflow. After successful execution of the workflow, there will be a file named after each sample, with an `.raw_output_rt_10_sn_2.5.to_excel` extension in the `data/input` folder. These files serve as input for report generation.

#### Report Generation

To generate the reports used for the corresponding figures, ensure that the pipeline has been run and then issue the following command:

```bash
python redime.py
```

The script extracts ratios corresponding to serine hydrolases from each input file, converts them into percent competition values, and outputs a CSV report for each tissue into the `data/output` directory. The script also normalizes the ratio for each file so that the median ratio for all proteins is equal to 1.

### TMT

TMT raw files were processed as described in the [ReDiMe section](#redime), with the exception of quantification, which performed using code in [`tmt.py`](tmt.py) and described below.

#### Pipeline

Assuming all input files are in `data/input` (see [data files](#data-files) section if not), you can run the pipeline by issuing the following command:

```bash
# edit the number of --cores if you want to go faster!
snakemake --snakefile tmt.smk --cores 1
```

Please inspect [`tmt.smk`](tmt.smk) to see exact parameters used for each step of the workflow. After successful execution of the workflow, the following two files will be in the `data/input` folder: `wat_dose_response_filtered_matches.idXML` and `wat_dose_response.mzML`. These files serve as input for report generation.

#### Report Generation

To generate the reports used for the corresponding figures, ensure that the pipeline has been run and then issue the following command:

```bash
python tmt.py
```

The script sums PSM level signal-to-noise (SNR) values for each protein, converts them into percent competition values relative to the mean of the protein-level SNR, and outputs a CSV report into `data/output` directory. Two kinds of normalization are applied: 

1. Total channel SNR normalization (by condition): The SNR of each PSM is adjusted so that each individual channel within that condition (a replicate) has a total SNR equal to the mean across all channels belonging to that condition.
2. Median normalization: The competition ratios of each protein in an individual channel (replicate), are normalized so that the median ratio for all proteins is equal to 1.

## General Notes

- If you have any questions about the code here, or have issues replicating the analyses, please feel free to open and issue.
- The code used for report generation has been extracted and modified from an internal Lundbeck application. For this reason it sometimes might do a bit more than would be strictly necessary for this publication.
- In both the TMT and ReDiMe workflows, there is a rule named `add_mod_descriptions`. The reason this exists is because I'm using `comet` instead of `CometAdapter` and the `pep.xml` output by `comet` needs a few tweaks so that the differential and static modifications can be properly translated to the OpenMS `idXML` format by `IDConvert` -- see [`utils.py`](utils.py#add_mod_descriptions_to_comet_output) for implementation details. The reason that I chose to do this is because `CometAdapter` had performance issues when combined with older versions of `comet`. These are solved now, so if you re-process with a newer version of `comet` you can likely remove this step.
- I've replaced the `unimod.xml` file that ships with OpenMS to add a differential modifications corresponding to the difference between the dimethyl modification (specified as a static mod), and the heavy modification used. The modification is named Dimethyl Heavy Delta and has a composition 13C(2) 2H(4) C(-2) H(-4) and monoisotopic mass of 6.031817.

## TODO

- add link to paper once published
- customize resources used by each rule
- review code and readme, adding comments/corrections as necessary

## Credits

Some of the code included here, particularly the bits related to processing Cimage data were originally authored by Kenneth Lum.
