# Mutational antigenic profiling of SARS-CoV-2 RBD against HAARVI cohort sera
Analysis of mutational antigenic profiling of barcoded codon variants of SARS-CoV-2 RBD against convalescent patient sera from the HAARVI cohort.

Study and analysis by Allie Greaney, Andrea Loes, Kate Crawford and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html) in collaboration with the [Helen Chu lab](https://www.chulab.org/). The experiments are described in this [pre-print](https://doi.org/10.1101/2020.12.31.425021).

## Summary of workflow and results
For a summary of the workflow and links to key results files, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis
The analysis consists of three components, all of which are contained in this repository:

 1. Instructions to build the computing environment.

 2. The required input data.

 3. The computer code and a [Snakemake](https://snakemake.readthedocs.io) file to run it.

### Configure `.git` to not track Jupyter notebook metadata
To simplify git tracking of Jupyter notebooks, we have added the filter described [here](https://stackoverflow.com/questions/28908319/how-to-clear-an-ipython-notebooks-output-in-all-cells-from-the-linux-terminal/58004619#58004619) to strip notebook metadata to [.gitattributes](.gitattributes) and [.gitconfig](.gitconfig).
The **first time** you check out this repo, run the following command to use this configuration (see [here](https://stackoverflow.com/a/18330114)):

    git config --local include.path ../.gitconfig

Then don't worry about it anymore.

### Build the computing environment
First, set up the computing environment, which is partially done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The fully pinned environment is specified in [environment.yml](environment.yml), and an unpinned version is in [environment_unpinned.yml](environment_unpinned.yml).
If the environment already exists, you can activate it with:

    conda activate SARS-CoV-2-RBD_MAP

If you need to build the environment, then first build it with:

    conda env create -f environment.yml

Then activate it as above.

### Input data
The input data are specified in [./data/](data); see the README in that subdirectory for more details.

### Running the code
The analysis consists of Jupyter notebooks in the top-level directory along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile), specifying the conda environment, as in:

    snakemake --use-conda --conda-prefix /fh/fast/bloom_j/software/miniconda3/envs/SARS-CoV-2-RBD_MAP -R make_summary -j 1

However, you probably want to using a cluster to help with computationally intensive parts of the analysis.
To run using the cluster configuration for the Fred Hutch server, simply run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash), which executes [Snakefile](Snakefile) in a way that takes advantage of the Hutch server resources.
You likely want to submit [run_Hutch_cluster.bash](run_Hutch_cluster.bash) itself to the cluster (since it takes a while to run) with:

    sbatch -t 7-0 run_Hutch_cluster.bash

## Configuring the analysis
The configuration for the analysis is specifed in [config.yaml](config.yaml).
This file defines key variables for the analysis, and should be relatively self-explanatory.
You should modify the analysis by changing this configuration file; do **not** hard-code crucial experiment-specific variables within the notebooks or `Snakefile`.

In general:
 - add new samples to [data/barcode_runs.csv](data/barcode_runs.csv)
 - specify new combinations of samples for escape-profile plotting via [data/escape_profiles_config.yaml](data/escape_profiles_config.yaml)
 - specify combinations of samles for multi-dimensional scaling via [data/mds_config.yaml](mds_config.yaml)
 - to output structural mappings, add information to [data/output_pdbs_config.yaml](data/output_pdbs_config.yaml) about which conditions should be mapped to which PDBs, and add information to [data/structural_annotation_config.yaml](data/structural_annotation_config.yaml) about which chains to analyze for defining structural contacts between RBD and antibody or other ligands.
 - to write supplementary data and `dms_view` input data, set the `make_supp_data` flag in [data/escape_profiles_config.yaml](data/escape_profiles_config.yaml) to `true` and `dms_view` input files will be written to [results/supp_data](results/supp_data) for the condition sets with `make_supp_data` as `true` for those conditions with PDB mappings in [data/output_pdbs_config.yaml](data/output_pdbs_config.yaml).
 - to update the GISAID sequence set used to look at natural mutations, update the file pointed to by `gisaid_spikes` in [config.yaml](config.yaml) as described in [data/README.md](data/README.md).
 - to plot information about viral escape selections, add them to [data/escape_selection_results.yaml](data/escape_selection_results.yaml).

## Cluster configuration
There is a cluster configuration file [cluster.yaml](cluster.yaml) that configures [Snakefile](Snakefile) for the Fred Hutch cluster.
The [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script uses this configuration to run [Snakefile](Snakefile).
If you are using a different cluster than the Fred Hutch one, you wll need to modify the cluster configuration file.

## Notebooks that perform the analysis
The Jupyter notebooks that perform most of the analysis are in this top-level directory with the extension `*.ipynb`.
These notebooks read the key configuration values from [config.yaml](config.yaml).

There is also a [./scripts/](scripts) subdirectory with related scripts.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them via [Snakefile](Snakefile) as described above.

## Results
Results are placed in the [./results/](results) subdirectory.
Many of the files created in this subdirectory are not tracked in the `git` repo as they are very large.
However, key results files are tracked as well as a summary that shows the code and results.
Click [here](./results/summary/summary.md) to see that summary.

The large results files are tracked via [git-lfs](https://git-lfs.github.com/).
This requires `git-lfs` to be installed, which it is in the `conda` environment specified by [environment.yml](environment.yml).
The following commands were then run:

    git lfs install

You may need to run this if you are tracking these files and haven't installed `git-lfs` in your user account.
Then the large results files were added for tracking with:

    git lfs track <FILENAME>

## Updating the conda environment
[environment.yml](environment.yml) contains a fully pinned conda environment.
An environment without all of the versions pinned is in [environment_unpinned.yml](environment_unpinned.yml).
If you need to update the environment, the suggested way to do it is add the new requirement to [environment_unpinned.yml](environment_unpinned.yml), then build, activate, and export that environment.
The last three commands can be done by the following commands:

    conda env create -f environment_unpinned.yml
    conda activate SARS-CoV-2-RBD_MAP
    conda env export > environment.yml

## Creating "subset" repos and uploading data to the SRA
Currently this repo contains analyses of many antibodies and sera, and should remain public since collaborators do not want all of these data to be public.

For papers, you can make a public "subset" repo by following the instructions in [./subset_data/](subset_data).
After making a subset repo, you can upload sequencing data to the Sequence Read Archive (SRA) following the instructions in [./SRA_upload/](SRA_upload).
