# How to make a new repo with just a subset of data
The master repository at [https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP](https://github.com/jbloomlab/SARS-CoV-2-RBD_MAP) for the SARS-CoV-2 mutational antigenic profiling is currently a private repository for the Bloom lab that contains data for lots of antibodies and sera.
Some of these data need to be kept private because they relate to antibodies or sera that have not yet been published and which collaborators may not want made public at this time.

Therefore, for specific papers you may want to make public versions of the repo that only contain a subset of the data for specific antibodies or sera.
Importantly, for these public subset repos you also need to clear the GitHub history relative the master repo, because otherwise someone could back through the history and see all the private data.

Here are the steps to do this:

 1. Edit the [samples_to_subset.csv](samples_to_subset.csv) file in this directory to contain just the names of the samples that you want to include in the new public subset-repo. These names should be specified in the same way as in the [../data/escape_profiles_config.yaml](../data/escape_profiles_config.yaml) file.

 2. Then run the Jupyter notebook [subset_barcode_runs.ipynb](subset_barcode_runs.ipynb), which will output a list of all of the Illumina barcode runs that are needed for those samples. This list of Illumina barcode runs will be in the created file [barcode_runs_subset.csv](barcode_runs_subset.csv). These runs will be a subset of the collection of all runs in [../data/barcode_runs.csv](../data/barcode_runs.csv).

 3. Now copy the **entire** repo to a new subdirectory with the name that you want to give the subset repo. For instance, if you want the new repo to be called `SARS-CoV2-RBD_MAP_Crowe_antibodies`, then you would do:

        cp -r SARS-CoV-2-RBD_MAP SARS-CoV2-RBD_MAP_Crowe_antibodies

 4. Now you want to copy the `barcode_runs_subset.csv` file that you have created in step (2) above to replace the `./data/barcode_runs.csv` file in your new subset repo subdirectory (e.g., `SARS-CoV2-RBD_MAP_Crowe_antibodies` in example in (3)). This means that for the subset repo in that subdirectory, you now have a `./data/barcode_runs.csv` file that only contains the barcode runs for the antibodies / sera of interest.

 5. You also want to remove extraneous antibody / sera sets that aren't relevant to this subset from [../data/escape_profiles_config.yaml](../data/escape_profiles_config.yaml).

 6. Next, you need to completely wipe the results and git history in the new subset repo. You do this by navigating to the subset repo ((e.g., `SARS-CoV2-RBD_MAP_Crowe_antibodies` in example in (3)) and doing:

        rm -rf .git
        rm -rf results

 7. Currently the general pipeline isn't totally general, and one of the steps (the `evolution_escape_XXX.Rmd` R markdown notebooks) are different for each antibody subset. So you need to delete the irrelevant ones of these and also remove them from the `Snakefile`. It also appears that the [../structural_views/](../structural_views) subdirectory isn't modular and so may need manual editing.

 9. Run the full `Snakemake` pipeline in the subset repo subdirectory.

 10. Now initialize a new GitHub repo in the subset directory you've made, commit it with a message something like "initial public version of repo for Crowe antibodies" and push it to GitHub.
