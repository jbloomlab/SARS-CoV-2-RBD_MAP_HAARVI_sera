# Custom analyses

This directory contains custom analyses for specific sets of antibodies or sera.

## Analyses in this directory:

#### Write pymol commmands to make PSEs for RBD bound to 6m0j
The notebook [pymol_commands_6moj.ipynb](./pymol_commands_6moj.ipynb) takes as input the pdb files output by `output_pdbs.ipynb`, the escape sites called by `call_strong_escape_sites.ipynb`, and the contact sites called by `annotate_structural_contacts.Rmd`.

It writes text files that contain commands that can be pasted into the PyMol command line to create a PSE that can be used for further structure-gazing and manual analysis.

This notebook might work for other antibodies, but I haven't tested it.
