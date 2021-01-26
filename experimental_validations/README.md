# RBD depletions of HAARVI sera to determine how much of each sample's neutralization potency is derived from RBD-targeting antibodies.


## Experimental procedure
I used commercially-available magnetic beads conjugated to SARS-CoV-2 RBD to deplete serum of RBD-binding antibodies. 
Control pre-pandemic serum spiked with either an RBD-binding monoclonal antibody, rREGN10987 (at 5, 25, or 50 ug/mL) or an NTD-binding antibody (at 50 ug/mL) were also incubated with the RBD beads to get an estimate of how well the depletion works. 

The process of pulling down the RBD antibodies involved diluting the serum 1:4 (because I added 50 uL serum + 150 uL of bead suspension at 1mg/mL). 
In one condition that we will not plot here, I used 2x as many beads, and thus added 50 uL serum + 300 uL of bead suspension, so this was a 1:7 dilution of the initial serum.

The pre-depletion serum was also diluted 1:4 (or 1:7, as appropriate) in PBS + 0.05% BSA (the buffer the beads were suspended in). 

Depletions were performed overnight at 4C in Eppendorf tubes in a tube rotator. 

## ELISAs to verify degree of RBD depletion
SARS-CoV-2 RBD and spike ELISAs were performed on these samples to verify the degree of RBD antibody depletion.

[Here](./rbd_depletions.ipynb), I calculate the AUC for pre- vs. post-depletion samples and the fold-change in AUC. 

This is really just a test of whether the depletion worked at all and I don't think it is very useful for anything beyond that (see notes at the end!).

## Neutralization assays to determine degree of neutralization derived from RBD antibodies
Andrea Loes then performed neutralization assays on these samples. 
She started each sample at a slightly different initial dilution to capture the full neutralization curve for each sample, given that each serum has a different neutralization potency (determined in Kate's original paper).

The data that are analyzed in [this notebook](./rbd_depletions.ipynb) were pre-processed by Kate's `excel_to_fracinfect` script. 

## Mutant virus neutralization assays
These pseudovirus neutralization assays were performed by me (Allie Greaney) and Andrea Loes to test how well the escape maps predict which mutations can escape neutralization by the sera. 

The data that are analyzed in [this notebook](./mutant_neuts.ipynb) were pre-processed by Kate's `excel_to_fracinfect` script.

## Mutant virus titering and p24 ELISAs
The pseudoviruses produced for this study (excluding Y369N and Y396T mutants, which did not grow well in initial titerng) were titered in 293T-ACE2 cells by Kate Crawford on January 25, 2021. p24 ELISAs were conducted according to manufacturer instructions [found here](http://ablinc.com/assets/HIVp24ELISA.pdf) on January 22, 2021 using a 1:100,000 dilution of lentivirual supernatants. 

Raw and analyzed plate reader data are in the `21-01-25....xlsx` and `21-01-22....xlsx` files in the [./data](./data) subdirectory. 

Plate reader results were manually analyzed and organized into the [2102125_p24_titers_calc.csv](./data/210125_p24_titers_calc.csv) file. These results are plotted in the [mutant_titers](./mutant_titers.ipynb) notebook and resulting plots are saved [here](./results/p24_titering).
