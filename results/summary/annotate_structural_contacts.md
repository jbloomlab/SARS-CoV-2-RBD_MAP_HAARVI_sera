Annotate antibody structural contacts
================
Tyler Starr
10/11/2020

-   [Annotate sites on the RBD that contact
    antibody](#annotate-sites-on-the-rbd-that-contact-antibody)
    -   [Important note:](#important-note)
-   [Annotate sites on the antibody that contact the
    RBD](#annotate-sites-on-the-antibody-that-contact-the-rbd)
-   [IMPORTANT DISCLAIMER!!!!](#important-disclaimer)
-   [Save list of structural
    contacts](#save-list-of-structural-contacts)

This notebook analyzes antibody-bound RBD crystal and cryo-EM structures
to annotate structural contacts. It generates a csv listing
structurally-defined contact residues for each ligand (i.e. Ab or ACE2).

    require("knitr")
    knitr::opts_chunk$set(echo = T)
    knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

    #list of packages to install/load
    packages = c("yaml","bio3d","tidyverse")
    #install any packages not already installed
    installed_packages <- packages %in% rownames(installed.packages())
    if(any(installed_packages == F)){
      install.packages(packages[!installed_packages])
    }
    #load packages
    invisible(lapply(packages, library, character.only=T))

    #read in config file
    config <- read_yaml("config.yaml")

    #read in config file for determining which structure and chains to determine as contacts
    contacts_config <- read_yaml(file=config$structural_contacts_config)

    #make output directory
    if(!file.exists(config$structural_contacts_dir)){
      dir.create(file.path(config$structural_contacts_dir))
    }

Session info for reproducing environment:

    sessionInfo()

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] forcats_0.4.0   stringr_1.4.0   dplyr_0.8.3     purrr_0.3.3    
    ##  [5] readr_1.3.1     tidyr_1.0.0     tibble_3.0.1    ggplot2_3.3.0  
    ##  [9] tidyverse_1.3.0 bio3d_2.4-0     yaml_2.2.0      knitr_1.26     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.11        haven_2.2.0      lattice_0.20-38 
    ##  [5] colorspace_1.4-1 vctrs_0.2.4      generics_0.0.2   htmltools_0.4.0 
    ##  [9] rlang_0.4.5      pillar_1.4.3     glue_1.3.1       withr_2.1.2     
    ## [13] DBI_1.1.0        dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1    
    ## [17] lifecycle_0.2.0  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_0.3.5      evaluate_0.14    parallel_3.6.2   fansi_0.4.0     
    ## [25] broom_0.5.6      Rcpp_1.0.3       scales_1.1.0     backports_1.1.5 
    ## [29] jsonlite_1.6     fs_1.3.1         hms_0.5.2        digest_0.6.23   
    ## [33] stringi_1.4.3    grid_3.6.2       cli_2.0.0        tools_3.6.2     
    ## [37] magrittr_1.5     crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0  
    ## [41] xml2_1.2.2       reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1
    ## [45] rmarkdown_2.0    httr_1.4.1       rstudioapi_0.10  R6_2.4.1        
    ## [49] nlme_3.1-143     compiler_3.6.2

Annotate sites on the RBD that contact antibody
-----------------------------------------------

We use the `binding.site` function from the `bio3d` package to identify
residues that are structural complexes within each pdb. We iterate
through the values in `structural_annotation_config.yaml` to perform
this calculation for each structure described by that configuration
file.

### Important note:

Chains cannot be called `N` in the config file. This will lead to an
error (read as not a character). Change all `N` chains to `"N"`.

    #empty data frame to append contacts to
    structural_contacts <- data.frame(name=character(0), pdb=character(0), chain=character(0),position=numeric(0))
    antibody_contacts <- data.frame(name=character(0), pdb=character(0), chain=character(0), position=numeric(0))

    for(entry in contacts_config){
      pdb <- read.pdb(file=entry$pdbfile)
      # get the name of the pdb file to be include in output CSV file
      pdb_short <- strsplit(entry$pdbfile, split = "/") %>% 
        unlist() %>% tail(n=1) %>% strsplit("\\.") %>% unlist() %>% head(n=1)
      contacts <- binding.site(pdb,
                             a.inds=atom.select(pdb,chain=entry$chains_RBD),
                             b.inds=atom.select(pdb,chain=entry$chains_ligand),
                             cutoff=entry$distance_cutoff, hydrogens=F)
      if (is.vector(contacts)){
        structural_contacts <- rbind(structural_contacts, 
                                     data.frame(name=entry$name, 
                                                pdb=pdb_short, 
                                                chain = contacts$chain, 
                                                position = contacts$resno)
                                     )
      }

    }

    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

Annotate sites on the antibody that contact the RBD
---------------------------------------------------

This is slightly more complicated because each antibody has multiple
chains. We examine each chain individually, and output a separate CSV
file that includes the chain designation. In addition, there are some
antibody chains that do not contact the RBD at all, so we must deal with
this edge case.

IMPORTANT DISCLAIMER!!!!
------------------------

Importantly, some PDB files (such as is the case for the CR3022
antibody, PDB 6W41) annotate the CDR loops with letters that are
stripped when read by the `bio3d` package so you should check each PDB
file individually to make sure this is not the case, otherwise there
will be problems with this script.

    # make a separate file for the sites on the antibody that contact RBD
    for(entry in contacts_config){
      pdb <- read.pdb(file=entry$pdbfile)
      # get the name of the pdb file to be include in output CSV file
      pdb_short <- strsplit(entry$pdbfile, split = "/") %>% 
        unlist() %>% tail(n=1) %>% strsplit("\\.") %>% unlist() %>% head(n=1)
      contacts <- binding.site(pdb,
                             a.inds=atom.select(pdb,chain=entry$chains_ligand),
                             b.inds=atom.select(pdb,chain=entry$chains_RBD),
                             cutoff=entry$distance_cutoff, hydrogens=F)
      if (is.vector(contacts)){
        antibody_contacts <- rbind(antibody_contacts, 
                                   data.frame(name=entry$name, 
                                              pdb=pdb_short, 
                                              chain = contacts$chain, 
                                              position = contacts$resno))

      }

     }

    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

Save list of structural contacts
--------------------------------

    write.csv(structural_contacts,file=config$structural_contacts,row.names=F,quote=F)
    write.csv(antibody_contacts,file=config$antibody_contacts,row.names=F,quote=F)
