# Antibody binding (or escape) to homologs.
This Python Jupyter notebook sees how well each antibody / sera binds to the homologs in the libraries.

## Set up analysis
Import Python modules:


```python
import math
import os

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read escape profiles config, which tells which sets to make plots for:


```python
with open(config['escape_profiles_config']) as f:
    escape_profiles_config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['homologs_dir'], exist_ok=True)
```

Read homolog escape fractions:


```python
print(f"Reading homolog escape fractions from {config['escape_fracs_homologs']}")
escape_fracs = (
    pd.read_csv(config['escape_fracs_homologs'])
    .query('library == "average"')
    .drop(columns=['library', 'nlibs'])
    .rename(columns={'selection': 'condition'})
    .assign(homolog=lambda x: x['homolog'].map(lambda h: 'SARS-CoV-1' if h == 'SARS-CoV' else h))
    .assign(homolog=lambda x: pd.Categorical(x['homolog'], x['homolog'].unique(), ordered=True))
    )
```

    Reading homolog escape fractions from results/escape_scores/escape_fracs_homologs.csv


## Escape from homologs for all antibodies / sera
Plot how well each homolg escapes each antibody or sera:


```python
p = (ggplot(escape_fracs) +
     aes('homolog', 'condition', fill='escape_frac') +
     geom_tile(color='black', size=0.3) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(0.3 * escape_fracs['homolog'].nunique(), 0.3 * escape_fracs['condition'].nunique()),
           ) +
     xlab('') +
     ylab('') +
     scale_fill_gradient(low='white', high='#353D41', limits=(0, 1))
     )
_ = p.draw()
```


    
![png](homolog_escape_files/homolog_escape_12_0.png)
    


## Escape on specified antibody / sera subsets
We perform analyses on all subsets in the escape profiles config:


```python
for name, specs in escape_profiles_config.items():

    plotfile = os.path.join(config['homologs_dir'], f"{name}_homolog_escape.pdf")
    print(f"\nAnalyzing homolog escape for {name}, saving plot to {plotfile}")
    
    conditions = specs['conditions']
    
    # get count for conditions of interest for this subset
    df = (escape_fracs
          .query('condition in @conditions')
          .assign(condition=lambda x: x['condition'].map(conditions))
          .assign(condition=lambda x: pd.Categorical(x['condition'], list(conditions.values()), ordered=True))
          )
    
    if len(df) == 0:
        raise RuntimeError(f"no homolog data for {name}")
        
    p = (ggplot(df) +
         aes('homolog', 'condition', fill='escape_frac') +
         geom_tile(color='black', size=0.3) +
         theme(axis_text_x=element_text(angle=90),
               figure_size=(0.3 * df['homolog'].nunique(), 0.3 * df['condition'].nunique()),
               panel_background=element_blank(),
           ) +
         xlab('') +
         ylab('') +
         scale_fill_gradient(low='white', high='#353D41', limits=(0, 1))
         )
    p.save(plotfile, verbose=False)
    fig = p.draw()
    display(fig)
    plt.close(fig)
```

    
    Analyzing homolog escape for human_sera_oldnames, saving plot to results/homologs/human_sera_oldnames_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_1.png)
    


    
    Analyzing homolog escape for human_sera, saving plot to results/homologs/human_sera_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_3.png)
    


    
    Analyzing homolog escape for human_sera_figS3, saving plot to results/homologs/human_sera_figS3_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_5.png)
    


    
    Analyzing homolog escape for human_sera_fig4A_1C, saving plot to results/homologs/human_sera_fig4A_1C_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_7.png)
    


    
    Analyzing homolog escape for human_sera_fig4B_24C, saving plot to results/homologs/human_sera_fig4B_24C_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_9.png)
    


    
    Analyzing homolog escape for human_sera_fig4C_25, saving plot to results/homologs/human_sera_fig4C_25_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_11.png)
    


    
    Analyzing homolog escape for human_sera_fig4D_7C, saving plot to results/homologs/human_sera_fig4D_7C_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_13.png)
    


    
    Analyzing homolog escape for human_sera_fig4E_25C, saving plot to results/homologs/human_sera_fig4E_25C_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_15.png)
    


    
    Analyzing homolog escape for human_sera_fig4X_6C, saving plot to results/homologs/human_sera_fig4X_6C_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_17.png)
    


    
    Analyzing homolog escape for human_sera_early_fig2, saving plot to results/homologs/human_sera_early_fig2_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_19.png)
    


    
    Analyzing homolog escape for 2165_validation, saving plot to results/homologs/2165_validation_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_21.png)
    


    
    Analyzing homolog escape for CB6_validation, saving plot to results/homologs/CB6_validation_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_23.png)
    


    
    Analyzing homolog escape for 2082_validation, saving plot to results/homologs/2082_validation_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_25.png)
    


    
    Analyzing homolog escape for 2832_validation, saving plot to results/homologs/2832_validation_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_27.png)
    


    
    Analyzing homolog escape for human_sera_pink, saving plot to results/homologs/human_sera_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_29.png)
    


    
    Analyzing homolog escape for human_sera_figS3_pink, saving plot to results/homologs/human_sera_figS3_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_31.png)
    


    
    Analyzing homolog escape for human_sera_fig4A_1C_pink, saving plot to results/homologs/human_sera_fig4A_1C_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_33.png)
    


    
    Analyzing homolog escape for human_sera_fig4B_24C_pink, saving plot to results/homologs/human_sera_fig4B_24C_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_35.png)
    


    
    Analyzing homolog escape for human_sera_fig4C_25_pink, saving plot to results/homologs/human_sera_fig4C_25_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_37.png)
    


    
    Analyzing homolog escape for human_sera_fig4D_7C_pink, saving plot to results/homologs/human_sera_fig4D_7C_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_39.png)
    


    
    Analyzing homolog escape for human_sera_fig4E_25C_pink, saving plot to results/homologs/human_sera_fig4E_25C_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_41.png)
    


    
    Analyzing homolog escape for human_sera_fig4X_6C_pink, saving plot to results/homologs/human_sera_fig4X_6C_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_43.png)
    


    
    Analyzing homolog escape for human_sera_early_fig2_pink, saving plot to results/homologs/human_sera_early_fig2_pink_homolog_escape.pdf



    
![png](homolog_escape_files/homolog_escape_14_45.png)
    



```python

```
