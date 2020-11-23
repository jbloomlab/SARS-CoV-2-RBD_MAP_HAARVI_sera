# Analyze counts and compute escape scores
This Python Jupyter notebook analyzes the variant counts and looks at mutation coverage and jackpotting.
It then computes an "escape scores" for each variant after grouping by barcode or substitutions as specified in the configuration.

## Set up analysis

This notebook primarily makes use of the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package, and uses [plotnine](https://github.com/has2k1/plotnine) for ggplot2-like plotting syntax:


```python
import collections
import math
import os
import warnings

import Bio.SeqIO

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import numpy

import pandas as pd

from plotnine import *

import seaborn

import yaml
```

Set [plotnine](https://github.com/has2k1/plotnine) theme to the gray-grid one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using dms_variants version 0.8.5


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['escape_scores_dir'], exist_ok=True)
```

Read information about the samples:


```python
samples_df = pd.read_csv(config['barcode_runs'])
```

## Initialize codon-variant table
Initialize [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) from wildtype gene sequence and the variant counts CSV file.
We will then use the plotting functions of this variant table to analyze the counts per sample:


```python
wt_seqrecord = Bio.SeqIO.read(config['wildtype_sequence'], 'fasta')
geneseq = str(wt_seqrecord.seq)
primary_target = wt_seqrecord.name
print(f"Read sequence of {len(geneseq)} nt for {primary_target} from {config['wildtype_sequence']}")
      
print(f"Initializing CodonVariantTable from gene sequence and {config['variant_counts']}")
      
variants = dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df(
                geneseq=geneseq,
                variant_count_df_file=config['variant_counts'],
                primary_target=primary_target)
      
print('Done initializing CodonVariantTable.')
```

    Read sequence of 603 nt for SARS-CoV-2 from data/wildtype_sequence.fasta
    Initializing CodonVariantTable from gene sequence and results/counts/variant_counts.csv
    Done initializing CodonVariantTable.


## Sequencing counts per sample
Average counts per variant for each sample.
Note that these are the **sequencing** counts, in some cases they may outstrip the actual number of sorted cells:


```python
p = variants.plotAvgCountsPerVariant(libraries=variants.libraries,
                                     by_target=False,
                                     orientation='v')
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_19_0.png)
    


And the numerical values plotted above:


```python
display(HTML(
 variants.avgCountsPerVariant(libraries=variants.libraries,
                               by_target=False)
 .pivot_table(index='sample',
              columns='library',
              values='avg_counts_per_variant')
 .round(1)
 .to_html()
 ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>lib1</th>
      <th>lib2</th>
    </tr>
    <tr>
      <th>sample</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>expt_34-41-none-0-reference</th>
      <td>548.6</td>
      <td>552.0</td>
    </tr>
    <tr>
      <th>expt_46-52-none-0-reference</th>
      <td>357.5</td>
      <td>367.3</td>
    </tr>
    <tr>
      <th>expt_24-33-none-0-reference</th>
      <td>541.5</td>
      <td>593.9</td>
    </tr>
    <tr>
      <th>expt_42-45-none-0-reference</th>
      <td>316.3</td>
      <td>353.8</td>
    </tr>
    <tr>
      <th>expt_13-16-none-0-reference</th>
      <td>250.4</td>
      <td>245.7</td>
    </tr>
    <tr>
      <th>expt_7-11-none-0-reference</th>
      <td>340.0</td>
      <td>358.2</td>
    </tr>
    <tr>
      <th>expt_22-23-none-0-reference</th>
      <td>1063.9</td>
      <td>1041.3</td>
    </tr>
    <tr>
      <th>expt_3-6-none-0-reference</th>
      <td>324.3</td>
      <td>364.0</td>
    </tr>
    <tr>
      <th>expt_40-12C_d152-80-escape</th>
      <td>34.7</td>
      <td>30.9</td>
    </tr>
    <tr>
      <th>expt_37-12C_d61-160-escape</th>
      <td>18.0</td>
      <td>21.1</td>
    </tr>
    <tr>
      <th>expt_52-13_d121-1250-escape</th>
      <td>20.0</td>
      <td>22.7</td>
    </tr>
    <tr>
      <th>expt_31-13_d15-200-escape</th>
      <td>27.3</td>
      <td>35.3</td>
    </tr>
    <tr>
      <th>expt_51-1C_d113-200-escape</th>
      <td>14.4</td>
      <td>21.3</td>
    </tr>
    <tr>
      <th>expt_41-1C_d26-200-escape</th>
      <td>20.0</td>
      <td>31.6</td>
    </tr>
    <tr>
      <th>expt_45-22C_d104-200-escape</th>
      <td>24.9</td>
      <td>30.1</td>
    </tr>
    <tr>
      <th>expt_38-22C_d28-200-escape</th>
      <td>29.8</td>
      <td>26.4</td>
    </tr>
    <tr>
      <th>expt_47-23C_d102-80-escape</th>
      <td>14.5</td>
      <td>32.0</td>
    </tr>
    <tr>
      <th>expt_46-23C_d26-80-escape</th>
      <td>21.4</td>
      <td>26.2</td>
    </tr>
    <tr>
      <th>expt_50-23_d120-500-escape</th>
      <td>20.4</td>
      <td>28.8</td>
    </tr>
    <tr>
      <th>expt_34-23_d21-1250-escape</th>
      <td>13.9</td>
      <td>25.2</td>
    </tr>
    <tr>
      <th>expt_39-23_d45-1250-escape</th>
      <td>25.1</td>
      <td>23.1</td>
    </tr>
    <tr>
      <th>expt_44-24C_d104-200-escape</th>
      <td>24.6</td>
      <td>28.4</td>
    </tr>
    <tr>
      <th>expt_35-24C_d32-200-escape</th>
      <td>30.5</td>
      <td>31.0</td>
    </tr>
    <tr>
      <th>expt_49-25C_d115-80-escape</th>
      <td>14.2</td>
      <td>28.3</td>
    </tr>
    <tr>
      <th>expt_48-25C_d48-200-escape</th>
      <td>18.9</td>
      <td>23.7</td>
    </tr>
    <tr>
      <th>expt_32-25_d18-500-escape</th>
      <td>38.0</td>
      <td>36.2</td>
    </tr>
    <tr>
      <th>expt_33-25_d94-200-escape</th>
      <td>21.5</td>
      <td>24.9</td>
    </tr>
    <tr>
      <th>expt_30-6C_d33-500-escape</th>
      <td>31.4</td>
      <td>29.7</td>
    </tr>
    <tr>
      <th>expt_42-6C_d76-500-escape</th>
      <td>28.7</td>
      <td>23.3</td>
    </tr>
    <tr>
      <th>expt_43-7C_d103-200-escape</th>
      <td>21.8</td>
      <td>27.2</td>
    </tr>
    <tr>
      <th>expt_36-7C_d29-500-escape</th>
      <td>41.8</td>
      <td>30.5</td>
    </tr>
    <tr>
      <th>expt_27-CB6-400-escape</th>
      <td>67.1</td>
      <td>37.0</td>
    </tr>
    <tr>
      <th>expt_13-COV2-2050-400-escape</th>
      <td>47.9</td>
      <td>58.2</td>
    </tr>
    <tr>
      <th>expt_11-COV2-2082-400-escape</th>
      <td>74.8</td>
      <td>43.5</td>
    </tr>
    <tr>
      <th>expt_14-COV2-2094-400-escape</th>
      <td>34.0</td>
      <td>29.6</td>
    </tr>
    <tr>
      <th>expt_23-COV2-2096-400-escape</th>
      <td>134.5</td>
      <td>133.4</td>
    </tr>
    <tr>
      <th>expt_9-COV2-2165-400-escape</th>
      <td>105.8</td>
      <td>34.2</td>
    </tr>
    <tr>
      <th>expt_5-COV2-2479-400-escape</th>
      <td>176.1</td>
      <td>182.3</td>
    </tr>
    <tr>
      <th>expt_15-COV2-2499-400-escape</th>
      <td>67.9</td>
      <td>88.5</td>
    </tr>
    <tr>
      <th>expt_16-COV2-2677-400-escape</th>
      <td>40.5</td>
      <td>40.8</td>
    </tr>
    <tr>
      <th>expt_10-COV2-2832-400-escape</th>
      <td>154.7</td>
      <td>89.7</td>
    </tr>
    <tr>
      <th>expt_3-CR3022-400-escape</th>
      <td>170.3</td>
      <td>78.2</td>
    </tr>
    <tr>
      <th>expt_24-REGN10933-400-escape</th>
      <td>44.7</td>
      <td>35.2</td>
    </tr>
    <tr>
      <th>expt_26-REGN10933+REGN10987-400-escape</th>
      <td>30.7</td>
      <td>24.8</td>
    </tr>
    <tr>
      <th>expt_25-REGN10987-400-escape</th>
      <td>40.1</td>
      <td>53.5</td>
    </tr>
  </tbody>
</table>


## Mutations per variant
Average number of mutations per gene among all variants of the primary target, separately for each date:


```python
# this plotting is very slow when lots of samples, so for now plots are commented out

#for date, date_df in samples_df.groupby('date', sort=False):
#    p = variants.plotNumCodonMutsByType(variant_type='all',
#                                        orientation='v',
#                                        libraries=variants.libraries,
#                                        samples=date_df['sample'].unique().tolist(),
#                                        widthscale=2)
#    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
#    fig = p.draw()
#    display(fig)
#    plt.close(fig)
```

Now similar plots but showing mutation frequency across the gene:


```python
# this plotting is very slow when lots of samples, so for now code commented out

# for date, date_df in samples_df.groupby('date', sort=False):
#    p = variants.plotMutFreqs(variant_type='all',
#                              mut_type='codon',
#                              orientation='v',
#                              libraries=variants.libraries,
#                              samples=date_df['sample'].unique().tolist(),
#                              widthscale=1.5)
#    fig = p.draw()
#    display(fig)
#    plt.close(fig)
```

## Jackpotting and mutation coverage in pre-selection libraries
We look at the distribution of counts in the "reference" (pre-selection) libraries to see if they seem jackpotted (a few variants at very high frequency):


```python
pre_samples_df = samples_df.query('selection == "reference"')
```

Distribution of mutations along the gene for the pre-selection samples; big spikes may indicate jackpotting:


```python
# this plotting is very slow when lots of samples, so for now code commented out

#p = variants.plotMutFreqs(variant_type='all',
#                          mut_type='codon',
#                          orientation='v',
#                          libraries=variants.libraries,
#                          samples=pre_samples_df['sample'].unique().tolist(),
#                          widthscale=1.5)
#_ = p.draw()
```

How many mutations are observed frequently in pre-selection libraries?
Note that the libraries have been pre-selected for ACE2 binding, so we expect stop variants to mostly be missing.
Make the plot both for all variants and just single-mutant variants:


```python
# this plotting is very slow when lots of samples, so for now code commented out

#for variant_type in ['all', 'single']:
#    p = variants.plotCumulMutCoverage(
#                          variant_type=variant_type,
#                          mut_type='aa',
#                          orientation='v',
#                          libraries=variants.libraries,
#                          samples=pre_samples_df['sample'].unique().tolist(),
#                          widthscale=1.8,
#                          heightscale=1.2)
#    _ = p.draw()
```

Now make a plot showing what number and fraction of counts are for each variant in each pre-selection sample / library.
If some variants constitute a very high fraction, then that indicates extensive bottlenecking:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_33_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_33_1.png)
    


Now make the same plot breaking down by variant class, which enables determination of which types of variants are at high and low frequencies.
For this plot (unlike one above not classified by category) we only show variants of the primary target (not the homologs), and also group synonymous with wildtype in order to reduce number of plotted categories to make more interpretable:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      by_variant_class=True,
                                      classifyVariants_kwargs={'syn_as_wt': True},
                                      primary_target_only=True,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_35_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_35_1.png)
    


We also directly look to see what is the variant in each reference library / sample with the highest fraction counts.
Knowing if the highest frequency variant is shared helps determine **where** in the experiment the jackpotting happened:


```python
frac_counts_per_variant = (
        variants.add_frac_counts(variants.variant_count_df)
        .query(f"sample in {pre_samples_df['sample'].unique().tolist()}")
        )

display(HTML(
    frac_counts_per_variant
    .sort_values('frac_counts', ascending=False)
    .groupby(['library', 'sample'])
    .head(n=1)
    .sort_values(['library', 'sample'])
    .set_index(['library', 'sample'])
    [['frac_counts', 'target', 'barcode', 'aa_substitutions', 'codon_substitutions']]
    .round(4)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>frac_counts</th>
      <th>target</th>
      <th>barcode</th>
      <th>aa_substitutions</th>
      <th>codon_substitutions</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="8" valign="top">lib1</th>
      <th>expt_34-41-none-0-reference</th>
      <td>0.0037</td>
      <td>SARS-CoV-2</td>
      <td>TAGTAACAATGCGGTA</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>expt_46-52-none-0-reference</th>
      <td>0.0051</td>
      <td>SARS-CoV-2</td>
      <td>TTCCAAAATATTGTCA</td>
      <td>D59N F156S</td>
      <td>GAT59AAC TTT156TCG</td>
    </tr>
    <tr>
      <th>expt_24-33-none-0-reference</th>
      <td>0.0049</td>
      <td>SARS-CoV-2</td>
      <td>TTCCAAAATATTGTCA</td>
      <td>D59N F156S</td>
      <td>GAT59AAC TTT156TCG</td>
    </tr>
    <tr>
      <th>expt_42-45-none-0-reference</th>
      <td>0.0046</td>
      <td>SARS-CoV-2</td>
      <td>TTCCAAAATATTGTCA</td>
      <td>D59N F156S</td>
      <td>GAT59AAC TTT156TCG</td>
    </tr>
    <tr>
      <th>expt_13-16-none-0-reference</th>
      <td>0.0646</td>
      <td>SARS-CoV-2</td>
      <td>TTCCAAAATATTGTCA</td>
      <td>D59N F156S</td>
      <td>GAT59AAC TTT156TCG</td>
    </tr>
    <tr>
      <th>expt_7-11-none-0-reference</th>
      <td>0.1332</td>
      <td>SARS-CoV-2</td>
      <td>TAGTAACAATGCGGTA</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>expt_22-23-none-0-reference</th>
      <td>0.0669</td>
      <td>SARS-CoV-2</td>
      <td>TAGTAACAATGCGGTA</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>expt_3-6-none-0-reference</th>
      <td>0.1310</td>
      <td>SARS-CoV-2</td>
      <td>TAGTAACAATGCGGTA</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th rowspan="8" valign="top">lib2</th>
      <th>expt_34-41-none-0-reference</th>
      <td>0.0012</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_46-52-none-0-reference</th>
      <td>0.0020</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_24-33-none-0-reference</th>
      <td>0.0011</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_42-45-none-0-reference</th>
      <td>0.0018</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_13-16-none-0-reference</th>
      <td>0.0754</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_7-11-none-0-reference</th>
      <td>0.1213</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_22-23-none-0-reference</th>
      <td>0.0619</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
    <tr>
      <th>expt_3-6-none-0-reference</th>
      <td>0.1219</td>
      <td>SARS-CoV-2</td>
      <td>CAGTACAAAAGTATAA</td>
      <td>K87N V173Y</td>
      <td>CTA60CTG AAG87AAC GTC173TAC</td>
    </tr>
  </tbody>
</table>


To further where the jackpotting relative to the generation of the reference samples, we plot the correlation among the fraction of counts for the different reference samples.
If the fractions are highly correlated, that indicates that the jackpotting occurred in some upstream step common to the reference samples:


```python
# this code makes a full matrix of scatter plots, but is REALLY SLOW. So for now,
# it is commented out in favor of code that just makes correlation matrix.
#for lib, lib_df in frac_counts_per_variant.groupby('library'):
#    wide_lib_df = lib_df.pivot_table(index=['target', 'barcode'],
#                                     columns='sample',
#                                     values='frac_counts')
#    g = seaborn.pairplot(wide_lib_df, corner=True, plot_kws={'alpha': 0.5}, diag_kind='kde')
#    _ = g.fig.suptitle(lib, size=18)
#    plt.show()
```

## Examine counts for wildtype variants
The type of score we use to quantify escape depends on how well represented wildtype is in the selected libraries.
If wildtype is still well represented, we can use a more conventional functional score that gives differential selection relative to wildtype.
If wildtype is not well represented, then we need an alternative score that does not involve normalizing frequencies to wildtype.

First get average fraction of counts per variant for each variant class:


```python
counts_by_class = (
    variants.variant_count_df
    .pipe(variants.add_frac_counts)
    .pipe(variants.classifyVariants,
          primary_target=variants.primary_target,
          non_primary_target_class='homolog',
          class_as_categorical=True)
    .groupby(['library', 'sample', 'variant_class'])
    .aggregate(avg_frac_counts=pd.NamedAgg('frac_counts', 'mean'))
    .reset_index()
    .merge(samples_df[['sample', 'library', 'date', 'antibody', 'concentration', 'selection']],
           on=['sample', 'library'], validate='many_to_one')
    )

display(HTML(counts_by_class.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>variant_class</th>
      <th>avg_frac_counts</th>
      <th>date</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>selection</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>wildtype</td>
      <td>0.000026</td>
      <td>200910</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>synonymous</td>
      <td>0.000024</td>
      <td>200910</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>1 nonsynonymous</td>
      <td>0.000020</td>
      <td>200910</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>&gt;1 nonsynonymous</td>
      <td>0.000008</td>
      <td>200910</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>stop</td>
      <td>0.000002</td>
      <td>200910</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
  </tbody>
</table>


Plot average fraction of all counts per variant for each variant class.
If the values for wildtype are low for the non-reference samples (such as more similar to stop the nonsynonymous), then normalizing by wildtype in calculating scores will probably not work well as wildtype is too depleted:


```python
min_frac = 1e-7  # plot values < this as this

p = (ggplot(counts_by_class
            .assign(avg_frac_counts=lambda x: numpy.clip(x['avg_frac_counts'], min_frac, None))
            ) +
     aes('avg_frac_counts', 'sample', color='selection') +
     geom_point(size=2) +
     scale_color_manual(values=CBPALETTE[1:]) +
     facet_grid('library ~ variant_class') +
     scale_x_log10() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.5 * counts_by_class['variant_class'].nunique(),
                        0.2 * counts_by_class['library'].nunique() * 
                        counts_by_class['sample'].nunique())
           ) +
     geom_vline(xintercept=min_frac, linetype='dotted', color=CBPALETTE[3])
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_43_0.png)
    


## Compute escape scores
We use the escape score metric, which does **not** involve normalizing to wildtype and so isn't strongly affected by low wildtype counts.
We compute the scores using the method [dms_variants.codonvarianttable.CodonVariantTable.escape_scores](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores).

First, define what samples to compare for each calculation, matching each post-selection (escape) to the pre-selection (reference) sample on the same date:


```python
score_sample_df = (
    samples_df
    .query('selection == "escape"')
    .rename(columns={'sample': 'post_sample',
                     'cells_sorted': 'pre_cells_sorted'})
    .merge(samples_df
           .query('selection == "reference"')
           [['sample', 'library', 'date', 'cells_sorted']]
           .rename(columns={'sample': 'pre_sample',
                            'cells_sorted': 'post_cells_sorted'}),
           how='left', on=['date', 'library'], validate='many_to_one',
           )
    .assign(name=lambda x: x['antibody'] + '_' + x['concentration'].astype(str))
    # add dates to names where needed to make unique
    .assign(n_libs=lambda x: x.groupby(['name', 'date'])['pre_sample'].transform('count'))
    .sort_values(['name', 'date', 'n_libs'], ascending=False)
    .assign(i_name=lambda x: x.groupby(['library', 'name'], sort=False)['pre_sample'].cumcount(),
            name=lambda x: x.apply(lambda r: r['name'] + '_' + str(r['date']) if r['i_name'] > 0 else r['name'],
                                   axis=1),
            )
    .sort_values(['antibody', 'concentration', 'library', 'i_name'])
    # get columns of interest
    [['name', 'library', 'antibody', 'concentration', 'concentration_units', 'date',
      'pre_sample', 'post_sample', 'frac_escape', 'pre_cells_sorted', 'post_cells_sorted']]
    )

assert len(score_sample_df.groupby(['name', 'library'])) == len(score_sample_df)

print(f"Writing samples used to compute escape scores to {config['escape_score_samples']}\n")
score_sample_df.to_csv(config['escape_score_samples'], index=False)

display(HTML(score_sample_df.to_html(index=False)))
```

    Writing samples used to compute escape scores to results/escape_scores/samples.csv
    



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>library</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>concentration_units</th>
      <th>date</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>frac_escape</th>
      <th>pre_cells_sorted</th>
      <th>post_cells_sorted</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>12C_d152_80</td>
      <td>lib1</td>
      <td>12C_d152</td>
      <td>80</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>0.040</td>
      <td>530574.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>lib2</td>
      <td>12C_d152</td>
      <td>80</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>0.052</td>
      <td>528440.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>12C_d61_160</td>
      <td>lib1</td>
      <td>12C_d61</td>
      <td>160</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_37-12C_d61-160-escape</td>
      <td>0.033</td>
      <td>289524.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>12C_d61_160</td>
      <td>lib2</td>
      <td>12C_d61</td>
      <td>160</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_37-12C_d61-160-escape</td>
      <td>0.031</td>
      <td>373381.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>13_d121_1250</td>
      <td>lib1</td>
      <td>13_d121</td>
      <td>1250</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_52-13_d121-1250-escape</td>
      <td>0.042</td>
      <td>521978.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>13_d121_1250</td>
      <td>lib2</td>
      <td>13_d121</td>
      <td>1250</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_52-13_d121-1250-escape</td>
      <td>0.048</td>
      <td>631476.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>13_d15_200</td>
      <td>lib1</td>
      <td>13_d15</td>
      <td>200</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_31-13_d15-200-escape</td>
      <td>0.044</td>
      <td>501386.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>13_d15_200</td>
      <td>lib2</td>
      <td>13_d15</td>
      <td>200</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_31-13_d15-200-escape</td>
      <td>0.057</td>
      <td>502894.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>1C_d113_200</td>
      <td>lib1</td>
      <td>1C_d113</td>
      <td>200</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_51-1C_d113-200-escape</td>
      <td>0.037</td>
      <td>445318.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>1C_d113_200</td>
      <td>lib2</td>
      <td>1C_d113</td>
      <td>200</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_51-1C_d113-200-escape</td>
      <td>0.047</td>
      <td>576459.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>1C_d26_200</td>
      <td>lib1</td>
      <td>1C_d26</td>
      <td>200</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_41-1C_d26-200-escape</td>
      <td>0.045</td>
      <td>599453.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>1C_d26_200</td>
      <td>lib2</td>
      <td>1C_d26</td>
      <td>200</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_41-1C_d26-200-escape</td>
      <td>0.036</td>
      <td>522412.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>22C_d104_200</td>
      <td>lib1</td>
      <td>22C_d104</td>
      <td>200</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_45-22C_d104-200-escape</td>
      <td>0.050</td>
      <td>645592.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>22C_d104_200</td>
      <td>lib2</td>
      <td>22C_d104</td>
      <td>200</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_45-22C_d104-200-escape</td>
      <td>0.044</td>
      <td>637035.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>22C_d28_200</td>
      <td>lib1</td>
      <td>22C_d28</td>
      <td>200</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_38-22C_d28-200-escape</td>
      <td>0.047</td>
      <td>532618.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>22C_d28_200</td>
      <td>lib2</td>
      <td>22C_d28</td>
      <td>200</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_38-22C_d28-200-escape</td>
      <td>0.029</td>
      <td>502292.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>23C_d102_80</td>
      <td>lib1</td>
      <td>23C_d102</td>
      <td>80</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_47-23C_d102-80-escape</td>
      <td>0.030</td>
      <td>370503.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>23C_d102_80</td>
      <td>lib2</td>
      <td>23C_d102</td>
      <td>80</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_47-23C_d102-80-escape</td>
      <td>0.058</td>
      <td>727991.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>23C_d26_80</td>
      <td>lib1</td>
      <td>23C_d26</td>
      <td>80</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_46-23C_d26-80-escape</td>
      <td>0.027</td>
      <td>450907.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>23C_d26_80</td>
      <td>lib2</td>
      <td>23C_d26</td>
      <td>80</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_46-23C_d26-80-escape</td>
      <td>0.046</td>
      <td>649814.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>23_d120_500</td>
      <td>lib1</td>
      <td>23_d120</td>
      <td>500</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_50-23_d120-500-escape</td>
      <td>0.045</td>
      <td>467275.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>23_d120_500</td>
      <td>lib2</td>
      <td>23_d120</td>
      <td>500</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_50-23_d120-500-escape</td>
      <td>0.064</td>
      <td>774180.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>23_d21_1250</td>
      <td>lib1</td>
      <td>23_d21</td>
      <td>1250</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_34-23_d21-1250-escape</td>
      <td>0.026</td>
      <td>203778.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>23_d21_1250</td>
      <td>lib2</td>
      <td>23_d21</td>
      <td>1250</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_34-23_d21-1250-escape</td>
      <td>0.017</td>
      <td>228514.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>23_d45_1250</td>
      <td>lib1</td>
      <td>23_d45</td>
      <td>1250</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_39-23_d45-1250-escape</td>
      <td>0.031</td>
      <td>459593.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>23_d45_1250</td>
      <td>lib2</td>
      <td>23_d45</td>
      <td>1250</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_39-23_d45-1250-escape</td>
      <td>0.025</td>
      <td>331160.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>24C_d104_200</td>
      <td>lib1</td>
      <td>24C_d104</td>
      <td>200</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_44-24C_d104-200-escape</td>
      <td>0.046</td>
      <td>641483.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>24C_d104_200</td>
      <td>lib2</td>
      <td>24C_d104</td>
      <td>200</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_44-24C_d104-200-escape</td>
      <td>0.055</td>
      <td>722686.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>24C_d32_200</td>
      <td>lib1</td>
      <td>24C_d32</td>
      <td>200</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_35-24C_d32-200-escape</td>
      <td>0.067</td>
      <td>569827.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>24C_d32_200</td>
      <td>lib2</td>
      <td>24C_d32</td>
      <td>200</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_35-24C_d32-200-escape</td>
      <td>0.065</td>
      <td>507246.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>lib1</td>
      <td>25C_d115</td>
      <td>80</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_49-25C_d115-80-escape</td>
      <td>0.044</td>
      <td>421799.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>lib2</td>
      <td>25C_d115</td>
      <td>80</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_49-25C_d115-80-escape</td>
      <td>0.059</td>
      <td>677714.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>25C_d48_200</td>
      <td>lib1</td>
      <td>25C_d48</td>
      <td>200</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_48-25C_d48-200-escape</td>
      <td>0.040</td>
      <td>466102.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>25C_d48_200</td>
      <td>lib2</td>
      <td>25C_d48</td>
      <td>200</td>
      <td>dilution</td>
      <td>200917</td>
      <td>expt_46-52-none-0-reference</td>
      <td>expt_48-25C_d48-200-escape</td>
      <td>0.050</td>
      <td>571563.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>lib1</td>
      <td>25_d18</td>
      <td>500</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_32-25_d18-500-escape</td>
      <td>0.051</td>
      <td>547164.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>lib2</td>
      <td>25_d18</td>
      <td>500</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_32-25_d18-500-escape</td>
      <td>0.059</td>
      <td>562520.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>25_d94_200</td>
      <td>lib1</td>
      <td>25_d94</td>
      <td>200</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_33-25_d94-200-escape</td>
      <td>0.052</td>
      <td>382787.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>25_d94_200</td>
      <td>lib2</td>
      <td>25_d94</td>
      <td>200</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_33-25_d94-200-escape</td>
      <td>0.042</td>
      <td>422273.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>6C_d33_500</td>
      <td>lib1</td>
      <td>6C_d33</td>
      <td>500</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_30-6C_d33-500-escape</td>
      <td>0.051</td>
      <td>589736.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>6C_d33_500</td>
      <td>lib2</td>
      <td>6C_d33</td>
      <td>500</td>
      <td>dilution</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_30-6C_d33-500-escape</td>
      <td>0.047</td>
      <td>506366.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>6C_d76_500</td>
      <td>lib1</td>
      <td>6C_d76</td>
      <td>500</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_42-6C_d76-500-escape</td>
      <td>0.043</td>
      <td>600854.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>6C_d76_500</td>
      <td>lib2</td>
      <td>6C_d76</td>
      <td>500</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_42-6C_d76-500-escape</td>
      <td>0.038</td>
      <td>601796.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>7C_d103_200</td>
      <td>lib1</td>
      <td>7C_d103</td>
      <td>200</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_43-7C_d103-200-escape</td>
      <td>0.049</td>
      <td>653383.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>7C_d103_200</td>
      <td>lib2</td>
      <td>7C_d103</td>
      <td>200</td>
      <td>dilution</td>
      <td>200915</td>
      <td>expt_42-45-none-0-reference</td>
      <td>expt_43-7C_d103-200-escape</td>
      <td>0.049</td>
      <td>654893.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>7C_d29_500</td>
      <td>lib1</td>
      <td>7C_d29</td>
      <td>500</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_36-7C_d29-500-escape</td>
      <td>0.061</td>
      <td>625811.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>7C_d29_500</td>
      <td>lib2</td>
      <td>7C_d29</td>
      <td>500</td>
      <td>dilution</td>
      <td>200910</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_36-7C_d29-500-escape</td>
      <td>0.061</td>
      <td>535256.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>lib1</td>
      <td>CB6</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_27-CB6-400-escape</td>
      <td>0.222</td>
      <td>1907893.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>lib2</td>
      <td>CB6</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_27-CB6-400-escape</td>
      <td>0.225</td>
      <td>927804.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>lib1</td>
      <td>COV2-2050</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_13-COV2-2050-400-escape</td>
      <td>0.078</td>
      <td>1253185.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>lib2</td>
      <td>COV2-2050</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_13-COV2-2050-400-escape</td>
      <td>0.104</td>
      <td>1588138.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2082_400</td>
      <td>lib1</td>
      <td>COV2-2082</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200616</td>
      <td>expt_7-11-none-0-reference</td>
      <td>expt_11-COV2-2082-400-escape</td>
      <td>0.057</td>
      <td>522474.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2082_400</td>
      <td>lib2</td>
      <td>COV2-2082</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200616</td>
      <td>expt_7-11-none-0-reference</td>
      <td>expt_11-COV2-2082-400-escape</td>
      <td>0.086</td>
      <td>512138.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2094_400</td>
      <td>lib1</td>
      <td>COV2-2094</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_14-COV2-2094-400-escape</td>
      <td>0.069</td>
      <td>853420.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2094_400</td>
      <td>lib2</td>
      <td>COV2-2094</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_14-COV2-2094-400-escape</td>
      <td>0.071</td>
      <td>954885.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2096_400</td>
      <td>lib1</td>
      <td>COV2-2096</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200810</td>
      <td>expt_22-23-none-0-reference</td>
      <td>expt_23-COV2-2096-400-escape</td>
      <td>0.155</td>
      <td>1018697.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2096_400</td>
      <td>lib2</td>
      <td>COV2-2096</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200810</td>
      <td>expt_22-23-none-0-reference</td>
      <td>expt_23-COV2-2096-400-escape</td>
      <td>0.110</td>
      <td>927621.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2165_400</td>
      <td>lib1</td>
      <td>COV2-2165</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200616</td>
      <td>expt_7-11-none-0-reference</td>
      <td>expt_9-COV2-2165-400-escape</td>
      <td>0.071</td>
      <td>911426.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2165_400</td>
      <td>lib2</td>
      <td>COV2-2165</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200616</td>
      <td>expt_7-11-none-0-reference</td>
      <td>expt_9-COV2-2165-400-escape</td>
      <td>0.029</td>
      <td>393068.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2479_400</td>
      <td>lib1</td>
      <td>COV2-2479</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200615</td>
      <td>expt_3-6-none-0-reference</td>
      <td>expt_5-COV2-2479-400-escape</td>
      <td>0.050</td>
      <td>684856.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2479_400</td>
      <td>lib2</td>
      <td>COV2-2479</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200615</td>
      <td>expt_3-6-none-0-reference</td>
      <td>expt_5-COV2-2479-400-escape</td>
      <td>0.044</td>
      <td>582089.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2499_400</td>
      <td>lib1</td>
      <td>COV2-2499</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_15-COV2-2499-400-escape</td>
      <td>0.177</td>
      <td>2593041.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2499_400</td>
      <td>lib2</td>
      <td>COV2-2499</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_15-COV2-2499-400-escape</td>
      <td>0.185</td>
      <td>2655138.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2677_400</td>
      <td>lib1</td>
      <td>COV2-2677</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_16-COV2-2677-400-escape</td>
      <td>0.086</td>
      <td>1243786.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2677_400</td>
      <td>lib2</td>
      <td>COV2-2677</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200720</td>
      <td>expt_13-16-none-0-reference</td>
      <td>expt_16-COV2-2677-400-escape</td>
      <td>0.088</td>
      <td>1457488.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2832_400</td>
      <td>lib1</td>
      <td>COV2-2832</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200616</td>
      <td>expt_7-11-none-0-reference</td>
      <td>expt_10-COV2-2832-400-escape</td>
      <td>0.116</td>
      <td>1525920.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>COV2-2832_400</td>
      <td>lib2</td>
      <td>COV2-2832</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200616</td>
      <td>expt_7-11-none-0-reference</td>
      <td>expt_10-COV2-2832-400-escape</td>
      <td>0.076</td>
      <td>975400.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>CR3022_400</td>
      <td>lib1</td>
      <td>CR3022</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200615</td>
      <td>expt_3-6-none-0-reference</td>
      <td>expt_3-CR3022-400-escape</td>
      <td>0.148</td>
      <td>1790556.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>CR3022_400</td>
      <td>lib2</td>
      <td>CR3022</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200615</td>
      <td>expt_3-6-none-0-reference</td>
      <td>expt_3-CR3022-400-escape</td>
      <td>0.091</td>
      <td>1175925.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>lib1</td>
      <td>REGN10933</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_24-REGN10933-400-escape</td>
      <td>0.174</td>
      <td>1252342.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>lib2</td>
      <td>REGN10933</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_24-REGN10933-400-escape</td>
      <td>0.154</td>
      <td>984237.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>REGN10933+REGN10987_400</td>
      <td>lib1</td>
      <td>REGN10933+REGN10987</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_26-REGN10933+REGN10987-400-escape</td>
      <td>0.101</td>
      <td>635898.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>REGN10933+REGN10987_400</td>
      <td>lib2</td>
      <td>REGN10933+REGN10987</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_26-REGN10933+REGN10987-400-escape</td>
      <td>0.083</td>
      <td>585188.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>lib1</td>
      <td>REGN10987</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_25-REGN10987-400-escape</td>
      <td>0.148</td>
      <td>1155004.0</td>
      <td>160000000.0</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>lib2</td>
      <td>REGN10987</td>
      <td>400</td>
      <td>ng_per_mL</td>
      <td>200904</td>
      <td>expt_24-33-none-0-reference</td>
      <td>expt_25-REGN10987-400-escape</td>
      <td>0.150</td>
      <td>1091502.0</td>
      <td>160000000.0</td>
    </tr>
  </tbody>
</table>


Compute the escape scores for variants of the primary target and classify the variants, then compute escape scores for homologs:


```python
print(f"Computing escape scores for {primary_target} variants using {config['escape_score_type']} "
      f"score type with a pseudocount of {config['escape_score_pseudocount']} and "
      f"an escape fraction floor {config['escape_score_floor_E']}, an escape fraction ceiling "
      f"{config['escape_score_ceil_E']}, and grouping variants by {config['escape_score_group_by']}.")

escape_scores = (variants.escape_scores(score_sample_df,
                                        score_type=config['escape_score_type'],
                                        pseudocount=config['escape_score_pseudocount'],
                                        floor_E=config['escape_score_floor_E'],
                                        ceil_E=config['escape_score_ceil_E'],
                                        by=config['escape_score_group_by'],
                                        )
                 .query('target == @primary_target')
                 .pipe(variants.classifyVariants,
                       primary_target=variants.primary_target,
                       syn_as_wt=(config['escape_score_group_by'] == 'aa_substitutions'),
                       )
                 )
print('Here are the first few lines of the resulting escape scores:')
display(HTML(escape_scores.head().to_html(index=False)))

print(f"\nComputing scores for homologs grouping by {config['escape_score_homolog_group_by']}")

escape_scores_homologs = (
        variants.escape_scores(score_sample_df,
                               score_type=config['escape_score_type'],
                               pseudocount=config['escape_score_pseudocount'],
                               floor_E=config['escape_score_floor_E'],
                               ceil_E=config['escape_score_ceil_E'],
                               by=config['escape_score_homolog_group_by'],
                               )
        .query('(target != @primary_target) | (n_aa_substitutions == 0)')
        )

print('Here are the first few lines of the resulting homolog escape scores:')
display(HTML(escape_scores_homologs.head().to_html(index=False)))
```

    Computing escape scores for SARS-CoV-2 variants using frac_escape score type with a pseudocount of 0.5 and an escape fraction floor 0, an escape fraction ceiling 1, and grouping variants by barcode.
    Here are the first few lines of the resulting escape scores:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>12C_d152_80</td>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>TAGTAACAATGCGGTA</td>
      <td>0.002678</td>
      <td>8.346082e-09</td>
      <td>200832</td>
      <td>862</td>
      <td></td>
      <td>0</td>
      <td></td>
      <td>0</td>
      <td>wildtype</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>TTCCAAAATATTGTCA</td>
      <td>0.104824</td>
      <td>4.563043e-07</td>
      <td>164928</td>
      <td>27725</td>
      <td>GAT59AAC TTT156TCG</td>
      <td>2</td>
      <td>D59N F156S</td>
      <td>2</td>
      <td>&gt;1 nonsynonymous</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>TTAATTAGTATCAGGT</td>
      <td>0.023192</td>
      <td>2.063353e-07</td>
      <td>72580</td>
      <td>2699</td>
      <td>GAC34CGC CTG125GTC</td>
      <td>2</td>
      <td>D34R L125V</td>
      <td>2</td>
      <td>&gt;1 nonsynonymous</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>TAATATGCAAGAGCCA</td>
      <td>0.044899</td>
      <td>6.533481e-07</td>
      <td>45851</td>
      <td>3301</td>
      <td>CCT96TGC</td>
      <td>1</td>
      <td>P96C</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>ATAAAAAGTCCATATG</td>
      <td>0.003782</td>
      <td>5.419200e-08</td>
      <td>43774</td>
      <td>265</td>
      <td>CTT5AAG</td>
      <td>1</td>
      <td>L5K</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
  </tbody>
</table>


    
    Computing scores for homologs grouping by aa_substitutions
    Here are the first few lines of the resulting homolog escape scores:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>aa_substitutions</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>12C_d152_80</td>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>BM48-31</td>
      <td>0.413966</td>
      <td>6.587496e-05</td>
      <td>6521</td>
      <td>4315</td>
      <td>0</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>GD-Pangolin</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>GD-Pangolin</td>
      <td>0.003058</td>
      <td>1.457457e-08</td>
      <td>131846</td>
      <td>644</td>
      <td>0</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>HKU3-1</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>HKU3-1</td>
      <td>0.479859</td>
      <td>7.505510e-05</td>
      <td>7053</td>
      <td>5410</td>
      <td>0</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>LYRa11</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>LYRa11</td>
      <td>0.522851</td>
      <td>1.786961e-05</td>
      <td>33292</td>
      <td>27825</td>
      <td>0</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>RaTG13</td>
      <td>lib1</td>
      <td>expt_34-41-none-0-reference</td>
      <td>expt_40-12C_d152-80-escape</td>
      <td>RaTG13</td>
      <td>0.330916</td>
      <td>2.678892e-06</td>
      <td>115308</td>
      <td>60995</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


## Apply pre-selection count filter to variant escape scores
Now determine a pre-selection count filter in order to flag for removal variants with counts that are so low that the estimated score is probably noise.
We know that stop codons should be largely purged pre-selection, and so the counts for them are a good indication of the "noise" threshold.
We therefore set the filter using the number of pre-selection counts for the stop codons.

To do this, we first compute the number of pre-selection counts for stop-codon variants at various quantiles and look at these.
We then take the number of pre-selection counts at the specified quantile as the filter cutoff, and filter scores for all variants with pre-selection counts less than this filter cutoff:


```python
filter_quantile = config['escape_score_stop_quantile_filter']
assert 0 <= filter_quantile <= 1

quantiles = sorted(set([0.5, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999] + [filter_quantile]))

stop_score_counts = (
    escape_scores
    .query('variant_class == "stop"')
    .groupby(['library', 'pre_sample'], observed=True)
    ['pre_count']
    .quantile(q=quantiles)
    .reset_index()
    .rename(columns={'level_2': 'quantile'})
    .pivot_table(index=['pre_sample', 'library'],
                 columns='quantile',
                 values='pre_count')
    )

print('Quantiles of the number of pre-selection counts per variant for stop variants:')
display(HTML(stop_score_counts.to_html(float_format='%.1f')))

print(f"\nSetting the pre-count filter cutoff to the {filter_quantile} quantile:")
pre_count_filter_cutoffs = (
    stop_score_counts
    [filter_quantile]
    .rename('pre_count_filter_cutoff')
    .reset_index()
    )
display(HTML(pre_count_filter_cutoffs.to_html(float_format='%.1f')))
```

    Quantiles of the number of pre-selection counts per variant for stop variants:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>quantile</th>
      <th>0.5</th>
      <th>0.9</th>
      <th>0.95</th>
      <th>0.98</th>
      <th>0.99</th>
      <th>0.995</th>
      <th>0.999</th>
    </tr>
    <tr>
      <th>pre_sample</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">expt_34-41-none-0-reference</th>
      <th>lib1</th>
      <td>61.0</td>
      <td>195.0</td>
      <td>254.0</td>
      <td>329.0</td>
      <td>388.0</td>
      <td>488.0</td>
      <td>1246.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>59.0</td>
      <td>203.0</td>
      <td>265.0</td>
      <td>349.0</td>
      <td>413.0</td>
      <td>488.0</td>
      <td>930.0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_46-52-none-0-reference</th>
      <th>lib1</th>
      <td>37.0</td>
      <td>118.0</td>
      <td>153.0</td>
      <td>203.0</td>
      <td>244.0</td>
      <td>304.0</td>
      <td>771.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>36.0</td>
      <td>125.0</td>
      <td>164.0</td>
      <td>218.0</td>
      <td>262.0</td>
      <td>330.0</td>
      <td>901.0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_24-33-none-0-reference</th>
      <th>lib1</th>
      <td>60.0</td>
      <td>195.0</td>
      <td>250.0</td>
      <td>322.0</td>
      <td>387.0</td>
      <td>484.0</td>
      <td>1096.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>63.0</td>
      <td>215.0</td>
      <td>283.0</td>
      <td>366.0</td>
      <td>446.0</td>
      <td>532.0</td>
      <td>1089.0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_42-45-none-0-reference</th>
      <th>lib1</th>
      <td>34.0</td>
      <td>111.0</td>
      <td>145.0</td>
      <td>191.0</td>
      <td>230.0</td>
      <td>284.9</td>
      <td>681.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>37.0</td>
      <td>128.0</td>
      <td>165.0</td>
      <td>220.0</td>
      <td>269.0</td>
      <td>320.7</td>
      <td>733.6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_13-16-none-0-reference</th>
      <th>lib1</th>
      <td>0.0</td>
      <td>10.0</td>
      <td>16.0</td>
      <td>27.0</td>
      <td>42.0</td>
      <td>74.0</td>
      <td>432.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0</td>
      <td>9.0</td>
      <td>15.0</td>
      <td>25.0</td>
      <td>37.0</td>
      <td>63.3</td>
      <td>505.2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_7-11-none-0-reference</th>
      <th>lib1</th>
      <td>0.0</td>
      <td>3.0</td>
      <td>7.0</td>
      <td>14.0</td>
      <td>21.0</td>
      <td>43.3</td>
      <td>870.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0</td>
      <td>3.0</td>
      <td>6.0</td>
      <td>12.0</td>
      <td>19.0</td>
      <td>38.5</td>
      <td>1129.4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_22-23-none-0-reference</th>
      <th>lib1</th>
      <td>0.0</td>
      <td>48.0</td>
      <td>76.0</td>
      <td>127.9</td>
      <td>174.4</td>
      <td>290.4</td>
      <td>1740.8</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0</td>
      <td>44.0</td>
      <td>72.0</td>
      <td>123.0</td>
      <td>176.4</td>
      <td>307.7</td>
      <td>2108.0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">expt_3-6-none-0-reference</th>
      <th>lib1</th>
      <td>0.0</td>
      <td>4.0</td>
      <td>7.0</td>
      <td>13.0</td>
      <td>21.0</td>
      <td>48.4</td>
      <td>813.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0</td>
      <td>3.0</td>
      <td>5.0</td>
      <td>10.0</td>
      <td>16.7</td>
      <td>32.1</td>
      <td>888.3</td>
    </tr>
  </tbody>
</table>


    
    Setting the pre-count filter cutoff to the 0.99 quantile:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>pre_sample</th>
      <th>library</th>
      <th>pre_count_filter_cutoff</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>expt_34-41-none-0-reference</td>
      <td>lib1</td>
      <td>388.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>expt_34-41-none-0-reference</td>
      <td>lib2</td>
      <td>413.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>expt_46-52-none-0-reference</td>
      <td>lib1</td>
      <td>244.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>expt_46-52-none-0-reference</td>
      <td>lib2</td>
      <td>262.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>expt_24-33-none-0-reference</td>
      <td>lib1</td>
      <td>387.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>expt_24-33-none-0-reference</td>
      <td>lib2</td>
      <td>446.0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>expt_42-45-none-0-reference</td>
      <td>lib1</td>
      <td>230.0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>expt_42-45-none-0-reference</td>
      <td>lib2</td>
      <td>269.0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>expt_13-16-none-0-reference</td>
      <td>lib1</td>
      <td>42.0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>expt_13-16-none-0-reference</td>
      <td>lib2</td>
      <td>37.0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>expt_7-11-none-0-reference</td>
      <td>lib1</td>
      <td>21.0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>expt_7-11-none-0-reference</td>
      <td>lib2</td>
      <td>19.0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>expt_22-23-none-0-reference</td>
      <td>lib1</td>
      <td>174.4</td>
    </tr>
    <tr>
      <th>13</th>
      <td>expt_22-23-none-0-reference</td>
      <td>lib2</td>
      <td>176.4</td>
    </tr>
    <tr>
      <th>14</th>
      <td>expt_3-6-none-0-reference</td>
      <td>lib1</td>
      <td>21.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>expt_3-6-none-0-reference</td>
      <td>lib2</td>
      <td>16.7</td>
    </tr>
  </tbody>
</table>


Apply the filter to the escape scores, so that scores that fail the pre-selection count filter are now marked with `pass_pre_count_filter` of `False`:


```python
escape_scores = (
    escape_scores
    .merge(pre_count_filter_cutoffs,
           on=['library', 'pre_sample'],
           how='left',
           validate='many_to_one')
    .assign(pass_pre_count_filter=lambda x: x['pre_count'] >= x['pre_count_filter_cutoff'])
    )

escape_scores_homologs = (
    escape_scores_homologs
    .merge(pre_count_filter_cutoffs,
           on=['library', 'pre_sample'],
           how='left',
           validate='many_to_one')
    .assign(pass_pre_count_filter=lambda x: x['pre_count'] >= x['pre_count_filter_cutoff'])
    )
```

Plot the fraction of variants of each type that pass the pre-selection count filter in each pre-selection sample.
The ideal filter would have the property such that no *stop* variants pass, all *wildtype* (or *synonymous*) variants pass, and some intermediate fraction of *nonsynonymous* variants pass.
However, if the variant composition in the pre-selection samples is already heavily skewed by jackpotting, there will be some deviation from this ideal behavior.
Here is what the plots actually look like:


```python
frac_pre_pass_filter = (
    escape_scores
    [['pre_sample', 'library', 'target', config['escape_score_group_by'],
      'pre_count', 'pass_pre_count_filter', 'variant_class']]
    .drop_duplicates()
    .groupby(['pre_sample', 'library', 'variant_class'], observed=True)
    .aggregate(n_variants=pd.NamedAgg('pass_pre_count_filter', 'count'),
               n_pass_filter=pd.NamedAgg('pass_pre_count_filter', 'sum')
               )
    .reset_index()
    .assign(frac_pass_filter=lambda x: x['n_pass_filter'] / x['n_variants'],
            pre_sample=lambda x: pd.Categorical(x['pre_sample'], x['pre_sample'].unique(), ordered=True))
    )

p = (ggplot(frac_pre_pass_filter) +
     aes('variant_class', 'frac_pass_filter', fill='variant_class') +
     geom_bar(stat='identity') +
     facet_grid('library ~ pre_sample') +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(3.3 * frac_pre_pass_filter['pre_sample'].nunique(),
                        2 * frac_pre_pass_filter['library'].nunique()),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_53_0.png)
    


## Apply ACE2-binding / expression filter to variant mutations
In [Starr et al (2020)](https://www.biorxiv.org/content/10.1101/2020.06.17.157982v1), we used deep mutational scanning to estimate how each mutation affected ACE2 binding and expression.
Here we flag for removal any variants of the primary target that have (or have mutations) that were measured to decrease ACE2-binding or expression beyond a minimal threshold, in order to avoid these variants muddying the signal as spurious escape mutants.

To do this, we first determine all mutations that do / do-not having binding that exceeds the thresholds.
We then flag variants as passing the ACE2-binding / expression filter if all of the mutations they contain exceed both thresholds, and failing as otherwise.
In addition, we look at the measured ACE2-binding / expression for each variant in each library, and flag as passing the filter any variants that have binding / expression that exeed the thresholds, and failing as otherwise.
If we are grouping by amino-acid substitutions, we use the average value for variants with those substitutions.

Note that these filters are only applied to mutants of the primary target; homologs are specified manually for filtering in the configuration file:


```python
print(f"Reading ACE2-binding and expression for mutations from {config['mut_bind_expr']}, "
      f"and for variants from {config['variant_bind']} and {config['variant_expr']}, "
      f"and filtering for variants with binding >={config['escape_score_min_bind_variant']}."
      f"and expression >= {config['escape_score_min_expr_variant']}, and also variants that "
      f"only have mutations with binding >={config['escape_score_min_bind_mut']} and "
      f"expression >={config['escape_score_min_expr_mut']}.")

# filter on mutations
mut_bind_expr = pd.read_csv(config['mut_bind_expr'])
assert mut_bind_expr['mutation_RBD'].nunique() == len(mut_bind_expr)
for prop in ['bind', 'expr']:
    muts_adequate = set(mut_bind_expr
                        .query(f"{prop}_avg >= {config[f'escape_score_min_{prop}_mut']}")
                        ['mutation_RBD']
                        )
    print(f"{len(muts_adequate)} of {len(mut_bind_expr)} mutations have adequate {prop}.")
    escape_scores[f"muts_pass_{prop}_filter"] = (
        escape_scores
        ['aa_substitutions']
        .map(lambda s: set(s.split()).issubset(muts_adequate))
        ) 

# filter on variants
for prop, col in [('bind', 'delta_log10Ka'), ('expr', 'delta_ML_meanF')]:
    filter_name = f"variant_pass_{prop}_filter"
    variant_pass_df = (
        pd.read_csv(config[f"variant_{prop}"], keep_default_na=False, na_values=['NA'])
        .groupby(['library', 'target', config['escape_score_group_by']])
        .aggregate(val=pd.NamedAgg(col, 'mean'))
        .reset_index()
        .assign(pass_filter=lambda x: x['val'] >= config[f"escape_score_min_{prop}_variant"])
        .rename(columns={'pass_filter': filter_name,
                         'val': f"variant_{prop}"})
        )
    print(f"\nTotal variants of {primary_target} that pass {prop} filter:")
    display(HTML(
        variant_pass_df
        .groupby(['library', filter_name])
        .aggregate(n_variants=pd.NamedAgg(config['escape_score_group_by'], 'count'))
        .to_html()
        ))
    escape_scores = (
        escape_scores
        .drop(columns=filter_name, errors='ignore')
        .merge(variant_pass_df,
               how='left',
               validate='many_to_one',
               on=['library', 'target', config['escape_score_group_by']],
               )
        )
    assert escape_scores[filter_name].notnull().all()

# annotate as passing overall filter if passes all mutation and binding filters:
escape_scores['pass_ACE2bind_expr_filter'] = (
        escape_scores['muts_pass_bind_filter'] &
        escape_scores['muts_pass_expr_filter'] &
        escape_scores['variant_pass_bind_filter'] &
        escape_scores['variant_pass_expr_filter']
        )
```

    Reading ACE2-binding and expression for mutations from results/prior_DMS_data/mutant_ACE2binding_expression.csv, and for variants from results/prior_DMS_data/variant_ACE2binding.csv and results/prior_DMS_data/variant_expression.csv, and filtering for variants with binding >=-2.35.and expression >= -1.0, and also variants that only have mutations with binding >=-2.35 and expression >=-1.0.
    3422 of 4221 mutations have adequate bind.
    2328 of 4221 mutations have adequate expr.
    
    Total variants of SARS-CoV-2 that pass bind filter:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>n_variants</th>
    </tr>
    <tr>
      <th>library</th>
      <th>variant_pass_bind_filter</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">lib1</th>
      <th>False</th>
      <td>56963</td>
    </tr>
    <tr>
      <th>True</th>
      <td>41606</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">lib2</th>
      <th>False</th>
      <td>55964</td>
    </tr>
    <tr>
      <th>True</th>
      <td>40548</td>
    </tr>
  </tbody>
</table>


    
    Total variants of SARS-CoV-2 that pass expr filter:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>n_variants</th>
    </tr>
    <tr>
      <th>library</th>
      <th>variant_pass_expr_filter</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">lib1</th>
      <th>False</th>
      <td>75269</td>
    </tr>
    <tr>
      <th>True</th>
      <td>23300</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">lib2</th>
      <th>False</th>
      <td>74439</td>
    </tr>
    <tr>
      <th>True</th>
      <td>22073</td>
    </tr>
  </tbody>
</table>


Plot the fraction of variants that **have already passed the pre-count filter** that are filtered by the ACE2-binding or expression thresholds:


```python
frac_ACE2bind_expr_pass_filter = (
    escape_scores
    .query('pass_pre_count_filter == True')
    [['pre_sample', 'library', 'target', config['escape_score_group_by'],
      'pre_count', 'pass_ACE2bind_expr_filter', 'variant_class']]
    .drop_duplicates()
    .groupby(['pre_sample', 'library', 'variant_class'], observed=True)
    .aggregate(n_variants=pd.NamedAgg('pass_ACE2bind_expr_filter', 'count'),
               n_pass_filter=pd.NamedAgg('pass_ACE2bind_expr_filter', 'sum')
               )
    .reset_index()
    .assign(frac_pass_filter=lambda x: x['n_pass_filter'] / x['n_variants'],
            pre_sample=lambda x: pd.Categorical(x['pre_sample'], x['pre_sample'].unique(), ordered=True))
    )

p = (ggplot(frac_ACE2bind_expr_pass_filter) +
     aes('variant_class', 'frac_pass_filter', fill='variant_class') +
     geom_bar(stat='identity') +
     facet_grid('library ~ pre_sample') +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(3.3 * frac_ACE2bind_expr_pass_filter['pre_sample'].nunique(),
                        2 * frac_ACE2bind_expr_pass_filter['library'].nunique()),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_57_0.png)
    


## Examine and write escape scores
Plot the distribution of escape scores across variants of different classes **among those that pass both the pre-selection count filter and the ACE2-binding / expression filter**.
If things are working correctly, we don't expect escape in wildtype (or synonymous variants), but do expect escape for some small fraction of nonsynymous variants.
Also, we do not plot the scores for the stop codon variant class, as most stop-codon variants have already been filtered out so this category is largely noise:


```python
nfacets = len(escape_scores.groupby(['library', 'name']).nunique())
ncol = min(8, nfacets)
nrow = math.ceil(nfacets / ncol)

df = (escape_scores
      .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
      .query('variant_class != "stop"')
      )
     
p = (ggplot(df) +
     aes('variant_class', 'score', color='variant_class') +
     geom_boxplot(outlier_size=1.5, outlier_alpha=0.1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 3 * nrow),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:])
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_59_0.png)
    


Also, we want to see how much the high escape scores are correlated with simple coverage.
To do this, we plot the correlation between escape score and pre-selection count just for the nonsynonymous variants (which are the ones that we expect to have true escape).
The plots below have a lot of overplotting, but are still sufficient to test of the score is simply correlated with the pre-selection counts or not.
The hoped for result is that the escape score doesn't appear to be strongly correlated with pre-selection counts:


```python
p = (ggplot(escape_scores
            .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
            .query('variant_class in ["1 nonsynonymous", ">1 nonsynonymous"]')
            ) +
     aes('pre_count', 'score') +
     geom_point(alpha=0.1, size=1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 2.35 * nrow),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:]) +
     scale_x_log10()
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_61_0.png)
    


Write the escape scores to a file:


```python
print(f"Writing escape scores for {primary_target} to {config['escape_scores']}")
escape_scores.to_csv(config['escape_scores'], index=False, float_format='%.4g')

print(f"Writing escape scores for homologs to {config['escape_scores_homologs']}")
escape_scores_homologs.to_csv(config['escape_scores_homologs'], index=False, float_format='%.4g')
```

    Writing escape scores for SARS-CoV-2 to results/escape_scores/scores.csv
    Writing escape scores for homologs to results/escape_scores/scores_homologs.csv

