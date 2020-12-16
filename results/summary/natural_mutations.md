# Analyze naturally occurring mutations at sites of strong escape
This Python Jupyter notebook sees how many naturally occuring mutations are observed at each site of strong escape

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
os.makedirs(config['gisaid_mutations_dir'], exist_ok=True)
```

Read counts of naturally ocurring mutations:


```python
print(f"Reading mutation counts from {config['gisaid_mutation_counts']}")

mut_counts = pd.read_csv(config['gisaid_mutation_counts'])
```

    Reading mutation counts from results/GISAID_mutations/mutation_counts.csv


Read sites of "strong escape" from all antibodies / sera:


```python
print(f"Reading sites of strong escape from {config['strong_escape_sites']}")

strong_sites = pd.read_csv(config['strong_escape_sites'])
```

    Reading sites of strong escape from results/escape_profiles/strong_escape_sites.csv


Read escape fractions for all antibodies / sera:


```python
print(f"Reading escape fractions from {config['escape_fracs']}")

escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .drop(columns='site')
    .rename(columns={'mutation': 'mutant',
                     'label_site': 'site'})
    [['condition', 'site', 'wildtype', 'mutant', config['mut_metric'], config['site_metric']]]
    )

escape_fracs
```

    Reading escape fractions from results/escape_scores/escape_fracs.csv





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>condition</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>mut_escape_frac_epistasis_model</th>
      <th>site_total_escape_frac_epistasis_model</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.002020</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>1</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.005616</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>2</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.002535</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>3</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.003032</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.003113</td>
      <td>0.04926</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>71990</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>R</td>
      <td>0.002448</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>71991</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>S</td>
      <td>0.002555</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>71992</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>V</td>
      <td>0.002227</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>71993</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>W</td>
      <td>0.001916</td>
      <td>0.04056</td>
    </tr>
    <tr>
      <th>71994</th>
      <td>REGN10987_400</td>
      <td>531</td>
      <td>T</td>
      <td>Y</td>
      <td>0.002564</td>
      <td>0.04056</td>
    </tr>
  </tbody>
</table>
<p>71995 rows × 6 columns</p>
</div>



## Counts of mutations at sites of escape
Get counts of naturally occurring mutations at sites of escape, along with the actual escape values:

First get mutation-level counts:


```python
mutcounts_strong_sites = (
    strong_sites[['condition', 'threshold', 'site']]
    .merge(mut_counts, how='inner', on='site')
    .merge(escape_fracs[['condition', 'site', 'wildtype', config['site_metric']]].drop_duplicates(),
           on=['condition', 'site', 'wildtype'],
           validate='many_to_one')
    .assign(mutation=lambda x: x['wildtype'] + x['site'].astype(str) + x['mutant'])
    .sort_values('count', ascending=False)
    )
```

Now get site-level counts (aggregating all mutations at a site):


```python
sitecounts_strong_sites = (
    mutcounts_strong_sites
    .assign(mut_count=lambda x: x['mutation'] + ' (' + x['count'].astype(str) + ')')
    .groupby(['condition', 'threshold', 'site', 'wildtype', config['site_metric']])
    .aggregate({'count': 'sum', 'mut_count': ', '.join})
    .rename(columns={'mut_count': 'counts_by_mutation'})
    .reset_index()
    .sort_values('count', ascending=False)
    )

print(f"Here are first few lines showing the most frequently mutated sites of escape:")
display(HTML(sitecounts_strong_sites.head(n=20).to_html(index=False)))
```

    Here are first few lines showing the most frequently mutated sites of escape:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>site_total_escape_frac_epistasis_model</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>25_d18_500</td>
      <td>sensitive_max_mut</td>
      <td>439</td>
      <td>N</td>
      <td>0.6198</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>sensitive_max_mut</td>
      <td>439</td>
      <td>N</td>
      <td>2.1660</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>sensitive</td>
      <td>439</td>
      <td>N</td>
      <td>0.6198</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>sensitive</td>
      <td>439</td>
      <td>N</td>
      <td>2.1660</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>REGN10987_400</td>
      <td>default</td>
      <td>439</td>
      <td>N</td>
      <td>2.1660</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>23C_d102_80</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>0.5154</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>23C_d102_80</td>
      <td>sensitive_max_mut</td>
      <td>453</td>
      <td>Y</td>
      <td>0.5154</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>23C_d26_80</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>0.4554</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>default</td>
      <td>453</td>
      <td>Y</td>
      <td>4.1080</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>23C_d26_80</td>
      <td>sensitive_max_mut</td>
      <td>453</td>
      <td>Y</td>
      <td>0.4554</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>4.1080</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>REGN10933_400</td>
      <td>sensitive_max_mut</td>
      <td>453</td>
      <td>Y</td>
      <td>4.1080</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>0.9038</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>25_d18_500</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.6856</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>CB6_400</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>1.5790</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>25C_d115_80</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.9038</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>COV2-2499_400</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>0.4765</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>22C_d28_200</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>1.2370</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>COV2-2499_400</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>0.4765</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>22C_d28_200</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>1.2370</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
  </tbody>
</table>


Now plot mutation counts (any mutation) at each site of escape for each antibody / sera:


```python
nconditions = sitecounts_strong_sites['condition'].nunique() * sitecounts_strong_sites['threshold'].nunique()
ncol = 8
nrow = math.ceil(nconditions / ncol)

p = (ggplot(sitecounts_strong_sites) +
     aes(config['site_metric'], 'count') +
     geom_point(alpha=0.5) +
     facet_wrap('~ condition + threshold', ncol=ncol) +
     theme(figure_size=(2 * ncol, 2 * nrow)) +
     xlab('site-level escape') +
     ylab('sequences with mutation at site')
     )

_ = p.draw()
```


    
![png](natural_mutations_files/natural_mutations_20_0.png)
    


## Perform analyses on subsets
We perform analyses on all subsets in the escape profiles config for which this is specified:


```python
for name, specs in escape_profiles_config.items():
    if 'analyze_natural_mutations' not in specs or not specs['analyze_natural_mutations']:
        continue
    print(f"\nAnalyzing natural mutations for {name}")
    
    conditions = specs['conditions']
    
    threshold = specs['plot_auto_identified_sites']
    if threshold not in sitecounts_strong_sites['threshold'].unique():
        raise ValueError(f"invalid threshold {threshold} for {name}")
    
    # get count for conditions of interest for this subset
    df = (sitecounts_strong_sites
          .query('condition in @conditions')
          .query('threshold == @threshold')
          .assign(condition=lambda x: x['condition'].map(conditions))
          .drop(columns=config['site_metric'])
          )
    countsfile = os.path.join(config['gisaid_mutations_dir'], f"{name}_mutation_counts.csv")
    print(f"Writing counts of mutations at sites of strong escape to {countsfile}. First few lines:")
    display(HTML(df.head(n=10).to_html(index=False)))
    df.to_csv(countsfile, index=False)
    
    # make plot showing escape sites with more than mincounts mutations
    if 'natural_mutations_mincounts' in specs:
        mincounts = specs['natural_mutations_mincounts']
    else:
        mincounts = 5
    plotsfile = os.path.join(config['gisaid_mutations_dir'], f"{name}_mutation_counts.pdf")
    print('Plotting which antibodies / sera are escaped by mutations at all sites of '
          f"escape with at least {mincounts} mutation counts and saving to {plotsfile}.")
    plot_df = (
        # data frame with all combinations of conditions and sites
        pd.DataFrame.from_records([(condition, site) for condition in conditions.values()
                                   for site in df['site'].unique()],
                                  columns=['condition', 'site'])
        # annotate sites of escape
        .merge(df.assign(escape=lambda x: x['count'] >= mincounts)[['condition', 'site', 'escape']],
               how='left',
               validate='one_to_one',
               on=['condition', 'site'])
        .assign(escape=lambda x: x['escape'].fillna(False))
        # add wildtype and counts of mutations at each site
        .merge(sitecounts_strong_sites[['site', 'wildtype', 'count']].drop_duplicates(),
               how='left',
               validate='many_to_one',
               on='site')
        # get only sites with sufficient mutation counts
        .query('count > @mincounts')
        # only get sites where at least one antibody escapes
        .assign(n_escape=lambda x: x.groupby('site')['escape'].transform('sum'))
        .query('n_escape > 0')
        # order conditions, and order sites by count after making nice label
        .assign(site_label=lambda x: x['wildtype'] + x['site'].astype(str) + ' (' + x['count'].astype(str) + ')')
        .sort_values('count')
        .assign(condition=lambda x: pd.Categorical(x['condition'], list(conditions.values()), ordered=True),
                site_label=lambda x: pd.Categorical(x['site_label'], x['site_label'].unique(), ordered=True)
                )
        )
    p = (ggplot(plot_df) +
         aes('condition', 'site_label', fill='escape') +
         geom_tile(color='black', size=0.3) +
         theme(axis_text_x=element_text(angle=90),
               figure_size=(0.3 * plot_df['condition'].nunique(), 0.3 * plot_df['site_label'].nunique()),
               panel_background=element_blank(),
               ) +
         xlab('') +
         ylab('') +
         scale_fill_manual(values=['white', 'dimgray'])
         )
    p.save(plotsfile, verbose=False)
    fig = p.draw()
    display(fig)
    plt.close(fig)
```

    
    Analyzing natural mutations for human_sera
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/human_sera_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>12C_d152</td>
      <td>default</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>12C_d152</td>
      <td>default</td>
      <td>494</td>
      <td>S</td>
      <td>91</td>
      <td>S494P (85), S494L (5), S494A (1)</td>
    </tr>
    <tr>
      <td>23C_d102</td>
      <td>default</td>
      <td>494</td>
      <td>S</td>
      <td>91</td>
      <td>S494P (85), S494L (5), S494A (1)</td>
    </tr>
    <tr>
      <td>13_d15</td>
      <td>default</td>
      <td>384</td>
      <td>P</td>
      <td>86</td>
      <td>P384L (66), P384S (20)</td>
    </tr>
    <tr>
      <td>23_d21</td>
      <td>default</td>
      <td>384</td>
      <td>P</td>
      <td>86</td>
      <td>P384L (66), P384S (20)</td>
    </tr>
    <tr>
      <td>1C_d26</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>66</td>
      <td>E484Q (29), E484K (27), E484A (6), E484D (3), E484R (1)</td>
    </tr>
    <tr>
      <td>23_d21</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>66</td>
      <td>E484Q (29), E484K (27), E484A (6), E484D (3), E484R (1)</td>
    </tr>
    <tr>
      <td>23C_d102</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>66</td>
      <td>E484Q (29), E484K (27), E484A (6), E484D (3), E484R (1)</td>
    </tr>
    <tr>
      <td>22C_d28</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>66</td>
      <td>E484Q (29), E484K (27), E484A (6), E484D (3), E484R (1)</td>
    </tr>
    <tr>
      <td>22C_d104</td>
      <td>default</td>
      <td>484</td>
      <td>E</td>
      <td>66</td>
      <td>E484Q (29), E484K (27), E484A (6), E484D (3), E484R (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/human_sera_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_22_3.png)
    


    
    Analyzing natural mutations for human_sera_sensitive
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/human_sera_sensitive_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>25_d18</td>
      <td>sensitive</td>
      <td>439</td>
      <td>N</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>23C_d102</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>23C_d26</td>
      <td>sensitive</td>
      <td>453</td>
      <td>Y</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>25_d18</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>25C_d115</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>22C_d28</td>
      <td>sensitive</td>
      <td>501</td>
      <td>N</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>23C_d102</td>
      <td>sensitive</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>12C_d152</td>
      <td>sensitive</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>23C_d26</td>
      <td>sensitive</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>22C_d104</td>
      <td>sensitive</td>
      <td>494</td>
      <td>S</td>
      <td>91</td>
      <td>S494P (85), S494L (5), S494A (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/human_sera_sensitive_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_22_7.png)
    


    
    Analyzing natural mutations for human_sera_sensitive_max_mut
    Writing counts of mutations at sites of strong escape to results/GISAID_mutations/human_sera_sensitive_max_mut_mutation_counts.csv. First few lines:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>condition</th>
      <th>threshold</th>
      <th>site</th>
      <th>wildtype</th>
      <th>count</th>
      <th>counts_by_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>25_d18</td>
      <td>sensitive_max_mut</td>
      <td>439</td>
      <td>N</td>
      <td>2242</td>
      <td>N439K (2242)</td>
    </tr>
    <tr>
      <td>23C_d102</td>
      <td>sensitive_max_mut</td>
      <td>453</td>
      <td>Y</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>23C_d26</td>
      <td>sensitive_max_mut</td>
      <td>453</td>
      <td>Y</td>
      <td>330</td>
      <td>Y453F (330)</td>
    </tr>
    <tr>
      <td>25C_d115</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>22C_d28</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>25_d18</td>
      <td>sensitive_max_mut</td>
      <td>501</td>
      <td>N</td>
      <td>225</td>
      <td>N501Y (210), N501T (10), N501S (5)</td>
    </tr>
    <tr>
      <td>23C_d26</td>
      <td>sensitive_max_mut</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>23C_d102</td>
      <td>sensitive_max_mut</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>12C_d152</td>
      <td>sensitive_max_mut</td>
      <td>483</td>
      <td>V</td>
      <td>106</td>
      <td>V483A (50), V483F (49), V483L (4), V483I (3)</td>
    </tr>
    <tr>
      <td>25C_d115</td>
      <td>sensitive_max_mut</td>
      <td>494</td>
      <td>S</td>
      <td>91</td>
      <td>S494P (85), S494L (5), S494A (1)</td>
    </tr>
  </tbody>
</table>


    Plotting which antibodies / sera are escaped by mutations at all sites of escape with at least 5 mutation counts and saving to results/GISAID_mutations/human_sera_sensitive_max_mut_mutation_counts.pdf.



    
![png](natural_mutations_files/natural_mutations_22_11.png)
    



```python

```
