# Analyze viral escape-mutant selections
Analyze results from viral escape-mutant selections.

Import Python modules:


```python
import collections
import os

import Bio.SeqIO

import dms_variants.constants
from dms_variants.constants import CBPALETTE
from dms_variants.utils import single_nt_accessible

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import mizani

import numpy

import pandas as pd

from plotnine import *

import yaml
```

Read in configuration and then escape-mutant selection results:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
    
print(f"Reading escape-mutant selection results from {config['escape_selection_results']}")
with open(config['escape_selection_results']) as f:
    selection_results = yaml.safe_load(f)
```

    Reading escape-mutant selection results from data/escape_selection_results.yaml


Make output directory:


```python
os.makedirs(config['escape_selections_dir'], exist_ok=True)
```

Read escape-mutation mapping and deep mutational scanning results, and then merge them:


```python
# read escape fractions
escape_fracs = (
    pd.read_csv(config['escape_fracs'])
    .query('library == "average"')
    .rename(columns={config['mut_metric']: 'mutation_escape',
                     config['site_metric']: 'site_escape',
                     'selection': 'antibody'})
    .assign(site=lambda x: x['label_site'])
    [['antibody', 'site', 'wildtype', 'mutation', 'mutation_escape', 'site_escape']]
    )

# read DMS data
bind_expr = (
    pd.read_csv(config['mut_bind_expr'])
    .rename(columns={'site_SARS2': 'site',
                     'bind_avg': 'ACE2 binding',
                     'expr_avg': 'RBD expression',
                     })
    .assign(mutation=lambda x: x['mutant'])
    [['site', 'mutation', 'ACE2 binding', 'RBD expression']]
    )

# merge escape and DMS data
escape_dms = (
    escape_fracs
    .merge(bind_expr,
           on=['site', 'mutation'],
           how='left',
           validate='many_to_one',
           )
    )

# first few lines of data frame
display(HTML(escape_dms.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutation_escape</th>
      <th>site_escape</th>
      <th>ACE2 binding</th>
      <th>RBD expression</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.002020</td>
      <td>0.04926</td>
      <td>-0.03</td>
      <td>-0.11</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.005616</td>
      <td>0.04926</td>
      <td>0.03</td>
      <td>-0.44</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.002535</td>
      <td>0.04926</td>
      <td>0.00</td>
      <td>-0.31</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.003032</td>
      <td>0.04926</td>
      <td>-0.10</td>
      <td>-0.70</td>
    </tr>
    <tr>
      <td>12C_d152_80</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.003113</td>
      <td>0.04926</td>
      <td>-0.04</td>
      <td>-0.25</td>
    </tr>
  </tbody>
</table>


Get data frame with just escape-selection counts and then add to the main data frame of escape / DMS data:


```python
records = []
antibody_order = {}
spike_nt_seq = {}
for selection_set, selection_set_d in selection_results.items():
    antibody_order[selection_set] = {}
    spike_nt_seq[selection_set] = str(Bio.SeqIO.read(selection_set_d['spike_sequence'],
                                                     'genbank',
                                                     ).seq)
    for antibody, antibody_d in selection_set_d['antibodies'].items():
        antibody_order[selection_set][antibody] = antibody_d['display_name']
        if 'mutations' in antibody_d:
            for mutation_str, n in antibody_d['mutations'].items():
                wt = mutation_str[0]
                site = int(mutation_str[1: -1])
                mutation = mutation_str[-1]
                if ('label_mutations' in antibody_d) and mutation_str in antibody_d['label_mutations']:
                    mutation_label = mutation_str
                else:
                    mutation_label = None
                records.append((selection_set, antibody, site, wt, mutation, n, mutation_label))
                assert 1 == len(escape_dms.query('antibody == @antibody')
                                          .query('wildtype == @wt')
                                          .query('site == @site')
                                          .query('mutation == @mutation')
                                          ), f"{mutation_str} not in `escape_dms` once for {antibody}"
        if 'label_mutations' in antibody_d:
            for mutation_str in antibody_d['label_mutations']:
                if 'mutations' not in antibody_d or mutation_str not in antibody_d['mutations']:
                    wt = mutation_str[0]
                    site = int(mutation_str[1: -1])
                    mutation = mutation_str[-1]
                    records.append((selection_set, antibody, site, wt, mutation, 0, mutation_str))
                    assert 1 == len(escape_dms.query('antibody == @antibody')
                                              .query('wildtype == @wt')
                                              .query('site == @site')
                                              .query('mutation == @mutation')
                                              ), f"{mutation_str} not in `escape_dms` once for {antibody}"
            
selection_df = pd.DataFrame.from_records(
                records,
                columns=['selection_set', 'antibody', 'site',
                         'wildtype', 'mutation', 'n_selected', 'mutation_label'],
                )

display(HTML(selection_df.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection_set</th>
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>n_selected</th>
      <th>mutation_label</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Crowe_selections</td>
      <td>COV2-2094_400</td>
      <td>378</td>
      <td>K</td>
      <td>E</td>
      <td>7</td>
      <td>None</td>
    </tr>
    <tr>
      <td>Crowe_selections</td>
      <td>COV2-2094_400</td>
      <td>378</td>
      <td>K</td>
      <td>N</td>
      <td>1</td>
      <td>None</td>
    </tr>
    <tr>
      <td>Crowe_selections</td>
      <td>COV2-2479_400</td>
      <td>484</td>
      <td>E</td>
      <td>K</td>
      <td>3</td>
      <td>None</td>
    </tr>
    <tr>
      <td>Crowe_selections</td>
      <td>COV2-2050_400</td>
      <td>484</td>
      <td>E</td>
      <td>K</td>
      <td>4</td>
      <td>None</td>
    </tr>
    <tr>
      <td>Crowe_selections</td>
      <td>COV2-2499_400</td>
      <td>446</td>
      <td>G</td>
      <td>D</td>
      <td>3</td>
      <td>None</td>
    </tr>
    <tr>
      <td>Crowe_selections</td>
      <td>COV2-2499_400</td>
      <td>498</td>
      <td>Q</td>
      <td>R</td>
      <td>2</td>
      <td>None</td>
    </tr>
  </tbody>
</table>


Add escape-selection counts and labels to main data frame with escape scores, and just keep antibodies of interest for each selection:


```python
escape_dms_selection = (
    pd.concat([escape_dms.merge(df,
                                on=['antibody', 'site', 'wildtype', 'mutation'],
                                how='left',
                                validate='one_to_one',
                                )
                          .assign(selection_set=selection_set,
                                  antibody_name=lambda x: x['antibody'].map(antibody_order[selection_set])
                                  )
                          .query('antibody_name.notnull()')
               for selection_set, df in selection_df.groupby('selection_set')
               ])
    .assign(n_selected=lambda x: x['n_selected'].fillna(0).astype(int),
            n_selected_total=lambda x: x.groupby(['selection_set', 'antibody'])
                                        ['n_selected'].transform('sum'),
            any_selected=lambda x: x['n_selected'] > 0,
            )
    )

display(HTML(escape_dms_selection.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutation_escape</th>
      <th>site_escape</th>
      <th>ACE2 binding</th>
      <th>RBD expression</th>
      <th>selection_set</th>
      <th>n_selected</th>
      <th>mutation_label</th>
      <th>antibody_name</th>
      <th>n_selected_total</th>
      <th>any_selected</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.000325</td>
      <td>0.02808</td>
      <td>-0.03</td>
      <td>-0.11</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.001232</td>
      <td>0.02808</td>
      <td>0.03</td>
      <td>-0.44</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.005871</td>
      <td>0.02808</td>
      <td>0.00</td>
      <td>-0.31</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.002217</td>
      <td>0.02808</td>
      <td>-0.10</td>
      <td>-0.70</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.001924</td>
      <td>0.02808</td>
      <td>-0.04</td>
      <td>-0.25</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Now add in the viral codon sequence and determine what amino-acid mutations are single-nucleotide accessible:


```python
codon_df = pd.DataFrame()
sites = escape_dms_selection['site'].unique()
for selection_set, spike in spike_nt_seq.items():
    codon_df = codon_df.append(
        pd.DataFrame({'selection_set': selection_set, 'site': sites})
          .assign(viral_codon=lambda x: x['site'].map(lambda r: spike[3 * (r - 1): 3 * r]),
                  viral_aa=lambda x: x['viral_codon'].map(dms_variants.constants.CODON_TO_AA))
        )
    
escape_dms_selection = (
    escape_dms_selection
    .drop(columns=['viral_codon', 'viral_aa'], errors='ignore')
    .merge(codon_df,
           how='left',
           on=['selection_set', 'site'],
           validate='many_to_one',
           )
    .assign(single_nt_accessible=lambda x: x.apply(lambda row: single_nt_accessible(row['viral_codon'],
                                                                                    row['mutation']),
                                                   axis=1)
            )
    )

# check viral spike is same one used in our DMS. If not, we need to somehow deal with that...
aa_mismatch = escape_dms_selection.query('wildtype != viral_aa')
if len(aa_mismatch):
    raise ValueError('mismatches at the following amino acids: ' + aa_mismatch)

display(HTML(escape_dms_selection.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>antibody</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutation_escape</th>
      <th>site_escape</th>
      <th>ACE2 binding</th>
      <th>RBD expression</th>
      <th>selection_set</th>
      <th>n_selected</th>
      <th>mutation_label</th>
      <th>antibody_name</th>
      <th>n_selected_total</th>
      <th>any_selected</th>
      <th>viral_codon</th>
      <th>viral_aa</th>
      <th>single_nt_accessible</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>0.000325</td>
      <td>0.02808</td>
      <td>-0.03</td>
      <td>-0.11</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>0.001232</td>
      <td>0.02808</td>
      <td>0.03</td>
      <td>-0.44</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
      <td>AAT</td>
      <td>N</td>
      <td>True</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>0.005871</td>
      <td>0.02808</td>
      <td>0.00</td>
      <td>-0.31</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>0.002217</td>
      <td>0.02808</td>
      <td>-0.10</td>
      <td>-0.70</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
    <tr>
      <td>COV2-2050_400</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>0.001924</td>
      <td>0.02808</td>
      <td>-0.04</td>
      <td>-0.25</td>
      <td>Crowe_selections</td>
      <td>0</td>
      <td>NaN</td>
      <td>COV2-2050</td>
      <td>4</td>
      <td>False</td>
      <td>AAT</td>
      <td>N</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Make plots showing effects of all mutations stratified by how many times observed in selections, and whether single-nucleotide accessible or not from the viral spike.
For mutations observed in selections point size, is proportional to times observed:


```python
for selection_set, df in escape_dms_selection.groupby('selection_set'):
    
    min_size = selection_results[selection_set]['min_size']
    max_size = selection_results[selection_set]['max_size']
    size_scale = selection_results[selection_set]['max_size']
    
    print(f"\nMaking plot for {selection_set}")
    
    selected_not_accessible = len(df.query('any_selected and not single_nt_accessible'))
    
    if 'custom_categories' in selection_results[selection_set]:
        custom_cats = selection_results[selection_set]['custom_categories']
        addtl_cats = list({cat: True for a_cats in custom_cats.values() for cat in a_cats.values()})
    else:
        custom_cats = {}
        addtl_cats = []
    
    def get_point_category(row):
        if row['antibody_name'] in custom_cats and (row['mutation_str'] in 
                                                    custom_cats[row['antibody_name']]):
            return custom_cats[row['antibody_name']][row['mutation_str']]
        elif row['any_selected'] and row['single_nt_accessible']:
            return 'selected (single-nucleotide)'
        elif row['any_selected']:
            return 'selected (multi-nucleotide)'
        elif row['single_nt_accessible']:
            return 'single-nucleotide'
        else:
            return 'multi-nucleotide'
        
    df = df.assign(
            mutation_str=lambda x: x['wildtype'] + x['site'].astype(str) + x['mutation'],
            antibody_name=lambda x: pd.Categorical(x['antibody_name'],
                                                   antibody_order[selection_set].values(),
                                                   ordered=True),
            point_area=lambda x: 0.5 * size_scale * numpy.clip(x['n_selected'].fillna(0),
                                                               min_size,
                                                               max_size),
            point_category=lambda x: x.apply(get_point_category, axis=1),
            )
    possible_point_categories = ['single-nucleotide', 'multi-nucleotide',
                                 'selected (single-nucleotide)',
                                 'selected (multi-nucleotide)'] + addtl_cats
    observed_point_categories = [cat for cat in possible_point_categories
                                 if cat in df['point_category'].unique()]
    df = df.assign(point_category=lambda x: pd.Categorical(x['point_category'],
                                                           observed_point_categories,
                                                           ordered=True)
                   )
    
    if 'shapes' in selection_results[selection_set]:
        cat_shapes = selection_results[selection_set]['shapes']
    else:
        cat_shapes = ['o', 'x', 'D', '^']
    assert len(cat_shapes) >= len(observed_point_categories), 'not enough shapes'
    if 'colors' in selection_results[selection_set]:
        cat_colors = selection_results[selection_set]['colors']
    else: 
        cat_colors = [CBPALETTE[2], CBPALETTE[0], CBPALETTE[1], CBPALETTE[1]]
    assert len(cat_colors) >= len(observed_point_categories), 'not enough colors'
    if 'alphas' in selection_results[selection_set]:
        cat_alphas = selection_results[selection_set]['alphas']
    else:  
        cat_alphas = [0.4, 0.3, 0.95, 0.95]
    assert len(cat_colors) >= len(observed_point_categories), 'not enough alphas'
    
    p = (ggplot(df) +
         aes('mutation_escape', 'ACE2 binding',
             color='point_category', alpha='point_category',
             size='point_area', shape='point_category', label='mutation_label') +
         geom_point() +
         geom_text(data=df.query('mutation_label.notnull()'),
                   size=8, va='top', ha='right', alpha=1, nudge_x=-0.025, nudge_y=-0.025,
                   show_legend=False,
                   # see here for adjust_text: https://stackoverflow.com/a/57717833
                   adjust_text={'avoid_text': True,
                                'avoid_points': False,
                                'avoid_self': True,
                                'expand_text': [1.1, 1.1]},
                   ) +
         facet_wrap('~ antibody_name', nrow=1, scales='free_x') +
         scale_color_manual(values=cat_colors) +
         scale_alpha_manual(values=cat_alphas) +
         scale_size_area(limits=(size_scale * min_size, size_scale * max_size),
                                 max_size=4) +
         scale_shape_manual(values=cat_shapes) +
         theme_classic() +
         scale_x_continuous(name='mutation escape fraction',
                            breaks=mizani.breaks.mpl_breaks(nbins=3),
                            
                            ) +
         scale_y_continuous(expand=(0.07, 0)) +
         guides(alpha=False, size=False,
                shape=guide_legend(override_aes={'size': 3},
                                   title='mutation type'),
                color=guide_legend(title='mutation type'),
                ) +
         theme(figure_size=(2.3 * df['antibody'].nunique(), 2.3),
               legend_position='top' if 'legend_position' not in selection_results[selection_set]
                                else selection_results[selection_set]['legend_position'])
         )
    
    # ad hoc code to put desired points on top: those with only a few of that color
    fig = p.draw()
    for child in fig.get_children():
        if 'AxesSubplot' in str(child):
            for c in child.collections:
                if len(c._linewidths) < 10:
                    c.zorder = 2
    
    plotfile = os.path.join(config['escape_selections_dir'], f"{selection_set}.pdf")
    svgfile = os.path.splitext(plotfile)[0] + '.svg'
    print(f"Saving plot to {plotfile} and {svgfile}")
    fig.savefig(plotfile, bbox_inches='tight')
    fig.savefig(svgfile, bbox_inches='tight')
    display(fig)
    plt.close(fig)
```

    
    Making plot for Crowe_selections
    Saving plot to results/escape_selections/Crowe_selections.pdf and results/escape_selections/Crowe_selections.svg



    
![png](escape_selections_files/escape_selections_15_1.png)
    



```python

```
