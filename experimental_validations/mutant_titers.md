# Titer analysis for mutant viruses

RLU titering and p24 ELISAs for all viruses made for this study that grew well in initial titering. 

This includes the double mutant viruses that were not actually run in neuts. 

Most of these viruses were grown by Andrea Loes. The December viruses were grown by Amin Addetia and the WT virus for neuts (September) was grown by Kate Crawford.

Initial titering was done by whoever grew the virus. Kate Crawford did final titering and p24 ELISAs.

Each virus was titered or run in the p24 ELISA in duplicate. Values were averaged prior to plotting and a single point is shown for each virus. 

Raw data is in the January 22 and 25, 2021 excel files in the `./data` subdirectory. Those excel files include titer and p24 calculations.


```python
import os
import warnings

from IPython.display import display, HTML
from IPython.display import display, SVG
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from plotnine import *

```


```python
warnings.simplefilter('ignore')
```


```python
CBP = ('#999999', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')
plt.style.use('seaborn-white')
theme_set(theme_seaborn(style='white', context='talk', font='FreeSans', font_scale=1))
```


```python
resultsdir='results/p24_titering/'
os.makedirs(resultsdir, exist_ok=True)
```

## Read in csv of cleaned titer and p24 data


```python
df = pd.read_csv('./data/210125_p24_titers_calc.csv').set_index('VirusNumber', drop=True)
display(df)
```


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
      <th>Sample</th>
      <th>VirusDate</th>
      <th>Used in neuts</th>
      <th>p24 pg/mL</th>
      <th>Avg RLU/mL</th>
      <th>Median RLU/mL</th>
      <th>TechRep</th>
      <th>Neut Dilution</th>
    </tr>
    <tr>
      <th>VirusNumber</th>
      <th></th>
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
      <th>1</th>
      <td>WT</td>
      <td>20-10-30</td>
      <td>N</td>
      <td>8.889640e+05</td>
      <td>401000000.0</td>
      <td>384134400</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>P384L</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.085360e+06</td>
      <td>209000000.0</td>
      <td>189574400</td>
      <td>1</td>
      <td>20.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>F456K</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.302027e+06</td>
      <td>304000000.0</td>
      <td>299872000</td>
      <td>1</td>
      <td>10.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>E484Q</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.022748e+06</td>
      <td>72900000.0</td>
      <td>76072000</td>
      <td>1</td>
      <td>20.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>G485R</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>7.605856e+05</td>
      <td>299541114.3</td>
      <td>294278200</td>
      <td>1</td>
      <td>32.0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>S494P</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.051577e+06</td>
      <td>106000000.0</td>
      <td>107374200</td>
      <td>1</td>
      <td>50.0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>G446V/F456V</td>
      <td>20-10-30</td>
      <td>N</td>
      <td>7.686937e+05</td>
      <td>39600000.0</td>
      <td>43172400</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>G446V/F456V/E484P</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>8.141892e+05</td>
      <td>48100000.0</td>
      <td>49775200</td>
      <td>1</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>WT</td>
      <td>20-10-16</td>
      <td>N</td>
      <td>1.787613e+06</td>
      <td>416523628.6</td>
      <td>425927200</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>10</th>
      <td>F456V</td>
      <td>20-10-16</td>
      <td>Y</td>
      <td>8.975225e+05</td>
      <td>478000000.0</td>
      <td>486155200</td>
      <td>1</td>
      <td>40.0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>E484P</td>
      <td>20-10-16</td>
      <td>Y</td>
      <td>4.173423e+05</td>
      <td>38400000.0</td>
      <td>38344000</td>
      <td>1</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>WT</td>
      <td>20-11-9</td>
      <td>N</td>
      <td>9.385135e+05</td>
      <td>383000000.0</td>
      <td>380955200</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>13</th>
      <td>G446V</td>
      <td>20-11-9</td>
      <td>Y</td>
      <td>8.425675e+05</td>
      <td>170000000.0</td>
      <td>166024000</td>
      <td>1</td>
      <td>25.0</td>
    </tr>
    <tr>
      <th>14</th>
      <td>E484P</td>
      <td>20-11-9</td>
      <td>Y</td>
      <td>5.024774e+05</td>
      <td>4120000.0</td>
      <td>4851200</td>
      <td>1</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>F456V/E484P</td>
      <td>20-11-9</td>
      <td>N</td>
      <td>6.569820e+05</td>
      <td>29500000.0</td>
      <td>29219200</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>16</th>
      <td>WT</td>
      <td>20-12-5</td>
      <td>N</td>
      <td>7.506757e+05</td>
      <td>181000000.0</td>
      <td>169088800</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>17</th>
      <td>C432D</td>
      <td>20-12-5</td>
      <td>N</td>
      <td>3.362613e+05</td>
      <td>32600.0</td>
      <td>3200</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>18</th>
      <td>E484K</td>
      <td>20-12-5</td>
      <td>Y</td>
      <td>6.479730e+05</td>
      <td>329000000.0</td>
      <td>327698000</td>
      <td>1</td>
      <td>50.0</td>
    </tr>
    <tr>
      <th>19</th>
      <td>WT</td>
      <td>20-11-13</td>
      <td>N</td>
      <td>1.681306e+06</td>
      <td>646000000.0</td>
      <td>601462400</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>20</th>
      <td>F456A</td>
      <td>20-11-13</td>
      <td>Y</td>
      <td>9.060810e+05</td>
      <td>32900000.0</td>
      <td>41859200</td>
      <td>1</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>21</th>
      <td>F456A/E484P</td>
      <td>20-11-13</td>
      <td>N</td>
      <td>1.334459e+06</td>
      <td>80400000.0</td>
      <td>79824000</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>22</th>
      <td>WT</td>
      <td>20-09-20</td>
      <td>Y</td>
      <td>9.632883e+05</td>
      <td>229000000.0</td>
      <td>195571600</td>
      <td>1</td>
      <td>50.0</td>
    </tr>
    <tr>
      <th>23</th>
      <td>NoVEP</td>
      <td>20-11-9</td>
      <td>N</td>
      <td>2.837162e+06</td>
      <td>3200.0</td>
      <td>800</td>
      <td>1</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>WT</td>
      <td>20-10-30</td>
      <td>N</td>
      <td>8.889640e+05</td>
      <td>372000000.0</td>
      <td>373344200</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>P384L</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.085360e+06</td>
      <td>193000000.0</td>
      <td>189942400</td>
      <td>2</td>
      <td>20.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>F456K</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.302027e+06</td>
      <td>207000000.0</td>
      <td>208743600</td>
      <td>2</td>
      <td>10.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>E484Q</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.022748e+06</td>
      <td>52900000.0</td>
      <td>49203200</td>
      <td>2</td>
      <td>20.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>G485R</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>7.605856e+05</td>
      <td>222040271.4</td>
      <td>230166400</td>
      <td>2</td>
      <td>32.0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>S494P</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.051577e+06</td>
      <td>119000000.0</td>
      <td>111993000</td>
      <td>2</td>
      <td>50.0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>G446V/F456V</td>
      <td>20-10-30</td>
      <td>N</td>
      <td>7.686937e+05</td>
      <td>48600000.0</td>
      <td>42281200</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>G446V/F456V/E484P</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>8.141892e+05</td>
      <td>38900000.0</td>
      <td>39028000</td>
      <td>2</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>WT</td>
      <td>20-10-16</td>
      <td>N</td>
      <td>1.787613e+06</td>
      <td>400264957.1</td>
      <td>382692000</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>10</th>
      <td>F456V</td>
      <td>20-10-16</td>
      <td>Y</td>
      <td>8.975225e+05</td>
      <td>475000000.0</td>
      <td>479437700</td>
      <td>2</td>
      <td>40.0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>E484P</td>
      <td>20-10-16</td>
      <td>Y</td>
      <td>4.173423e+05</td>
      <td>46200000.0</td>
      <td>43204200</td>
      <td>2</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>WT</td>
      <td>20-11-9</td>
      <td>N</td>
      <td>9.385135e+05</td>
      <td>400000000.0</td>
      <td>383332800</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>13</th>
      <td>G446V</td>
      <td>20-11-9</td>
      <td>Y</td>
      <td>8.425675e+05</td>
      <td>161000000.0</td>
      <td>145967100</td>
      <td>2</td>
      <td>25.0</td>
    </tr>
    <tr>
      <th>14</th>
      <td>E484P</td>
      <td>20-11-9</td>
      <td>Y</td>
      <td>5.024774e+05</td>
      <td>8240000.0</td>
      <td>6706600</td>
      <td>2</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>F456V/E484P</td>
      <td>20-11-9</td>
      <td>N</td>
      <td>6.569820e+05</td>
      <td>32500000.0</td>
      <td>31820800</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>16</th>
      <td>WT</td>
      <td>20-12-5</td>
      <td>N</td>
      <td>7.506757e+05</td>
      <td>200000000.0</td>
      <td>179872000</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>17</th>
      <td>C432D</td>
      <td>20-12-5</td>
      <td>N</td>
      <td>3.362613e+05</td>
      <td>3990.0</td>
      <td>1200</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>18</th>
      <td>E484K</td>
      <td>20-12-5</td>
      <td>Y</td>
      <td>6.479730e+05</td>
      <td>372000000.0</td>
      <td>355718400</td>
      <td>2</td>
      <td>50.0</td>
    </tr>
    <tr>
      <th>19</th>
      <td>WT</td>
      <td>20-11-13</td>
      <td>N</td>
      <td>1.681306e+06</td>
      <td>634000000.0</td>
      <td>574771200</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>20</th>
      <td>F456A</td>
      <td>20-11-13</td>
      <td>Y</td>
      <td>9.060810e+05</td>
      <td>52700000.0</td>
      <td>43949400</td>
      <td>2</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>21</th>
      <td>F456A/E484P</td>
      <td>20-11-13</td>
      <td>N</td>
      <td>1.334459e+06</td>
      <td>86300000.0</td>
      <td>84700100</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>22</th>
      <td>WT</td>
      <td>20-09-20</td>
      <td>Y</td>
      <td>9.632883e+05</td>
      <td>223000000.0</td>
      <td>222003200</td>
      <td>2</td>
      <td>50.0</td>
    </tr>
    <tr>
      <th>23</th>
      <td>NoVEP</td>
      <td>20-11-9</td>
      <td>N</td>
      <td>2.837162e+06</td>
      <td>6300.0</td>
      <td>800</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>


## Add columns for avg and median RLUs normalized by p24


```python
# add columns for Avg and median RLUs normalized by p24
df['Avg RLU/pg p24'] = df['Avg RLU/mL'] / df['p24 pg/mL']
df['Med RLU/pg p24'] = df['Median RLU/mL'] / df['p24 pg/mL']
display(df.head())
```


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
      <th>Sample</th>
      <th>VirusDate</th>
      <th>Used in neuts</th>
      <th>p24 pg/mL</th>
      <th>Avg RLU/mL</th>
      <th>Median RLU/mL</th>
      <th>TechRep</th>
      <th>Neut Dilution</th>
      <th>Avg RLU/pg p24</th>
      <th>Med RLU/pg p24</th>
    </tr>
    <tr>
      <th>VirusNumber</th>
      <th></th>
      <th></th>
      <th></th>
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
      <th>1</th>
      <td>WT</td>
      <td>20-10-30</td>
      <td>N</td>
      <td>8.889640e+05</td>
      <td>401000000.0</td>
      <td>384134400</td>
      <td>1</td>
      <td>NaN</td>
      <td>451.086900</td>
      <td>432.114703</td>
    </tr>
    <tr>
      <th>2</th>
      <td>P384L</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.085360e+06</td>
      <td>209000000.0</td>
      <td>189574400</td>
      <td>1</td>
      <td>20.0</td>
      <td>192.562766</td>
      <td>174.664932</td>
    </tr>
    <tr>
      <th>3</th>
      <td>F456K</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.302027e+06</td>
      <td>304000000.0</td>
      <td>299872000</td>
      <td>1</td>
      <td>10.0</td>
      <td>233.482106</td>
      <td>230.311665</td>
    </tr>
    <tr>
      <th>4</th>
      <td>E484Q</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>1.022748e+06</td>
      <td>72900000.0</td>
      <td>76072000</td>
      <td>1</td>
      <td>20.0</td>
      <td>71.278574</td>
      <td>74.380023</td>
    </tr>
    <tr>
      <th>5</th>
      <td>G485R</td>
      <td>20-10-30</td>
      <td>Y</td>
      <td>7.605856e+05</td>
      <td>299541114.3</td>
      <td>294278200</td>
      <td>1</td>
      <td>32.0</td>
      <td>393.829574</td>
      <td>386.910018</td>
    </tr>
  </tbody>
</table>
</div>


## Order samples by spike mutant

Order is NoVEP/None, WT, single mutants, double mutants, triple mutant.


```python
# change 'NoVEP' to 'None'
df['Sample'] = df['Sample'].replace('NoVEP', 'None')

# order muts in a logical way
mut_list = ['None', 'WT', 'C432D', 'E484K', 'E484P', 'E484Q', 'F456A', 'F456K',
            'F456V', 'G446V', 'G485R', 'P384L', 'S494P', 'F456A/E484P', 'F456V/E484P', 'G446V/F456V',
            'G446V/F456V/E484P']
mut_cat = pd.Categorical(df['Sample'], categories=mut_list)
df = df.assign(mut_order = mut_cat)

# add linebreak for triple mutant
df['mut_order'] = df['mut_order'].replace('G446V/F456V/E484P', 'G446V/F456V/\nE484P')
```

## Plots

### Plot avg RLUs/pg p24


```python
avgrlu_p24_all = (ggplot(df, aes(x='mut_order', y='Avg RLU/pg p24', fill='Used in neuts')) +
                  geom_point(size=4, alpha=0.8, position=position_jitterdodge(jitter_width=0.2, jitter_height=0, dodge_width=0.2, random_state=123)) +
#             geom_crossbar(data=luc_titers.groupby('Spike', as_index=False).aggregate(Titer=pd.NamedAgg('Avg RLU per mL', 'mean')),
#                           mapping=aes(x='Spike', y='Titer', ymin='Titer', ymax='Titer'), color='black') +
                  scale_fill_manual(values=[CBP[0], CBP[0]], guide=False) +
                  theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                        figure_size=(10,5)) +
                  geom_hline(yintercept=np.mean(df[df['mut_order']=='WT']['Avg RLU/pg p24']), 
                        linetype='dashed', color='grey') +
                  scale_y_continuous(trans='log10', limits=[0.0001, 1e3]) +
                  ylab('relative luciferase units per pg p24') + 
                  xlab('spike variant')
            )

_ = avgrlu_p24_all.draw()
```


![png](mutant_titers_files/mutant_titers_12_0.png)


### Plot median RLUs/pg p24


```python
medrlu_p24_all = (ggplot(df, aes(x='mut_order', y='Med RLU/pg p24', fill='Used in neuts')) +
            geom_point(size=4, alpha=0.8, position=position_jitterdodge(jitter_width=0.2, jitter_height=0, dodge_width=0.2, random_state=123)) +
#             geom_crossbar(data=luc_titers.groupby('Spike', as_index=False).aggregate(Titer=pd.NamedAgg('Avg RLU per mL', 'mean')),
#                           mapping=aes(x='Spike', y='Titer', ymin='Titer', ymax='Titer'), color='black') +
            scale_fill_manual(values=[CBP[0], CBP[0]], guide=False) +
            theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                  figure_size=(10,5)) +
            geom_hline(yintercept=np.median(df[df['mut_order']=='WT']['Med RLU/pg p24']), 
                       linetype='dashed', color='grey') +
            scale_y_continuous(trans='log10', limits=[0.0001, 1e3]) +
            ylab('relative luciferase units per pg p24') + 
            xlab('spike variant')
            )

_ = medrlu_p24_all.draw()

plotfile = f'{resultsdir}/median_p24_norm_titers.pdf'
print(f"Saving to {plotfile}")
medrlu_p24_all.save(plotfile, verbose=False)
```

    Saving to results/p24_titering//median_p24_norm_titers.pdf



![png](mutant_titers_files/mutant_titers_14_1.png)


### Plot median RLUs/pg p24 without showing double mutants


```python
doubles = ['F456A/E484P', 'F456V/E484P', 'G446V/F456V']

medrlu_p24_nodoubles = (ggplot(df[~df['mut_order'].isin(doubles)], aes(x='mut_order', y='Med RLU/pg p24', fill='Used in neuts')) +
            geom_point(size=4,alpha=0.8, position=position_jitterdodge(jitter_width=0.2, jitter_height=0, dodge_width=0.2, random_state=123)) +
#             geom_crossbar(data=luc_titers.groupby('Spike', as_index=False).aggregate(Titer=pd.NamedAgg('Avg RLU per mL', 'mean')),
#                           mapping=aes(x='Spike', y='Titer', ymin='Titer', ymax='Titer'), color='black') +
            scale_fill_manual(values=[CBP[0], CBP[0]], guide=False) +
            theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                  figure_size=(8,3)) +
            geom_hline(yintercept=np.median(df[df['mut_order']=='WT']['Med RLU/pg p24']), 
                       linetype='dashed', color='grey') +
            scale_y_continuous(trans='log10', limits=[0.0001, 1e3]) +
            ylab('relative luciferase units\nper pg p24') +
            xlab('spike variant')
            )

_ = medrlu_p24_nodoubles.draw()

plotfile = f'{resultsdir}/median_p24_norm_titers_nodoubles.pdf'
print(f"Saving to {plotfile}")
medrlu_p24_nodoubles.save(plotfile, verbose=False)
```

    Saving to results/p24_titering//median_p24_norm_titers_nodoubles.pdf



![png](mutant_titers_files/mutant_titers_16_1.png)


## Convert to markdown


```python
!jupyter nbconvert mutant_titers.ipynb --to markdown
```

    [NbConvertApp] Converting notebook mutant_titers.ipynb to markdown
    [NbConvertApp] Support files will be in mutant_titers_files/
    [NbConvertApp] Making directory mutant_titers_files
    [NbConvertApp] Making directory mutant_titers_files
    [NbConvertApp] Making directory mutant_titers_files
    [NbConvertApp] Writing 19314 bytes to mutant_titers.md



```python

```
