'''
IKI Bangladesh (MIOASI): Plot boxplot summaries

Plot boxplot summaries of differences between downscaled RA2 data and IBTrACS and ERA5.

Note: Python 3 compatible only

Author: HS
Created: 26/5/20
'''

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

var = 'mslp'
# var = 'gust'

df = pd.read_csv('validation_diff.csv')
# Convert timedeltas
# df.TIMEDIFF = pd.to_timedelta(df.TIMEDIFF).dt.total_seconds()/3600
# Find median intensity order
stormorder = df[df.DIFF_TYPE == var].groupby(['STORM']).median().sort_values(by='INTDIFF').index.values
# Plot
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(7, 5), sharey=True)

for axnum, xvar in enumerate(['INTDIFF', 'TIMEDIFF']):
    sns.boxplot(data=df[df.DIFF_TYPE == var], x=xvar, y='STORM', ax=axs[axnum],
                whis=(5, 95), color='lightgrey', showfliers=False, order=stormorder, linewidth=1)
    
    sns.swarmplot(data=df[df.DIFF_TYPE == var], x=xvar, y='STORM', ax=axs[axnum],
                  hue='DIFF_WRT', size=3, order=stormorder)

    if axnum==1:
        sns.despine(left=True)


sns.despine(fig=fig, offset=10, trim=True)
axs[0].set(xlabel='Difference (hPa)')
axs[1].set(xlabel='Difference (hours)')
plt.tight_layout()
# plt.show()
fig.savefig(f'compare_{var}.eps')

# Data summaries can also be found via pandas, eg.
# df[df.DIFF_TYPE == 'mslp'].groupby(['STORM']).median().sort_values(by=['INTDIFF'])
# or df[df.DIFF_TYPE == 'gust'].groupby(['STORM', 'DIFF_WRT']).max().sort_values(by=['DIFF_WRT', 'INTDIFF'])
# etc.
