'''
Fisher's exact test to identify, following on footprinting analysis, the set of atac peaks with 
significant differences in motif bases bound per peak across FC (Care) and NC (NoCare) conditions.
'''

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np

# total motifs bases per peak
total = pd.read_csv('potential_motifs_bases_per_peak.bed', sep='\t', header=None, usecols=[3, 5])

#motifs bases that TOBIAS predicted bound in Care
a_bound = pd.read_csv('FC_bound_bases_per_peak.bed', sep='\t', header=None, usecols=[3, 5])

#motifs bases that TOBIAS predicted bound in noCare
b_bound = pd.read_csv('NC_bound_bases_per_peak.bed', sep='\t', header=None, usecols=[3, 5])

print(total)

df = pd.DataFrame({
    'Peak_ID': total[3],
    'Potential': total[5],
    'Bound_FC': a_bound[5],
    'Bound_NC': b_bound[5]
})

df = df[(df['Bound_FC'] > 0) | (df['Bound_NC'] > 0)].copy()

df = df[(df['Bound_FC'] != df['Bound_NC'])].copy()

# Fisher's Exact Test
def run_fisher(row):
    # Potential bases vs Bound bases
    # Condition A: [Bound, Not_Bound]
    # Condition B: [Bound, Not_Bound]
    a, b = row['Bound_FC'], row['Potential'] - row['Bound_FC']
    c, d = row['Bound_NC'], row['Potential'] - row['Bound_NC']
    
    # Return 1.0 if no binding in either to avoid errors
    if (a + c) == 0: return 1.0
    
    _, p = fisher_exact([[a, b], [c, d]])
    return p

df['p_value'] = df.apply(run_fisher, axis=1)

# Correct for multiple testing
_, df['padj'], _, _ = multipletests(df['p_value'], method='fdr_bh')

df['log2FC'] = np.log2((df['Bound_FC'] + 1) / (df['Bound_NC'] + 1))

df.to_csv('DB_peaks_all_elise.tsv', index=False, sep='\t')

significant = df[df['padj'] < 0.01]
print(significant)
care = significant[significant['log2FC']>2]
nocare = significant[significant['log2FC']<-2]

nocare.to_csv('DB_noCarepeaks.tsv', index=False, sep='\t')
care.to_csv('DB_Carepeaks.tsv', index=False, sep='\t')