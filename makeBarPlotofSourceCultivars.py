

chrom_dict = {}

gff3 = open('BOLEPan.genes.13062016.gff3')
for line in gff3:
    ll = line.split()
    # jcf7180000680165_Cabbage1       maker   gene    95      1696    .       -       .       ID=BOLEPAN_00000379;Name=BOLEPAN_00000379;Alias=maker-jcf7180000680165_Cabbage1-snap-gene-0.2;
    try:
        if ll[2] != 'gene': continue
    except:
        continue
    names = ll[-1].split(';')
    thisid = names[0].replace('ID=','')
    chrom = ll[0]
    source = chrom.split('_')[-1]
    chrom_dict[thisid] = source


import glob

from collections import Counter

import pandas as pd

all_sources = []
all_classes = []
all_counts = []

for l in glob.glob('*lst'):
    if 'RGA' in l: continue
    if 'TM' in l: continue
    c = Counter()
    fh = open(l)
    for line in fh:
        ll = line.split()
        gene = ll[0][:-1]
        c[chrom_dict[gene]] += 1
    this_class = l.split('.')[1]
    
    for x in c:
        all_sources.append(x)
        all_classes.append(this_class)
        all_counts.append(c[x])
big_table = pd.DataFrame({'Source':all_sources, 'Class':all_classes, 'Count':all_counts})
big_table.sort_values(['Class', 'Count'],axis=0,inplace=True, ascending=[True, False])

# I wan tto sort the individuals in the graph by their total number of R-genes, not by the counts of classes
grouped = big_table.groupby('Source').sum()
names = list(grouped.index)
counts = list(grouped['Count'])
classes = len(names) * ['Total']
small_table = pd.DataFrame({'Source':names, 'Class':classes, 'Count':counts})
big_table = big_table.append(small_table)
big_table.sort_values('Count', inplace=True, ascending=False)
# now delete the Total rows
my_order = big_table[big_table['Class']=='Total']['Source']
print(big_table.groupby('Class').sum())
big_table = big_table[big_table['Class'] != 'Total']

import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns
sns.set(style="ticks")
#sns.swarmplot(x='Class',y='Count', hue='Source', data = big_table)
print(big_table)
p = sns.barplot(hue='Class',y='Count', x='Source', data = big_table, palette='colorblind', order=my_order)
plt.xticks(rotation=45)

sns.despine()
p.legend(loc='upper right')
#sns.boxplot(x='Class',y='Count', hue='Source', data = big_table)
plt.tight_layout()
plt.savefig('Figure_Barplot.png')
plt.savefig('Figure_Barplot.eps', format='eps', dpi=1000)
