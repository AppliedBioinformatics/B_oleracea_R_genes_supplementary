

def merge(lsts):
  # from https://stackoverflow.com/a/9112588
  sets = [set(lst) for lst in lsts if lst]
  merged = 1
  while merged:
    merged = 0
    results = []
    while sets:
      common, rest = sets[0], sets[1:]
      sets = []
      for x in rest:
        if x.isdisjoint(common):
          sets.append(x)
        else:
          merged = 1
          common |= x
      results.append(common)
    sets = results
  return sets

def get_names(lines, contig):
    names = []
    for l in lines:
        ll = l.split()
        chrom = ll[0]
        if chrom != contig:
            continue
        name = ll[3]
        names.append(name)
    return names


import os
# For example - runs any type of RGA
fh = open('Brassica.TMCC.candidates.lst')
genes_to_get = set()
for line in fh:
    ll = line.split()
    name = ll[0]
    genes_to_get.add(name)


cluster_counter = 0
cluster_dict = {}
# bed file generated via gff2bed
for line in open('Brassica_oleracea.v2.1.34.bed'):
    ll = line.split()
    if not ll: continue
    name = ll[3]
    this_contig = ll[0]
    #name = get_name(ll)
    if name not in genes_to_get:
        continue
    # slow solution, but bed files are small
    previous_ten_lines = os.popen('LC_ALL=C grep -B 10 "%s" Brassica_oleracea.v2.1.34.bed'%(name)).read().rstrip().split('\n')[1:]
    next_ten_lines = os.popen('LC_ALL=C grep -A 10 "%s" Brassica_oleracea.v2.1.34.bed'%(name)).read().rstrip().split('\n')[1:]
    
    previous_ten_names = get_names(previous_ten_lines, this_contig)
    next_ten_names = get_names(next_ten_lines, this_contig)
    this_cluster = set([name])
    for l in previous_ten_names:
        if l in genes_to_get:
            this_cluster.add(l)
    for l in next_ten_names:
        if l in genes_to_get:
            this_cluster.add(l)

    cluster_dict[cluster_counter] = this_cluster
    cluster_counter += 1

list_of_lists = []

for k in cluster_dict:
    list_of_lists.append(cluster_dict[k])

list_of_lists = merge(list_of_lists)

for m in list_of_lists:
    print(m)


