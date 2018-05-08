
This folder shows how I made the Waterfall plots using VEP, vcftools, and GeneVisR.

GeenVisR needs one VEP output file per individual.

Input files are hosted at http://www.brassicagenome.net/databases.php

Using blastn I made input tables like this:

```
comp7993:BoRSdcaps2-10
False   Bo2g127270.1    CN      core
False   Bo2g127290.1    CN      core
False   Bo2g126980.1    NBS     core
False   Bo2g061100.1    NBS     core
False   Bo2g079110.1    NBS     variable
False   Bo2g131620.1    NL      core
False   Bo2g126860.1    NL      core
False   Bo2g127000.1    OTHER   core
False   Bo2g079150.1    TN      core
False   Bo2g079130.1    TN      core
False   Bo2g131540.1    TN      core
````
Where the first line is the name of the QTL (in terms of markers) and the other lines shows which RGA candidates are within that QTL, whether they are in another QTL (False/True), their class, and whether they're variable or core.

A simple script to make one folder for each QTL:

```python
import os
fh = open('All_QTL_nonrepetitive_Overlap_gene_names.csv')
for line in fh:
        if ':' in line:
                foldername = line.rstrip().replace(":",'_')
                os.popen('mkdir %s'%foldername)
                out = foldername + '/' + 'R_gene_names.txt'
                out = open(out, 'w')
        else:
                ll = line.split()
                name = ll[1].split('.')[0]
                out.write('%s\n'%name)
```

and a bash script to subset the SNPs and annotation for each folder:

```bash
for l in $(find . -name 'R_gene_names.txt');
do
    dirs=`dirname $l`
    cd $dirs
    grep -f R_gene_names.txt ../Brassica_oleracea.v2.1.38.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c  > Brassica_oleracea.v2.1.38_QTLs_sorted.gff3.gz
    grep -f R_gene_names.txt ../Brassica_oleracea.v2.1.34.bed > This_Brassica_oleracea.v2.1.34.bed
    vcftools --vcf ../BOLEPan.snps.13062016.vcf --bed This_Brassica_oleracea.v2.1.34.bed --out BOLEPan.snps.13062016_only_R.vcf --recode --keep-INFO-all
    cd ..
done
```

Then another Python script to subsample the vcf files (bcftools view -s should also work):

```python
from glob import glob
import os
for x in glob('*/BOLEPan.snps.13062016_only_R.vcf.recode.vcf'):
        dirname = os.path.dirname(x)
        header = []
        outs = []
        fh = open(x)
        for line  in fh:
                ll = line.split()
                if line.startswith('#'):
                        if ll[0] != '#CHROM':
                                header.append(line)
                        else:
                                # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Broccoli        Brussels        Cabbage1        Cabbage2        Cauliflower1    Cauliflower2    Kale    Kohlrabi        Macrocarpa      TO1000

                                snp_start = ll.index('FORMAT')+1
                                inds = ll[ll.index('FORMAT')+1:]
                                for i in inds:
                                        this_fh = open(fh.name + '_' + i + '.vcf', 'w')
                                        outs.append( this_fh )
                                        for h in header:
                                                this_fh.write(h)
                                        this_fh.write('\t'.join( ll[:ll.index('FORMAT')+1] + [i]) + '\n')
                        continue
                alleles = ll[snp_start:]
                start_lines = ll[:snp_start]
                for i, out in zip(alleles, outs):
                        thisll = list(start_lines)# make a copy
                        thisll.append(i)
                        out.write('\t'.join(thisll) + '\n')
```

Then to filter those resulting vcfs using bcftools:

```bash
for l in */*recode.vcf_*vcf; do bcftools view --trim-alt-alleles $l > ${l}_trimmedalt.vcf; done
for l in */*trimmedalt.vcf; do bcftools view -m 2 $l > ${l}snps_only.vcf; done
```

Now to run VEP, it's important that at the time of this writing, if you use a VEP version other than 88 GenVisR will crash.

```bash
   wget https://github.com/Ensembl/ensembl-vep/archive/release/88.13.zip
   unzip 88.13.zip
   cd ensembl-vep-release-88.13
   perl INSTALL.pl
```

and to run VEP:

```bash
for l in */*vcf_*snps_only.vcf; do
   out=${l/BOLEPan.snps.13062016_only_R.vcf.recode.vcf_/}
   out=${out/.vcf_trimmedalt.vcfsnps_only.vcf/}.vep
   dirs=`dirname $l`;
   ./ensembl-vep-release-88.13/vep -i $l -v -gff ${dirs}/Brassica_oleracea.v2.1.38_QTLs_sorted.gff3.gz -fasta Brassica_oleracea.v2.1.dna.toplevel.fa.gz -species boleracea --output_file ${out} --fork 8 --force_overwrite --no_intergenic
done
```

Now the annoying thing is that GenVisR expects us to have a SYMBOL column everywhere, since apparently if you run VEP with human data you get that. So let's add that:

```python
from glob import glob

for x in glob('*/*.vep'):

        out = open(x + 'fixed', 'w')
        for line in open(x):
                if line.startswith('#'):
                        out.write(line)
                        continue
                ll = line.split()
                if not ll: continue
                print(ll)
                line = line.rstrip()
                line += ';SYMBOL=%s'%(ll[3])
                out.write(line + '\n')
```

All this does is take the gene/protein name as SYMBOL.

Now rename the output files to overwrite the input if everything worked:

    for l in */*vepfixed; do mv $l ${l//fixed}; done

Now I'm adding the PAV status to each VEP file:

```python
pav_dict = {} # key: gene, value: list of LOSS individuals

with open('BOLEPan.pav.13062016.vcf') as fh:
        for line in fh:
                ll = line.split()
                if ll[0].startswith('#CHROM'):
                        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Broccoli        Brussels        Cabbage1        Cabbage2        Cauliflower1    Cauliflower2    Kale    Kohohlrabi      TO1000  Macrocarpa
                        start_pos = ll.index('Broccoli')
                        names = ll[start_pos:]
                        continue
                gene = ll[2]

                lost_inds = set()
                alleles = ll[start_pos:]
                for allele, ind in zip(alleles, names):
                        if allele == '0/0':
                                lost_inds.add(ind)
                if lost_inds:
                        pav_dict[gene] = lost_inds

from glob import glob

for f in glob('*/*.vep'):
        # BOLEPan.snps.13062016_only_R.vcf.recode.vcf_Cabbage1.vcf_trimmedalt.vcfsnps_only.vcf_variant_effect_outputfixed
        print(f)
        thisind = f.split('/')[-1].replace('.vep','')
        print(thisind)
        out = open(f + '_wih_pav.vep', 'w')
        added = set()
        with open(f) as fh:
                for line in fh:
                        if line.startswith('#'):
                                out.write(line)
                                continue
                        ll = line.split()
                        # .       C1:26272519     A       Bo1g087950      Bo1g087950.1    Transcript      intron_variant  -       -       -       -       -       -       IMPACT=MODIFIER;STRAND=1;SOURCE=Brassica_oleracea.v2.1.38_QTLs_sorted.gff3.gz;SYMBOL=Bo1g087950
                        name = ll[3]
                        if name in pav_dict:
                                if thisind in pav_dict[name]:
                                        # add PAV
                                        if name in added:
                                                continue
                                        added.add(name)
                                        newll = ll[:6] + ['gene_lost']+ll[7:]
                                        out.write('\t'.join(newll) + '\n')
                        out.write(line)
```

and after checking move the 'wih_pav' files over:

   for l in */*_wih_pav.vep; do mv $l ${l%%_wih_pav.vep}; done


AND FINALLY we can plot:

```R
# set a seed
set.seed(426)
library(RColorBrewer)
mypal <- c('black', brewer.pal(11, 'Paired'))
tol21rainbow= c("#AA4488", "#CC99BB", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
mypal <- c('black', tol21rainbow)
mypal <- mypal[1:12]
print(mypal)
library(data.table)
mutations <- unique(getMutation(vepObject)$Conse)
getMutation(vepObject$Consequence)
mymutations <- data.table(mutation=rev(c('downstream_gene_variant', 'upstream_gene_variant', 'intron_variant', 'synonymous_variant', 'splice_region_variant',
                 'missense_variant','start_lost','stop_lost','stop_gained','splice_donor_variant',
                 'splice_acceptor_variant', 'gene_lost')), color=mypal)

# load GenVisR into R
library(GenVisR)
for (i in list.files('Plot', full.names=T)){
  if(dir.exists(i)) {
    print(paste0('Plotting in ', i))
    filesdir <- i
    paste0(filesdir, '*vep')
    my_vep <- Sys.glob(paste0(filesdir, "/*vep"))
    vepObject <- VEP(my_vep)

    png(paste0(filesdir, 'test.png'), height=480*3, width=480*3)

    drawPlot(Waterfall(vepObject,recurrence = 0.4, mutationHierarchy=mymutations ))
    dev.off()

  }
}
```

