![](media/logo.png)
---
Bulk2SNPs is an automated pipeline to speed up SNP discovery from next-generation sequencing data, i.e., DNA and RNA, generated from contrasting bulks to map quantitative trait loci (QTLs) for desired traits. 

## USAGE: 
#### STEP I: Prepare your data
```md
a) make a new folder for your analysis. E.g. Bulk2SNPs
b) Inside this folder create two more folder:
  i) data
  ii) genome or index (if you have index already generated)
c) paste your illumina data in `data` folder and your genome files, such as fasta and gff3 (only for RNAseq data) files in `genome` folder
```


### Example usage:
```bash
nextflow run main.nf \
  --bulk1 'Bulk1Name' \
  --bulk2 'Bulk2Name' \
  --genome 'genome.fasta' \
  --gff3 'annotation.gff3
```

## Here's the flowdiagram depicting pipeline's workflow
![](media/flowdiagram.png)
