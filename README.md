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
c) paste your illumina data in `data` folder and your genome files,
   such as fasta and gff3 (only for RNAseq data) files in `genome` folder

Now, your folder structure will be as follows:
- Bulk2SNPs
  - data
    - Bulk1_1.fastq
    - Bulk1_2.fastq
    - Bulk2_1.fastq
    - Bulk2_2.fastq
  - genome
    - genome.fasta
    - annotation.gff
```

#### STEP II: Install docker and keep it running in the background
```md
Download the docker software based on your operating system from
https://www.docker.com
```

#### STEP III: Clone the Bulk2SNPs repo inside the folder created in STEP I 
```bash
git clone https://github.com/NDSUrustlab/Bulk2SNPs.git
```

#### STEP IV: Run the pipeline 
```bash
nextflow run main.nf \
  --bulk1 'Bulk1' \
  --bulk2 'Bulk2' \
  --genome 'genome.fasta' \
  --gff3 'annotation.gff3
```

## Here's the flowdiagram depicting pipeline's workflow
![](media/flowdiagram.png)
