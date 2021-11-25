# Variant-calling-with-CCS-PacBio
----
### OUTLINE
#### 1) Import data
#### 2) Filter data
#### 3) Mapping
#### 4) Variation calling
3.1)  Small variantion using Deepvariant \
3.2) Structural variantion using pbsv
#### 5) Accuratecy variant
----
#### 1) Import data

Cite: https://github.com/rvalieris/parallel-fastq-dump
```
$ conda install parallel-fastq-dump
```
https://www.ncbi.nlm.nih.gov/sra/?term=SRX5327410
[SraAccList_39runs_CCSPacBio.txt](https://github.com/Piyanut-Rat/Variant-calling-with-CCS-PacBio/blob/main/SraAccList_39runs_CCSPacBio.txt)
```
$ parallel-fastq-dump --sra-id {file name} –threads 16 –outdir ccs_input/ --gzip
```
