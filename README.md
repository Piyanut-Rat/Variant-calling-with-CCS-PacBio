# Variant-calling-with-CCS-PacBio
----
### OUTLINE
#### 1) Import data
1.1) CCS (sample)
1.2) Download reference genome
1.3) Download Benchmark
  1.3.1) Benchmark for small variant
  1.3.2) Benchmark for structural variant
1.4)  Tandem repeat annotations 
#### 2) Filter data
#### 3) Mapping
#### 4) Variation calling
4.1) Small variantion using NanoCaller \
4.2) Structural variantion using pbsv
#### 5) Accuratecy variant
5.1)  Small variantion using hap.py
5.2)  Structural variantion using truvari

----
#### 1) Import data
1.1) CCS (sample)

Cite: https://github.com/rvalieris/parallel-fastq-dump
```
$ conda install parallel-fastq-dump
```
[CCS reads are available on NCBI SRA with accession code SRX5327410](https://www.ncbi.nlm.nih.gov/sra/?term=SRX5327410)
[SraAccList_39runs_CCSPacBio.txt](https://github.com/Piyanut-Rat/Variant-calling-with-CCS-PacBio/blob/main/SraAccList_39runs_CCSPacBio.txt)
```
$ parallel-fastq-dump --sra-id {file name} –threads 16 –outdir ccs_input/ --gzip
```
[read more example for file import](https://raw.githubusercontent.com/Piyanut-Rat/Import-data-from-Humam-genome-database/main/README.md?token=ASQS4ZRP3ACBK5NQTCSK7LLBT65OC) 

1.2) Download reference genome
  * 37
cite: [GRCh37](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GRCh37_seqs_for_alignment_pipelines/)
```
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GRCh37_seqs_for_alignment_pipelines/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz \
-O -| gunzip -c > GRCh37.fa
```
* Creating the fasta index file using samtools
cite: https://anaconda.org/bioconda/samtools
```
$ conda install -c bioconda samtools
```
cite: https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
```
$ samtools faidx ref.fasta 
```
Example
```
$ samtools faidx GRCh37.fa 
```
  * 38 
cite: [GRCh38](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)));right click [open link in the new tap](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta;tab=live_object?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))) for go to link copy to wget. (You can study and
[read more information about human genome reference](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) more.

```
$ wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
```
```
$ samtools faidx Homo_sapiens_assembly38.fasta
```
1.3) Download Benchmark 
1.3.1) Benchmark for small variant
cite: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/
```
$ mkdir -p benchmark
# GRCh38
$ FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```
```
# GRCh37
$ FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/
curl ${FTPDIR}/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
```
1.3.2) Benchmark for structural variant
download GIAB sv_benchmark\
cite: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/
```
mkdir -p giab
# GRCh37
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh37_difficult_medical_gene_SV_benchmark_v0.01.bed
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh37_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh37_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz.tbi
# GRCh38
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz.tbi
```

1.4) Tandem repeat annotations 
The pbsv discover stage was run separately per chromosome with tandem repeat annotations.
cite: https://github.com/PacificBiosciences/pbsv/tree/master/annotations
```
# GRCh38
$ wget https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
# GRCh37
$ wget https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_hs37d5.trf.bed
```

#### 2) Filter data (CCS)
###### nanoplot
https://bioconda.github.io/recipes/nanoplot/README.html
```
$ conda install nanoplot
$ conda update nanoplot
$ conda install seaborn==0.10.1
```
```
$ NanoPlot -t 16 --fastq *.fastq.gz --maxlength 40000 --plots {dot,kde} -o summary-plots-dotkde
```
###### NanoFilt
```
$ for i in /home_bif2/monwipha.milin/pj_622/data/outputfull/*.fastq.gz; do \
base=$i; \
output=${base: -19}; \
gunzip -c $base | NanoFilt -q 20 -l 10000 | gzip > /home_bif2/monwipha.milin/pj_622/data/outputfull/filt_$output; \
done
```
#### 3) Mapping
## pbmm2

* This is used for GRCh38
```
pbmm2 align /home_bif2/piyanut.ra/pj_622/data/refseq/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
  ccs_map.fofn ccs_align_sorted_38.bam --preset CCS --rg '@RG\tID:myid\tSM:HG38' --num-threads 15 --sort
```  
* This is used for GRCh37
```
$ pbmm2 align /home_bif2/piyanut.ra/pj_622/data/ref/GRCh37.fa \
  ccs_map.fofn ccs_align_sorted_GRCh37.bam --preset CCS --rg '@RG\tID:myid\tSM:HG37' --num-threads 18 --sort
```

#### 4) Variation calling
4.1)  Small variantion using NanoCaller
cite: https://github.com/WGLab/NanoCaller
```
# create the working directory

$ mkdir NanoCaller_ONT_Case_Study
$ cd NanoCaller_ONT_Case_Study
```
```
# Clone NanoCaller repository, install all the packages needed for the case study

$ git clone https://github.com/WGLab/NanoCaller.git
$ conda env create -f NanoCaller/environment.yml
$ conda activate NanoCaller
$ pip install awscli
$ conda install -y -c bioconda minimap2 bedtools 
```
```
# run nanocaller.   #! Should use screen before work because run a time per 2 days
## GRCh38
python NanoCaller/scripts/NanoCaller_WGS.py \
-bam /home_bif2/piyanut.ra/pj_622/data/pbmm2/ccs_align_sorted_38.bam -ref /home_bif2/piyanut.ra/pj_622/data/refseq/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
-prefix HG002 \
-p ccs \
--sequencing pacbio \
--snp_model CCS-HG002 \
--indel_model CCS-HG002 \
-nbr_t 0.3,0.7 \
--ins_threshold=0.4 --del_threshold=0.4 -o Nanocalls_HG38 \
--exclude_bed hg38 \
-cpu 16
```
```
## GRCh37  #! should use screen; cpu 18 used time elappsed about 12 hr.
$ python NanoCaller/scripts/NanoCaller_WGS.py \
-bam /home_bif2/piyanut.ra/pj_622/data/pbmm2/ccs_align_sorted_GRCh37.bam -ref /home_bif2/piyanut.ra/pj_622/data/ref/GRCh37.fa \
-prefix HG002_HG37 \
-p ccs \
--sequencing pacbio \
--snp_model CCS-HG002 \
--indel_model CCS-HG002 \
-nbr_t 0.3,0.7 \
-o Nanocalls_hg_37 \
--exclude_bed hg19 \
-cpu 18
```

4.2) Structural variantion using pbsv

installation: https://anaconda.org/bioconda/pbsv
```
$ conda install -c bioconda pbsv 
```
```
## This is used for HG38
$ pbsv discover ccs_align_sorted_38.bam HG38.svsig.gz --tandem-repeats /home_bif2/piyanut.ra/pj_622/data/refseq/human_GRCh38_no_alt_analysis_set.trf.bed

$ pbsv call /home_bif2/piyanut.ra/pj_622/data/refseq/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
  HG38.svsig.gz HG38.pbsv.vcf --ccs -t INS,DEL

$ bgzip HG38.pbsv.vcf
$ tabix HG38.pbsv.vcf.gz
```
```  
## This is used for HG37
$ pbsv discover ccs_align_sorted_37.bam HG37.svsig.gz --tandem-repeats /home_bif2/piyanut.ra/pj_622/data/ref/human_hs37d5.trf.bed

$ pbsv call /home_bif2/piyanut.ra/pj_622/data/ref/human_hs37d5.fasta HG37.svsig.gz HG37.pbsv.vcf --ccs -t INS,DEL

$ bgzip HG37.pbsv.vcf
$ tabix HG37.pbsv.vcf.gz
```
#### 5) Accuratecy variant
5.1)  Small variantion using hap.py \
Using cite: https://github.com/illumina/hap.py
Installed cite: https://anaconda.org/bioconda/hap.py
```
$ conda create --name hap.py
$ conda activate hap.py
$ conda install -c bioconda hap.py
```
```
$ hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa

```
```
## GRCh38
$ hap.py /home_bif2/piyanut.ra/pj_622/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf \
/home_bif2/piyanut.ra/pj_622/NanoCaller_ONT_Case_Study/Nanocalls_HG38/HG002.final.vcf.gz \
-f /home_bif2/piyanut.ra/pj_622/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-o output_prefix_hap_py_HG38 \
-r /home_bif2/piyanut.ra/pj_622/data/refseq/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta
```
```
## GRCh37
$ hap.py /home_bif2/piyanut.ra/pj_622/data/benchmark/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf \
/home_bif2/piyanut.ra/pj_622/NanoCaller_ONT_Case_Study/Nanocalls_hg_37/HG002_HG37.final.vcf.gz \
-f /home_bif2/piyanut.ra/pj_622/data/benchmark/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed \
-o output_prefix_hap_py_HG37 \
-r /home_bif2/piyanut.ra/pj_622/data/ref/GRCh37.fa
```
5.2)  Structural variantion using truvari
Installation: https://github.com/spiralgenetics/truvari
Using cite: https://github.com/PacificBiosciences/sv-benchmark
```
$ conda deactivate
$ python3 -m pip install Truvari # should base and pyth #not screen
```
```
# GRCh38
$ truvari bench -f /home_bif2/piyanut.ra/pj_622/data/refseq/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
 -b /home_bif2/piyanut.ra/pj_622/data/giab/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz \
 --includebed /home_bif2/piyanut.ra/pj_622/data/giab/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed \
 -o /home_bif2/piyanut.ra/pj_622/data/refseq/sv2_bench_pbsv38nanoFlit --passonly \
 --giabreport -r 1000 -p 0.01 --multimatch -c /home_bif2/piyanut.ra/pj_622/data/pbmm2/HG38.pbsv.vcf.gz
 # GRCh37
 $ truvari bench -f /home_bif2/piyanut.ra/pj_622/data/ref/human_hs37d5.fasta \
 -b /home_bif2/piyanut.ra/pj_622/data/giab/HG002_GRCh37_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz \
 --includebed /home_bif2/piyanut.ra/pj_622/data/giab/HG002_GRCh37_difficult_medical_gene_SV_benchmark_v0.01.bed \
 -o /home_bif2/piyanut.ra/pj_622/data/refseq/sv_bench_pbsv37nanoFlit --passonly \
 --giabreport -r 1000 -p 0.01 --multimatch -c /home_bif2/piyanut.ra/pj_622/data/pbmm2/HG37.pbsv.vcf.gz
```
Read more: https://github.com/spiralgenetics/truvari/wiki/bench
