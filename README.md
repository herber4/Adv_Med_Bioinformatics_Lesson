# Adv_Med_Bioinformatics_Lesson
This repo contains a lesson plan created for as an assignment for Advanced Medical Bioinformatics @ Clemson 040623

# Isoform expression analysis with STAR-RSEM-DESeq2 pipeline

## Outline of pipeline

1. fetch reads sratoolkit fastq-dump
2. generate genome index with star
3. align reads with star
4. generate RSEM genome index
5. quantify isoform level expression with RSEM
6. combine expression files into GEM
7. log2 transform and quantile normalize GEM for building GCN (Not applicable here)
7. use isoform expression files for DESeq2

### 1. here we will use fastq-dump to download our example data set from ncbi

fastq-dump SRR1660308
fastq-dump SRR1660309
fastq-dump SRR1660310
fastq-dump SRR1660311
fastq-dump SRR1660312
fastq-dump SRR1660313
fastq-dump SRR1660314
fastq-dump SRR1660315
fastq-dump SRR1660316
fastq-dump SRR1660317
fastq-dump SRR1660318
fastq-dump SRR1660319
fastq-dump SRR1660320
fastq-dump SRR1660321
fastq-dump SRR1660322
fastq-dump SRR1660323
fastq-dump SRR1660324

### download ref sequences including genome sequence, annotation, and transcript sequences

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz


### 2. generate genome index. This command generates an index for star to use. It takes in a primary assembly and an annotation file. create a star index in its own directory, do not create one in the working dir


STAR --runThreadN 24 --runMode genomeGenerate \
--genomeDir star/ \
--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v39.annotation.gtf \
--sjdbOverhang 99


### 3. align reads using star, this command aligns paired end reads using star

module load star/2.7.10a

for i in *_1.fastq; do name=$(basename ${i} _1.fastq);

STAR --runThreadN 24 --runMode alignReads \
--outSAMtype BAM Unsorted \
--genomeDir /data2/lackey_lab/austin/dbs/star/ \
--outFileNamePrefix ${name}_ \
--readFilesIn ${name}_1.fastq ${name}_2.fastq \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--quantMode TranscriptomeSAM \
--alignSJoverhangMin 1 \
--alignSJDBoverhangMin 8 ;
done


### 4. prepare the reference for rsem to use. This is necessary, though it seems like rsem just uses star to create an index, it actually creates its own rsem specific index


rsem-prepare-reference -p 24 --star \
--gtf gencode.v39.annotation.gtf \
GRCh38.primary_assembly.genome.fa \
rsem/


### 5. rsem command structure

-p specificies a number of threads to use [8, 16, 24, 32].
The input order is the transcriptome bam from STAR, a path to the rsem DB, and a prefix you wish the output to have

for b in *_Aligned.toTranscriptome.out.bam; do name=$(basename ${b} _Aligned.toT
ranscriptome.out.bam)
rsem-calculate-expression --alignments --paired-end -p 24 --no-bam-output --appe
nd-names ${b} /data2/lackey_lab/austin/dbs/rsem/ ${name} ;
done

### 6. rsem will output sample.isoform.results and sample.genes.results for each sample. These contain gene abundances, gene identifiers

We will use an RSEM script to combine the isoform or gene results from all samples into corresponding gene expression matrix

rsem-generate-data-matrix *.isoforms.results > rsem_isoform.results
rsem-generate-data-matrix *.genes.results > rsem_gene.results

### 7. now, proceeding with the GEM created from the sample.genes.results from all samples, we will use GEMprep python scripts to log2 transform and quantile this gene expression matrix

python ../gemprep/GEMprep/bin/normalize.py rsem_genes.txt rsem_log2.txt --log2
python ../gemprep/GEMprep/bin/normalize.py rsem_log2.txt rsem_log_quant.txt --quantile

### 8. follow commands in DESeq2 for conducting differential expression analyses

Items you will need:
1. dir with EITHER all samples .isoform.results OR .genes.results file
2. a df for meta data and setting conditions for analysis, sample names must match that of the SAMPLE.isoform.results file

| Sample | Condition |
|-----:|---------------|
|     Sample1|Mutant|
|     Sample2|Mutant|
|     Sample3|WT|
|     Sample4|WT|
