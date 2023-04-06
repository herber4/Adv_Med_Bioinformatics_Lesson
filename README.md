# Adv_Med_Bioinformatics_Lesson
This repo contains a lesson plan created as an assignment for Advanced Medical Bioinformatics @ Clemson 040623

# Isoform expression analysis with STAR-RSEM-DESeq2 pipeline, lab overview

Here we are interested in myeolodysplastic syndrome cases with SF3B1 K700E mutations. SF3B1 is a core componenet of the U2 snRNP complex in the spliceosome and is involved in branch point identification and 3' splice site selection. It is known that mutation of SF3B1 results in global increases to alternative splicing. So, how can splicing factor mutations impact gene expression? 

1. Alternative splicing, through the use of different splice sites, can directly result in differential expression of transcript isoforms.
2. Alternative splicing may inactivate/activate transcription factors which regulate expression of groups of genes. 

To test this hypothesis, we will utilize bulk RNA-sequencing data from patient CD34+ B cell precursors diagnosed with MDS and positive or negative for the SF3B1 K700E mutation. With this data, we will perform transcriptome mapping (STAR) and quantification (RSEM) and use DESeq2 to identify differentially expressed genes between WT and Mutant samples. 

GEO link to data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63569
Link to data source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5029572/

## Outline of pipeline, learning objectives

1. fetch reads sratoolkit fastq-dump
2. generate genome index with star
3. align reads with star
4. generate RSEM genome index
5. quantify isoform level expression with RSEM
6. combine expression files into GEM
7. log2 transform and quantile normalize GEM for building GCN (Not applicable here)
7. use isoform expression files for DESeq2

### 1. here we will use fastq-dump to download our example data set from ncbi

Here is one example of the fastq-dump command from the sra toolkit. The fetch_fqs.sh script contains commands for all samples in this study. 

fastq-dump SRR1660308

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

for b in *_Aligned.toTranscriptome.out.bam; do name=$(basename ${b} _Aligned.toTranscriptome.out.bam)
rsem-calculate-expression --alignments --paired-end -p 24 --no-bam-output --append-names ${b} /data2/lackey_lab/austin/dbs/rsem/ ${name} ;
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
2. a df for meta data and setting conditions for analysis, sample names must match that of the SAMPLE.isoform.results file. An example of this file is included as samp_meta.txt

| Sample | Condition |
|-----:|---------------|
|     Sample1|Mutant|
|     Sample2|Mutant|
|     Sample3|WT|
|     Sample4|WT|
