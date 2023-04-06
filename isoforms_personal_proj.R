library(ggplot2)
library(DESeq2)
library(tximportData)
library(tximport)
library(stringr)
library(magrittr)

#set dir
data <- "/path/to/isoform_result_files/"
#assign list of files
fls <- dir(data, "*isoforms.results$", full.names = TRUE)
#import sample meta data
samp_meta <- read.table(file = "sample_meta.txt",
                        header = T, sep = "\t")
samp_meta$Sample <- sub('genes', 'isoforms', samp_meta$Sample)
samp_meta$mutation[samp_meta$sf3b1_status == "wildtype" | samp_meta$sf3b1_status == "healthy"] <- "WT"
samp_meta$mutation[samp_meta$sf3b1_status == "mutated"] <- "Mutant"
rownames(samp_meta) <- samp_meta$Sample
#assigns sample names to data
names(fls) <- samp_meta$Sample
#import test sample to fetch transcript to gene IDs
test <- read.table(file = "SRR1660308.isoforms.results",
                   sep = "\t", header = TRUE)
#take the first two columns transcript_ID and gene_ID into a new df
tx2gene <- test[,1:2]

#use tximport to fetch counts
isos <- tximport(files = fls, type = "rsem",
                 txIn = TRUE, txOut = FALSE,
                 tx2gene = tx2gene)
#convert txi object to deseq object
ddstxi <- DESeqDataSetFromTximport(isos,
                                   colData = samp_meta,
                                   design = ~ mutation)


#apply prefiltering
keep <- rowSums(counts(ddstxi) >= 5) >= 6
ddstxi <- ddstxi[keep,]
ddstxi$mutation <- relevel(ddstxi$mutation, ref = "WT")

des_seq <- DESeq(ddstxi)

res <- results(des_seq)

plotMA(res, ylim=c(-10, 10))

write.table(as.data.frame(res), file = "/Users/herber4/Desktop/Med_Bioinfo/personal_project/des_output.txt", sep = "\t")

