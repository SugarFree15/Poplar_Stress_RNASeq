# Poplar_Stress_RNASeq
Source code for M.S. report involving RNA-Seq of poplar stress.


## Downloading SRA Files From NCBI

fasterq-dump ERR1864411

> For this I have demonstrated with one of the 81 samples used (ERR1864411). To use this package, you must have the sra-toolkit from NCBI installed along your PATH, otherwise provide the path to the program as demonstrated below. This command will download the SRA RNA-Seq data of the accession number provided to the current working directory in fastq format. In this case, the samples were paired-end reads, resulting in downloading both ERR1864411_1.fastq and ERR1864411_2.fastq into the working directory.

~/Desktop/Software/sra-toolkit/bin/fasterq-dump ERR1864411


## Quality Control of Raw Reads

> To ensure high confidence in the downloaded sequencing data, we first trimmed weak reads and any remaining Illumina adapters with the following command. For this, cutadapt must be installed along your PATH. Again, this was run with all samples used, but demonstrated with only one.

cutadapt -q 20,20 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 50 --max-n 0 -o ERR1864411_1.cutadapt.fastq -p ERR1864411_2.cutadapt.fastq ERR1864411_1.fastq ERR1864411_2.fastq

> From this, two new files are outputted containing the new, trimmed reads named ERR1864411_1.cutadapt.fastq and ERR1864411_2.cutadapt.fastq. To verify these new file hold high confidence, we assess the quality using FastQC, which will need to be installed along your PATH.

fastqc ERR1864411_1.cutadapt.fastq
fastqc EER1864411_2.cudadapt.fastq


## Developing Poplar Genome Reference

> For the eventual alignment and quantification of transcripts, references are needed. For the command below, you will need RSEM and STAR installed to your PATH and the Phytozome P. trichocarpa v4.1 genome and accompanying grr3 annotation file downloaded into the working directory.

rsem-prepare-reference --star --gff3 Ptrichocarpa_533_v4.1.gene.gff3 Ptrichocarpa_533_v4.0.fa Potri_4.0


## Developing de novo Transcriptome Reference

> For the unannotated pri-miRNA genes, we required a de novo transcriptome as a reference. Using the control samples, we developed the transcriptome with rnaSPAdes installed to our PATH.

rnaspades.py ERR1864411_1.cutadapt.fastq ERR1864411_2.cutadapt.fastq ERR1864412_1.cutadapt.fastq ERR1864412_2.cutadapt.fastq ERR1864413_1.cutadapt.fastq ERR1864413_2.cutadapt.fastq

> From there, we developed the reference similar the the section above.

rsem-prepare-reference --star hard_filtered_transcripts.fasta transcriptome_de_novo


## Aligning and Quantifying Transcripts

> With the references from above, we then ran each sample through the following commands: one for use in studying GRFs and the other for pri-miR396.

rsem-calculate-expression --star --no-bam-output --paired-end ERR1864411_1.cutadapt.fastq ERR1864411_2.cutadapt.fastq Potri_4.0 ERR1864411_gen
rsem-calculate-expression --star --no-bam-output --paired-end ERR1864411_1.cutadapt.fastq ERR1864411_2.cutadapt.fastq transcriptome_de_novo ERR1864411_novo

> From these will be .results files with the names given at the end of the command. All of these outputs were separated into folders labelled GRF_data and miR_data within the working directory to be used for the R environment below.


## Differential Expression Analysis of GRFs in R

> Once all the output files from above were placed into their respective folders within the working directory for the R environment (~/Potri_stress_RNASeq in this case), we could begin importing this data into the R environment.

library(readr)
library(tidyverse)
library(DESeq2)
library(tximport)

> These prepared some packages that will be needed later on. 
> For the next part, a metadata file is needed (named metadata.tsv) containing a column labelled "id" filled with the name of every sample in the order they are present in the GRF_data and miR_data folders and a column labeled "dex" containing a 2-3 letter code for the type of sample (i.e. LSD = leaf, short-term, drought or LC = leaf control), which will be used to clump the replicates together.

setwd("~/Potri_stress_RNASeq")

GRFcommonNames = c("GRF01", "GRF02", "GRF03", "GRF04", "GRF05", "GRF06", "GRF07", "GRF08", "GRF09", "GRF10", "GRF11", "GRF12", "GRF13", "GRF14", "GRF15", "GRF16", "GRF17", "GRF18", "GRF19")

metadata = read_tsv("metadata.tsv")

sample_files = c(list.files(path="GRF_data"))

> After all the above commands prepare items for later steps, next is to import the RSEM data into R.

setwd("GRF_data")

GRF_data = tximport(files = sample_files,
                           type = "rsem",
                           geneIdCol = "gene_id",
                           countsCol = "expected_count",
                           abundanceCol = "FPKM",
                           ignoreTxVersion = TRUE)

setwd("~/Potri_stress_RNASeq")

> Next, we used the metadata to name the data columns and set up the transcript counts for the DESeq software. This ensures the program will run smoothly with no errors.

colnames(GRF_data$abundance) <- metadata$id
colnames(GRF_data$counts) <- metadata$id
colnames(GRF_data$length) <- metadata$id

GRF_counts = GRF_data$counts
GRF_counts = as.data.frame(GRF_counts)
GRF_counts = GRF_counts+1
GRF_counts=rownames_to_column(GRF_counts, var = "gene_id")
GRF_counts[-1] = lapply(GRF_counts[-1], as.integer)

> Here, the data is formatted for the calculations and then run through DESeq

dds_GRF.data = DESeqDataSetFromMatrix(countData = GRF_counts,
                                      colData = metadata,
                                      design = ~dex,
                                      tidy = TRUE)
dds_GRF = DESeq(dds_GRF.data)

> Now, there is a DESeq object containing all the data, which is then filtered out using the following commands, repeated for each condition tested. This example demonstrated the stem xylem, long-term, cold stress comparison to the control.

res_SLC <- results(dds_GRF, contrast = c("dex","SLC", "SC"), tidy = TRUE)
res_SLC <- tbl_df(res_SLC)

resSLC_GRFs <- filter(res_SLC, row=="Potri.001G082700.v4.1" | row=="Potri.001G114000.v4.1" | row=="Potri.001G132600.v4.1" | row=="Potri.001G169100.v4.1" | row=="Potri.002G115100.v4.1" | row=="Potri.003G065000.v4.1" | row=="Potri.003G100800.v4.1" | row=="Potri.003G118100.v4.1" | row=="Potri.006G115200.v4.1" | row=="Potri.006G143200.v4.1" | row=="Potri.007G007100.v4.1" | row=="Potri.012G022600.v4.1" | row=="Potri.013G077500.v4.1" | row=="Potri.014G007200.v4.1" | row=="Potri.014G012800.v4.1" | row=="Potri.014G071800.v4.1" | row=="Potri.015G006200.v4.1" | row=="Potri.018G065400.v4.1" | row=="Potri.019G042300.v4.1" )
resSLC_GRFs$row = GRFcommonNames
resSLC_GRFs <- mutate(resSLC_GRFs,Significance=ifelse(padj<=0.05,"Significant","Non-significant"))

> At this point, there is now two table, one of all genes present (res_SLC) and only GRFs (resSLC_GRFs) which we then plot the total transcriptome and the GRFs alone in volcano plots.

ggplot(res_SLC)+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), color=GRF))+
  ggtitle("Differential Expression of Poplar Stem Xylem Genes after Long-term Cold Stress")

ggplot(resSLC_GRFs)+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), color=Significance))+
  geom_text(check_overlap = TRUE, aes(x=log2FoldChange, y=-log10(padj), label=resSLC_GRFs$row), nudge_x = .15)+
  ggtitle("Long-term Effects of Cold Stress on GRF Expression in Poplar Stem Xylem")

> This process was then repeated through all the sample types, generating the GRF plot present in the report.


## Differential Expression of pri-miR396 in R

> This process is near identical to the section previous at the start, with the exception of file and folder titles. The same metadata file can be used for this process, also.

library(readr)
library(tidyverse)
library(DESeq2)
library(tximport)

miRcommonNames = c("miR396a", "miR396e", "miR396d", "miR396b")

setwd("~/Potri_stress_RNASeq")

metadata = read_tsv("metadata.tsv")

sample_files = c(list.files("miR_data"))

setwd("miR_data")

miR_data = tximport(files = sample_files,
                     type = "rsem",
                     geneIdCol = "gene_id",
                     countsCol = "expected_count",
                     abundanceCol = "FPKM",
                     ignoreTxVersion = TRUE)

setwd("~/Potri_stress_RNASeq")

colnames(miR_data$abundance) <- metadata$id
colnames(miR_data$counts) <- metadata$id
colnames(miR_data$length) <- metadata$id
miR_abundance = as.data.frame(miR_data$abundance)

miR_counts = miR_data$counts
miR_counts = as.data.frame(miR_counts)
miR_counts=rownames_to_column(miR_counts, var = "gene_id")
miR_counts[-1] = lapply(miR_counts[-1], as.integer)

dds_miR.data = DESeqDataSetFromMatrix(countData = miR_counts,
                                       colData = metadata,
                                       design = ~dex,
                                       tidy = TRUE)
dds_miR = DESeq(dds_miR.data)

> From here, I will again demonstrate one of the samples, though this process was repeated for every sample type.

res_SLCr <- results(dds_miR, contrast = c("dex","SLC", "stem_control"), tidy = TRUE)
res_SLCr <- tbl_df(res_SLCr)

res_SLCr <- filter(res_SLCr, row=="NODE_13714_length_1435_cov_125.716450_g9817_i0" | row=="NODE_6976_length_1993_cov_41.498971_g5023_i0" | row=="NODE_29758_length_645_cov_4.515101_g21784_i0" | row=="NODE_20572_length_1054_cov_5.062687_g14590_i0" )
res_SLCr$row = miRcommonNames
res_SLCr <- mutate(res_SLCr,Significance=ifelse(padj<=0.05,"Significant","Non-significant"))

ggplot(res_SLCr)+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), color=Significance))+
  geom_text(check_overlap = TRUE, aes(x=log2FoldChange, y=-log10(padj), label=row), nudge_x = .2)+
  ggtitle("Long-term Effects of Cold Stress on miR396 Expression in Poplar Stem Xylem")
  
