library(scuttle)
library(scran)
library(irlba)
library(Rtsne)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(viridisLite)
library(viridis)
library(scDblFinder)

library(umap)
library(reticulate)
#use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

path2data   <- '/data2/hanna/synaptogenesis/newvolume/analysis/combined_h/all-well/DGE_unfiltered'
sample_info <- read.table('/data2/ivanir/Feline2023/ParseBS/newvolume/analysis/sample_info.tab',
  sep = "\t", header = TRUE)

#read the spare matrix into counts
#read the geneIDs, names and genome-of-origine into genes 
counts    <- t(readMM(paste0(path2data, "/DGE.mtx")))
genes     <- read.csv(paste0(path2data, "/all_genes.csv"))
metadata  <- read.csv(paste0(path2data, "/cell_metadata.csv"))

lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

genes_human <- genes[genes$genome == "hg38",]
dim(counts)
[1]   62704 2214461
dim(counts[,ngenes > 400 & lib.sizes > 500])
[1] 62704 56066

counts   <- counts[,ngenes > 400 & lib.sizes > 500]
metadata <- metadata[ngenes > 400 & lib.sizes > 500,]
lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

hist(ngenes/lib.sizes)

counts   <- counts[,ngenes/lib.sizes < 0.9]
metadata <- metadata[ngenes/lib.sizes < 0.9,]
lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

sample_bc1_well <- rep(NA, nrow(metadata))        
sample_number   <- rep(NA, nrow(metadata))
sample_name    <- rep(NA, nrow(metadata))

#editing sample info
sample_info$H_day <- sample_info$H_Timepoint 
sample_info$H_day <- gsub("55\\+","",sample_info$H_day)
sample_info$H_day <- as.integer(sample_info$H_day)
sample_info$H_day <-  sample_info$H_day +55
sample_info$Sample_name_H <- paste(sample_info$H_Batch, sample_info$H_day, sample_info$H_Replicate, sep="_")

sample_info$M_day <- sample_info$M_Timepoint
sample_info$M_day <- gsub("8\\+","",sample_info$M_day)
sample_info$M_day <- as.integer(sample_info$M_day)
sample_info$M_day <-  sample_info$M_day +8
sample_info$Sample_name_M <- paste(sample_info$M_Batch, sample_info$M_day, sample_info$M_Replicate, sep="_")

samples <- unique(sample_info$Sample_well)
for (i in 1:length(samples)){
  sample_bc1_well[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))] <- sample_info$Sample_well[i]
  sample_number[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]   <- sample_info$Sample_Number[i]
  sample_name[metadata$bc1_well %in% unlist(strsplit(samples[i],split=","))]     <- sample_info$Sample_name_H[i]
}

submeta <- data.frame(rlist::list.rbind(strsplit(sample_name, split="_")))
colnames(submeta) <- c("batch", "day", "replicate")
submeta$day <- gsub("d","",submeta$day)

metadata <- data.frame(cbind(metadata, lib.sizes, sample_number, sample_bc1_well, sample_name, submeta))
plot_df <- metadata
setwd('/data2/hanna/synaptogenesis/newvolume/analysis/QC_h')

ggplot(plot_df, aes (x = factor(sample_name), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_beforeQC.pdf")

pdf("cell_complexity.pdf")
qplot(lib.sizes, ngenes, col = ifelse(ngenes > 400 & lib.sizes > 500 , "drop", "keep")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  labs(x = "UMI count", y = "Number of expressed genes") +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()

ensembl <- useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl",mirror = "useast")
gene_map  <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol", values = genes$gene_name, mart = ensembl)
  
mt.index    <- gene_map$chromosome_name == "MT"
mt.counts   <- counts[which(genes$gene_name %in% gene_map$hgnc_symbol[mt.index]), ]
mt.fraction <- colSums(mt.counts/lib.sizes)
dim(mt.counts)
[1]    37 56018

mt.p  <- pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim<- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.05)])
mt.lim
[1] 0.05610907
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.001)])
mt.lim
[1] 0.07825714

metadata <- data.frame(cbind(metadata,mt.fraction))

pdf("Mtreadfraction.pdf")
qplot(lib.sizes, mt.fraction, col = ifelse(mt.fraction>mt.lim, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()
print(plot)

dim(counts[,mt.fraction < mt.lim])
[1] 62704 47004
dim(counts[,mt.fraction < 0.2])
[1] 62704 51985

mtlim <- 0.2

sce <- SingleCellExperiment(list(counts=counts[,mt.fraction < mt.lim]),
  colData=DataFrame(metadata[mt.fraction < mt.lim,]))
rownames(sce) <- genes$gene_id

rownames(genes) <- rownames(sce)
rowData(sce) <- DataFrame(genes)

colnames(sce) <- metadata$bc_wells[mt.fraction  < mt.lim]
colData(sce)  <- DataFrame(metadata[mt.fraction < mt.lim,])

lib.sizes <- colSums(counts(sce))
sce_filt  <- sce[calculateAverage(sce)>0.05,]

#class: SingleCellExperiment 
#dim: 6 47004 
#metadata(0):
#assays(1): counts
#rownames(6): ENSG00000000003 ENSG00000000005 ... ENSG00000000460
  ENSG00000000938
#rowData names(3): gene_id gene_name genome
#colnames: NULL
#colData names(20): bc_wells sample ... replicate mt.fraction
#reducedDimNames(0):
#mainExpName: NULL
#altExpNames(0):

#Cluster similar cells based on their expression profiles, using either log-expression values or ranks.
clusts <- as.numeric(quickCluster(sce_filt, method = "igraph", min.size = 100))
min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce_filtn <- computeSumFactors(sce_filt, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

#A size factor is a scaling factor used to divide the raw counts of a particular cell to obtain normalized expression values/
#Gets or sets the size factors for all cells in a SingleCellExperiment object.
sizeFactors(sce) <- sizeFactors(sce_filtn)

pdf("sizefactors.pdf")
ggplot(data = data.frame(X = lib.sizes, Y = sizeFactors(sce)), mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(500, 2000, 5000, 10000, 30000), labels = c("5,00", "2,000", "5,000", "10,000", "30,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  theme_minimal() +
  theme(text = element_text(size=20))  +
  labs(x = "Number of UMIs", y = "Size Factor")
dev.off()

ggplot(data.frame(colData(sce)), aes (x = factor(sample_name), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() +  coord_flip() +
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(100, 1000, 5000, 10000, 50000, 100000),
    labels = c("100","1,000", "5,000", "10,000", "50,000", "100,000"))
ggsave("UMIsBySample_afterQC.pdf")


library(BiocParallel)
bp <- MulticoreParam(12, RNGseed=1234)
bpstart(bp)

#detection and evaluation of doublets/multiplets
sce <- scDblFinder(sce, samples="bc1_well", dbr=.03, dims=30, BPPARAM=bp)
bpstop(bp)
table(sce$scDblFinder.class)
#singlet     doublet 
#43663        3341 
  
#normalisation
sce_filtn <- sce[calculateAverage(sce)>0.05,]
sce_filtn <- logNormCounts(sce_filtn)

decomp <- modelGeneVar(sce_filtn)
hvgs   <- rownames(decomp)[decomp$FDR < 0.5]
pca    <- prcomp_irlba(t(logcounts(sce_filtn[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_filtn)
tsne <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)



library(umap)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

df_plot <- data.frame(
 colData(sce),
 doublet  = colData(sce)$scDblFinder.class,
 tSNE1    = tsne$Y[, 1],
 tSNE2    = tsne$Y[, 2], 
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

plot.index <- order(df_plot$doublet)
ggplot(df_plot[plot.index,], aes(x = tSNE1, y = tSNE2, col = factor(doublet))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=c("gray","#0169c1"), name = "") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("tsne_doublets.pdf")

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(doublet))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=c("gray","#0169c1"), name = "") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))

colData(sce) <- DataFrame(df_plot)

saveRDS(sce,paste0(path2data,"sce.rds"))
