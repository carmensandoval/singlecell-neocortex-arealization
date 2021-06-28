# From abhaduri homefiles on 2020-11-23.
# General neocortex sample processing.

library(Seurat)
library(dplyr)
library(Matrix)
library(RANN)
library(igraph)
library(matrixStats)
library(qlcMatrix)
source("doFastPCA.R")

load(file="/kriegsteinlab/data1/aparna/homefiles/cleanedobjects_bysample_June2019/neocortex.RData")
Neocortex <- NormalizeData(object = Neocortex, normalization.method = "LogNormalize")
Neocortex <- FindVariableGenes(object = Neocortex, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
ribo <- read.table("/kriegsteinlab/data1/aparna/homefiles/ribogenes.txt", sep=" ", row.names=1)
ribo.genes <- intersect(rownames(Neocortex@data),rownames(ribo))
percent.ribo <- colSums(expm1(Neocortex@data[ribo.genes,]))/colSums(expm1(Neocortex@data))
Neocortex <- AddMetaData(Neocortex, percent.ribo, "percent.ribo")
Neocortex <- ScaleData(object = Neocortex, genes.use = Neocortex@var.genes, display.progress = FALSE, do.center=TRUE, vars.to.regress="Individual")
matrix_scaled <- Neocortex@scale.data
matrix_pca <- doFastPCA(t(matrix_scaled),50,center.use = F, scale.use = F,iterative = F)
ev <- matrix_pca$sdev^2
max(which(matrix_pca$sdev^2>(sqrt(length(row.names(Neocortex@data))/length(colnames(Neocortex@data))) + 1)^2))
u <- max(which(matrix_pca$sdev^2>(sqrt(length(row.names(Neocortex@data))/length(colnames(Neocortex@data))) + 1)^2))
sig_PCAs <- matrix_pca$rotation[,1:u]
cells_projected_sig_PCAs <- t(matrix_scaled) %*% sig_PCAs
nearest <- nn2(cells_projected_sig_PCAs, k=10)
rownames(nearest$nn.idx) <- rownames(cells_projected_sig_PCAs)
write.table(nearest$nn.idx, file="nn2_output_neighbors_sigPCAs_Neocortex.txt", sep="\t", quote=F, col.names=NA)
write.table(nearest$nn.dists, file="nn2_output_distance_sigPCAs_Neocortex.txt", sep="\t", quote=F, col.names=NA)
system("perl jaccard_toedges2.pl nn2_output_neighbors_sigPCAs_Neocortex.txt nn2_output_distance_sigPCAs_Neocortex.txt > jaccard_weighted_edgesNeocortex.txt")
edgedata <- read.table("jaccard_weighted_edgesNeocortex.txt", sep="\t", header=T)
edges <- as.data.frame(edgedata)
input <- graph_from_data_frame(edges, directed=FALSE)
clusters <- cluster_louvain(input)
clustered <- membership(clusters)
write(clustered, file="names_incol1Neocortex.txt", ncol=1)
names <- edgedata[,1]
uniquenames <- unique(names)
meta <- read.table("names_incol1Neocortex.txt", sep="\t")
rownames(meta) <- uniquenames
sizes(clusters)
Neocortex <- SetDimReduction(Neocortex, reduction.type = "fastPCA", slot = "cell.embeddings", new.data = cells_projected_sig_PCAs)
Neocortex <- SetDimReduction(Neocortex, reduction.type = "fastPCA", slot = "gene.loadings", new.data = matrix_pca$rotation)
Neocortex <- SetDimReduction(Neocortex, reduction.type = "fastPCA", slot = "sdev", new.data = matrix_pca$sdev)
Neocortex <- RunTSNE(Neocortex, reduction.use = "fastPCA", perplexity=30, dims.use = 1:ncol(cells_projected_sig_PCAs), do.fast = T, check_duplicates=FALSE)
Neocortex <- AddMetaData(Neocortex, meta, "groups")
Neocortex <- SetAllIdent(Neocortex, id="V1")

pdf("/kriegsteinlab/data1/aparna/homefiles/Neocortex_tsne.pdf")
TSNEPlot(Neocortex, do.label=TRUE)
dev.off()

write.table(Neocortex@ident, "/kriegsteinlab/data1/aparna/homefiles/Neocortex_clusteridentity.txt", 
            sep="\t", quote=F, col.names=NA)

save(Neocortex, cells_projected_sig_PCAs, matrix_pca, ev, matrix_scaled, 
     file="/kriegsteinlab/data1/aparna/homefiles/Neocortex.RData")

pdf("/kriegsteinlab/data1/aparna/homefiles/Neocortex_featureplots.pdf")

FeaturePlot(Neocortex, c("HOPX", "VIM","GLI3","BCL11B","SATB2","EOMES", "PPP1R17", "DLX6-AS1", "PDGFRA"),
            cols.use = c("grey","blue"), pt.size=0.1)
FeaturePlot(Neocortex, c("CRYM", "MKI67","APOE","CRYAB","FN1","CCL3", "HES6", "STMN2", "MEF2C"),
            cols.use = c("grey","blue"), pt.size=0.1)
dev.off()

cluster.markers <- FindAllMarkers(Neocortex, min.pct = 0.25, thresh.use = 0.25, 
                                  only.pos=TRUE, max.cells.per.ident=2000)
write.table(cluster.markers, "/kriegsteinlab/data1/aparna/homefiles/Neocortex_clustermarkers.txt", 
            sep="\t", quote=F, col.names=T)
