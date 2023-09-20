#### 1. Analysis of single nucleus RNA-seq in adult heart sample(GSE128628)###
#### 2. Analysis of single nucleus RNA-seq in P2, P4, P9, P11 Heart Samples (GSE130699)###
#### The steps involve Quality control, Data Normalization, Scaling, Data Integration, and Cell Clustering.
#### We filtered out non cardiomyocytes and preserved the cardiomyocytes single nucleus data matrix.
#### 3. Analysis of single nucleus RNA-seq of cardiomyocytes (GSE130699 and GSE128628) ###
#### Loading Single nucleus RNA-seq data of cardiomyocytes at different stages ####
library(devtools)
library(Seurat)
library(cowplot)
library(harmony)
library(tidyverse)
library(patchwork)
library(dplyr) 
library(scRNAtoolVis)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(scales)
#### laod data which was integrated and clustered ####
snRNA <- readRDS("sn_RNA_CM_ALL_clustered.rds")

#### Normalization > FindVariableFeatures > scale > PCA ####
# snRNA <- NormalizeData(snRNA) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
#### integration and remove batch effects ####
# system.time({snRNA <- RunHarmony(snRNA, group.by.vars = "orig.ident")})
#### FindNeighbors ####
# snRNA <- FindNeighbors(snRNA, dims = 1:50)
# snRNA <- FindClusters(snRNA, resolution = 0.5)
# table(snRNA$seurat_clusters)

#### RunUMAP ####
# snRNA <- RunUMAP(snRNA, reduction = "harmony", dims = 1:30)
p <- DimPlot(snRNA, reduction = "umap", label = T)
ggsave(file.path( "1-UMAP-Cluster.pdf"),width =7, height =5,p)
#group_by_sample
p = DimPlot(snRNA, reduction = "umap", group.by = 'orig.ident') 
ggsave(file.path(Redout, "1-UMAP-Group.pdf"),width =7, height =5,p)

##### Quality control ####
# DefaultAssay(snRNA) <- "RNA"
# snRNA[["percent.mt"]] <- PercentageFeatureSet(snRNA, pattern = "^MT-|^mt-")
# HB.genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz"
#               ,"HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
# HB_m <- match(HB.genes, rownames(snRNA@assays$RNA))
# HB.genes <- rownames(snRNA@assays$RNA)[HB_m]
# HB.genes <- HB.genes[!is.na(HB.genes)]
# snRNA[["percent.HB"]]<-PercentageFeatureSet(snRNA, features=HB.genes)
# 
# p <-VlnPlot(snRNA,group.by = "orig.ident",
#             features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
#             pt.size = 0,
#             ncol = 4) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave(file.path(qcout, "1-vlnplot.QC.pdf"),width =7, height =5,p)

##### Tnnt2 expression #####
P <- FeaturePlot(snRNA,features = "Tnnt2",reduction = "umap")
ggsave(file.path(Redout, "TNNT2_umap.pdf"),width =7, height =5,P)
P <- VlnPlot(snRNA,features = "Tnnt2",pt.size = 0)
ggsave(file.path(Redout, "TNNT2_VN.pdf"),width =7, height =5,P)
# save reduction data
saveRDS(snRNA,"snRNA_CM_reduction.rds")

##### find markers of each cluster #####
snRNA.markers <- FindAllMarkers(snRNA,logfc.threshold= 0.25 ,test.use= "wilcox" ,
                                    min.pct=0.1,only.pos=TRUE)
snRNA.markers <- snRNA.markers %>% group_by(cluster)
markegene <- as.data.frame(snRNA.markers)
head(snRNA.markers)
write.csv(markegene,file="markerGenesCLU.csv",quote = F)

# ##### Heatmap show markers of each cluster #####
# top4 <- snRNA.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
# TOP4 <- top4$gene

##### selected marker genes according to related literatures #####
marker_CM_sle <- c("Ryr2","Myh6","Myh7","Tnni1","Actc1","Dab2","Runx1","Mki67","Cenpp",
               "Cav1","Pecam1","Dcn","Lum","Xirp2","Nppa","Nppb")
##### set color #####
identcolors = hue_pal()(length(levels(snRNA)))
P <- AverageHeatmap(object = snRNA,annoCol = T,myanCol = identcolors,
                    markerGene = marker_CM_sle)
ggsave(file.path(Redout, "maker_clu.pdf"),width =7, height =5,P)
##### cluster correlation analysis #####
av <-AverageExpression(snRNA,
                       group.by = "seurat_clusters",
                       assays = "RNA")
av=av[[1]]
#Select the 2000 genes with the highest standard deviation
cg=names(tail(sort(apply(av, 1, sd)),2000))
cor_clu<- cor(av[cg,],method = 'spearman')
write.csv(cor_clu,"cor_clu.csv")
#pheatmap view
P <- pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),
                        color=colorRampPalette (c ("#08519C","#F9FBFC","#DE2D26")) (100)) 
ggsave(file.path(Redout, "CLU_COR.pdf"),width =7, height =5,P)

##### manually annotated 6 classes of cardiomyocyte subtypes according to selected marker genes
# and based on selected marker genes and the correlation between each of these clusters #####
new.cluster.ids <- c("CM1","CM2","CM1","CM2","CM1","CM3",
                     "CM1","CM4","CM5","CM1","CM3","CM3",
                     "CM1","CM1","CM6")
snRNA_CELL <- snRNA
names(new.cluster.ids) <- levels(snRNA_CELL)
snRNA_CELL <- RenameIdents(snRNA_CELL, new.cluster.ids)
# # save data 
# saveRDS(snRNA_CELL,"snRNA_CELL.rds")
# UMAP visulization
P <- DimPlot(snRNA_CELL, reduction = 'umap', 
             label = TRUE)
ggsave(file.path(CELL1out, "CM_CELL_UMAP.pdf"),width =7, height =5,P)

# find markers of each CM cell
snRNA_CELL.markers <- FindAllMarkers(snRNA_CELL,logfc.threshold= 0.25 ,test.use= "wilcox" ,
                                     min.pct=0.1,only.pos=TRUE)
snRNA_CELL.markers <- snRNA_CELL.markers %>% group_by(cluster)
markegene <- as.data.frame(snRNA_CELL.markers)
head(snRNA_CELL.markers)
write.csv(markegene,file="markerGenes_CM_CELL.csv",quote = F)

##### Heatmap show markers of each cluster #####
# top6 <- markegene %>% group_by(cluster) %>% top_n(6, avg_log2FC)
# TOP6 <- top6$gene
# identcolors = hue_pal()(length(levels(snRNA_CELL)))
# P <- AverageHeatmap(object = snRNA_CELL,annoCol = T,myanCol = identcolors,
#                     markerGene = TOP6)

#### pheatmap view of selected marker genes of each CMs cluster #####
mark_CM_SLE<- read.csv("Selected_Marker_Genes_CM_Cluster.CSV",header = T,sep=",")
mark_CM_SLE <- mark_CM_SLE$GENE
identcolors = hue_pal()(length(levels(snRNA_CELL)))
P <- AverageHeatmap(object = snRNA_CELL,annoCol = T,myanCol = identcolors,
                    markerGene = mark_CM_SLE)
ggsave(file.path(CELL1out, "maker_CELL.pdf"),width =7, height =5,P)

#### CELLPORT anlysis ####
library(stringr)
snRNA_CELL1 <- snRNA_CELL
snRNA_CELL1$celltype.P <- paste(Idents(snRNA_CELL1), snRNA_CELL1$orig.ident, sep = "_") 

CELL_PART<- as.data.frame(table(snRNA_CELL1$celltype.P))
CELL_PART_CELL <- unlist(lapply(X = as.character(CELL_PART$Var1), FUN = function(x) 
{return(strsplit(x, split = "_")[[1]][1])}))
CELL_PART_GROUP <- unlist(lapply(X = as.character(CELL_PART$Var1), FUN = function(x) 
{return(strsplit(x, split = "_")[[1]][2])}))

CELL_PART$GROUP <- CELL_PART_GROUP
CELL_PART$CELL <- CELL_PART_CELL
head(CELL_PART)
CELL_PART <- CELL_PART[,-1]
# calculate percent
CP_sum <- aggregate(CELL_PART$Freq,by=list(CELL_PART$GROUP),sum)
colnames(CP_sum) <- c("GROUP","ALL")
CP_FIN <- merge(CELL_PART,CP_sum,by.x="GROUP",by.y="GROUP")
head(CP_FIN)
CP_FIN$per <- CP_FIN$Freq/CP_FIN$ALL
CP_FIN$per <- CP_FIN$per*100
head(CP_FIN)
write.csv(CP_FIN,"CELL_PORT_PER.csv")
CELLPORT <- CP_FIN[,c(1,3,5)]
head(CELLPORT)
colnames(CELLPORT)=c("Group","cluster","proportion")
# factor order

CELLPORT$Group <- ordered(CELLPORT$Group, 
                          levels = c("P2","P4","P9","P11","Adult"))
CELLPORT$cluster <- ordered(CELLPORT$cluster, 
                            levels = names(table(Idents(snRNA_CELL))))
# cellport view
p <- ggplot(CELLPORT,aes(Group,proportion,fill=cluster))+
  geom_bar(stat="identity",width = 0.6,position="fill")+
  ggtitle("")+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme_classic()
# set color
identcolors = hue_pal()(length(levels(snRNA_CELL)))
p <- p+scale_fill_manual(values = identcolors)
ggsave(file.path(CELL1out, "CELL_PORT.pdf"),width =7, height =5,p)

#### monocle analysis ####
library(monocle)
library(getopt)
library(Seurat)
library(clusterProfiler)
source("order_cells.R")
dir_file <- getwd()
MONO = file.path(dir_file, "MONO")
if (!dir.exists(MONO)) dir.create(MONO)
# load data
snRNA_CELL$cell_type <- Idents(snRNA_CELL)
# parameter set
options(stringsAsFactors=F)
str(snRNA_CELL)
# Extract expression matrixs, phenotypes, and genes from Seurat objects
expr_matrix <- as(as.matrix(snRNA_CELL@assays$RNA@counts), 'sparseMatrix')
print(head(expr_matrix[,1:4]))
p_data <- snRNA_CELL@meta.data
p_data$celltype <- snRNA_CELL@active.ident
head(p_data)
f_data <- data.frame(gene_short_name = row.names(snRNA_CELL), row.names = row.names(snRNA_CELL))
head(f_data)
## create monocle2 object
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
# estimate  SizeFactors and Dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
# sort genes
cds <- detectGenes(cds , min_expr =0.1)
head(fData(cds))
## remove genes expressed less 5 cells
express_genes <- row.names(subset(fData(cds),num_cells_expressed>=5))
## DifferentialGeneTest function
diff <- differentialGeneTest(cds[express_genes,],fullModelFormulaStr="~celltype",cores=16)
## remove genes with qvalue >= 0.01; sort genes
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
# write.table(deg,file="snRNA_CELL.monocle.DEG.xls",col.names=T,row.names=F,sep="t",quote=F)
write.csv(deg,file.path(MONO,"snRNA_CELL.monocle.DEG.csv"))
# select ordered genes 
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
# reduction “DDRTree”
cds <- reduceDimension(cds,max_components = 2, method = DDRTree)
# Trajectory construction
cds <- orderCells(cds)
# definite root state
cds <- orderCells(cds,root_state =3)
## Visualization of trajectories
types <- c("Pseudotime","State","celltype","orig.ident")
for (type in types) {
  plot_cell_traj <- plot_cell_trajectory(cds,color_by=type, cell_size=1,show_backbone=TRUE)
  pdf(file.path(MONO,paste("monocle",type,"pdf",sep=".")))
  print(plot_cell_traj)
  dev.off()
}
#### define state ####
State_info <- data.frame(ID=colnames(cds),state=cds$State)
head(State_info)
#### corrected ####
state_cor <- c()
for (state in State_info$state) {
  if (state == "1") {
    state_cor <- c(state_cor, 'state2')
  } else if (state == "2") {
    state_cor <- c(state_cor, 'state3')
  } else {
    state_cor <- c(state_cor, 'state1')
  }
}
State_info$state_cor <- state_cor
head(State_info)  
#### Add state information to snRNA_CELL ####
snRNA_CELL$state <- State_info$state_cor
# umap visulization
DimPlot(snRNA_CELL,reduction = "umap",group.by = "state")

##### find markers of each state #####
snRNA_state = SetIdent(snRNA_CELL, value = snRNA_CELL$state)
state.markers <- FindAllMarkers(snRNA_state,threshold= 0.25 ,test.use= "wilcox" ,
                                min.pct=0.05,only.pos=TRUE)
state.markers <- state.markers %>% group_by(cluster)
markegene_state <- as.data.frame(state.markers)
head(markegene_state)
markegene_state <- filter(markegene_state,p_val_adj<0.05)
write.csv(markegene_state,file="markegene_state.csv",quote = F)

## selected marker genes of each state to heatmap visualization 
top30 <- markegene_state %>% group_by(cluster) %>% top_n(30, avg_log2FC)
top30_mark <- top30$gene
br_heat <- plot_genes_branched_heatmap(cds[top30_mark,],
                                       branch_point = 1,
                                       num_clusters = 3,
                                       cores = 16,
                                       branch_labels = c("Cell fate1","Cell fate2"),
                                       # hmcols = colorRampPalette()
                                       hmcols = NULL,
                                       branch_colors =c("#979797","#F05662","#7990C8"),
                                       use_gene_short_name = T,
                                       show_rownames = F,
                                       return_heatmap = T)
##### GO analysis of marker genes of each state #####
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
# Turn symbol to Entrez ID
ids=bitr(markegene_state$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
# merge symbol and Entrez ID
sce.markers=merge(markegene_state,ids,by.x='gene',by.y='SYMBOL')
# split bt cluster(state)
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 
## GO analysis
go_state <- compareCluster(gcSample,
                           fun = "enrichGO",
                           OrgDb = "org.Mm.eg.db",#org.Hs.eg.db
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.01,
                           qvalueCutoff = 0.05)
# 将EntrezID转换为Symbol
go_state = setReadable(go_state,
                       OrgDb = "org.Mm.eg.db",
                       keyType = "ENTREZID") 
write.csv(go_state,file="go_state.csv",quote=F,row.names = F)

##### obatain state1 high expression RBPs #####
state1_mark <- subset(markegene_state,cluster=="state1")
sate1_RBPs <- intersect(state1_mark$gene,RBP2GO$GENE)
# GO analysis
gene_info <- bitr(as.character(sate1_RBPs), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
target_gene_id <- as.character(gene_info[, 2])

GO_BP <- enrichGO(OrgDb="org.Mm.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable = TRUE)
GO_BP <- as.data.frame(GO_BP)
write.table(GO_BP, "GO_BP_State1_RBPs.xls", quote = F, sep = "\t", col.names = T, row.names = F)

##### Gene trend analysis of state1 high expression RBPs #####
library(Mfuzz)
# get expression matrix
RBP_matrix <- AverageExpression(snRNA_CELL,features = sate1_RBPs,
                                slot = 'data' ,verbose = FALSE,group.by = "orig.ident")$RNA

#create Mfuzz analysis object
mfuzz_class <- new('ExpressionSet',exprs = RBP_matrix)
#Processing loss or abnormal value
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
#normalization
mfuzz_class <- standardise(mfuzz_class)
# divide into 4 clusters
set.seed(123)
cluster_num <- 4
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
#Visulization
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 2), time.labels = colnames(RBP_matrix))
#merge
RBP_cluster <- mfuzz_cluster$cluster
RBP_cluster <- cbind(RBP_matrix[names(RBP_cluster ), ], RBP_cluster )
RBP_cluster <- as.data.frame(RBP_cluster )
RBP_cluster$GENE <- rownames(RBP_cluster )
write.csv(RBP_cluster , 'RBP_cluster.csv', col.names = T, quote = FALSE)

##### The interaction network of proteins was analyzed using the STRING #####
# Networks with a combined score exceeding 0.8 were retained for further analysis.
##### Centrality measures of network were evaluated by igragh package. #####
library(igraph)
edges <- read.csv("network_filterd.CSV",header = T,sep=",",check.names = F)
colnames(edges) <- c("from","to")
g <- make_graph(t(edges),directed = FALSE)
# caculate degree
g_degree <- degree(g,v=V(g),mode = "all",loops = TRUE)
g_degree_name<- names(g_degree)
g_degree_value <- as.numeric(g_degree)
g_degree_final <- data.frame(GENE=g_degree_name,degree=g_degree_value)
# caculate closeness
g_close <- closeness(g,v=V(g),mode = "all",normalized = T)
g_close_name<- names(g_close)
g_close_value <- as.numeric(g_close)
g_close_final <- data.frame(GENE=g_close_name,close=g_close_value)
# caculate betweeness
g_bet <- betweenness(g,v=V(g),directed = TRUE,weights = NULL,
                     nobigint = TRUE,
                     normalized = FALSE)
g_bet_name <- names(g_bet)
g_bet_value <- as.numeric(g_bet)
g_bet_final <- data.frame(GENE=g_bet_name,between=g_bet_value)
# caculate closeness
g_eig <- eigen_centrality(g)
# caculate Eigenvector
g_eig_name <- names(g_eig$vector)
g_eig_value <- as.numeric(g_eig$vector)
g_eig_final <- data.frame(GENE=g_eig_name,eigen_centrality=g_eig_value)
# caculate PageRank
g_page <- page.rank(g, algo = "prpack", 
                    vids = V(g), directed = TRUE, damping = 0.85, 
                    personalized = NULL, weights = NULL, options = NULL)
g_page_name <- names(g_page$vector)
g_page_value <- as.numeric(g_page$vector)
g_page_final <- data.frame(GENE=g_page_name,page=g_page_value)
# merge
data_list <- list(g_degree_final,g_close_final,g_bet_final,g_eig_final,g_page_final)
my_merge <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "GENE")
}
g_arg_net <- Reduce(my_merge, data_list)
write.csv(g_arg_net,"g_arg_net.csv")