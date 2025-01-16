library(patchwork)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
options(warn=-1)
set.seed(1)
library(Hmisc)

GCB1.data<-Read10X("D:/GCB1_2.2.filtered_feature_bc_matrix/")
GCB2.data<-Read10X("D:/GCB2_2.2.filtered_feature_bc_matrix/")

GCB1.data<-as.data.frame(GCB1.data)
GCB2.data<-as.data.frame(GCB2.data)

# GCB1.metadata<-read.table("d:/GCB1metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# GCB2.metadata<-read.table("d:/GCB2metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# 
# GCB1.data<-GCB1.data[,rownames(GCB1.metadata)]
# GCB2.data<-GCB2.data[,rownames(GCB2.metadata)]

for (i in 1:length(colnames(GCB1.data))) {
  colnames(GCB1.data)[i] <- paste(colnames(GCB1.data)[i],"GCB1",i,sep = "-")  
}

for (i in 1:length(colnames(GCB2.data))) {
  colnames(GCB2.data)[i] <- paste(colnames(GCB2.data)[i],"GCB2",i,sep = "-")  
}

# rownames(GCB1.metadata)<-colnames(GCB1.data)
# rownames(GCB2.metadata)<-colnames(GCB2.data)
# 
# GCB1.metadata$group<-rep("GCB1",length(colnames(GCB1.data)))
# GCB2.metadata$group<-rep("GCB2",length(colnames(GCB2.data)))

# pbmc.metadata<-rbind(GCB1.metadata,GCB2.metadata)

pbmc.data<-cbind(GCB1.data,GCB2.data)
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB",min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F,vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F,npcs = 100)
ElbowPlot(pbmc,ndims = 100)

library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref.data<-read.table("D:/CJJ immunity/GSE109732/GSE109732_GeneExpression_DZ_Fraction1-3_GC-PB.txt",sep = "\t",header = TRUE)
ref.data<-ref.data[!duplicated(ref.data$Gene),]
rownames(ref.data)<-ref.data$Gene
ref.data<-ref.data[,-1]
ref.metadata<-data.frame(colnames(ref.data))
ref.metadata$Clustername<-ref.metadata[,1]
colnames(ref.metadata)[1]<-"ID"
rownames(ref.metadata)<-ref.metadata[,1]
ref <- CreateSeuratObject(counts = ref.data, project = "ref3k",meta.data = ref.metadata,min.cells = 3)

ds_seurat<-SummarizedExperiment(assays=list(counts=ref@assays$RNA@counts),colData = ref.metadata[,c("ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:14)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:14)

Idents(pbmc)<-pbmc$Cluster
DimPlot(pbmc,reduction = "ldaumap",label = TRUE)
DimPlot(pbmc,reduction = "ldatsne")

ref.data2<-ref.data[,1:12]
ref.metadata2<-ref.metadata[1:12,]
ds_seurat<-SummarizedExperiment(assays=list(counts=ref.data2),colData = ref.metadata2[,c("ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster2<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster2,reduction.name = "lda2")
pbmc<-RunUMAP(pbmc,reduction = "lda2",reduction.name = "lda2umap",dims = 1:11)
pbmc<-RunTSNE(pbmc,reduction = "lda2",reduction.name = "lda2tsne",dims = 1:11)

Idents(pbmc)<-pbmc$Cluster2
DimPlot(pbmc,reduction = "lda2umap",label = TRUE)
DimPlot(pbmc,reduction = "lda2tsne")

VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
Idents(pbmc)<-pbmc$Cluster
new.cluster.ids <- c("DZ_2", "Fraction3_2", "Fraction3_3",  "DZ_1","Fraction3_1", "DZ_3",
                     "Fraction2_3", "Fraction2_2", "prePB","Fraction2_1","prePB","prePB","prePB","prePB","prePB")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, label = TRUE, pt.size = 0.5) + NoLegend()
pbmc@meta.data$newcluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$newcluster,reduction.name = "lda3")
pbmc<-RunUMAP(pbmc,reduction = "lda3",reduction.name = "lda3umap",dims = 1:9)
pbmc<-RunTSNE(pbmc,reduction = "lda3",reduction.name = "lda3tsne",dims = 1:9)

Idents(pbmc)<-pbmc$newcluster
DimPlot(pbmc,reduction = "lda3umap")
DimPlot(pbmc,reduction = "lda3tsne")

Idents(pbmc)<-pbmc$Cluster2
new.cluster.ids <- c("DZ_2", "Fraction3_2", "Fraction3_1", "DZ_1", "DZ_3",
                     "Fraction2_3", "Fraction3_3","Fraction2_2", "prePB","Fraction2_1","prePB","prePB")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, label = TRUE, pt.size = 0.5) + NoLegend()
pbmc@meta.data$newcluster2<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$newcluster2,reduction.name = "lda4")
pbmc<-RunUMAP(pbmc,reduction = "lda4",reduction.name = "lda4umap",dims = 1:9)
pbmc<-RunTSNE(pbmc,reduction = "lda4",reduction.name = "lda4tsne",dims = 1:9)

Idents(pbmc)<-pbmc$newcluster2
DimPlot(pbmc,reduction = "lda4umap")
DimPlot(pbmc,reduction = "lda4tsne")
remove(ds_seurat,GCB1.data,GCB2.data,pbmc.data,pbmc.data_t,pred.hesc,ref.data,ref.data2,ref.metadata,ref.metadata2,all.genes,i,new.cluster.ids)
gc()

FeaturePlot(pbmc,features =c("Irf4"),order = TRUE,cols = c("lightgrey", "Blue"),reduction = "lda3umap")
FeaturePlot(pbmc,features =c("Kdm6b"),order = TRUE,cols = c("lightgrey", "Blue"),reduction = "lda3umap")
FeaturePlot(pbmc,features =c("nFeature_RNA"),order = TRUE,cols = c("lightgrey", "Blue"),reduction = "lda4umap",min.cutoff = "q25",max.cutoff = "q75")

FeaturePlot(pbmc,features =c("Ly75"),order = TRUE,cols = c("lightgrey", "Blue"),reduction = "lda4umap")
FeaturePlot(pbmc,features =c("plasmacell"),max.cutoff = "q75",reduction = "lda4umap",order = TRUE,min.cutoff = "q25") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Irf4"),cols = c("lightgrey", "Blue"),reduction = "lda4umap")
FeaturePlot(pbmc,features =c("Kdm6b"),cols = c("lightgrey", "Blue"),reduction = "lda4umap")
FeaturePlot(pbmc,features =c("Ly75"),cols = c("lightgrey", "Blue"),reduction = "lda4umap")
FeaturePlot(pbmc,features =c("plasmacell"),reduction = "lda4umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Glutamine1"),max.cutoff = "q75",reduction = "lda4umap",order = TRUE,min.cutoff = "q25") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Glutamine2"),max.cutoff = "q75",reduction = "lda4umap",order = TRUE,min.cutoff = "q25") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Glutamine3"),max.cutoff = "q75",reduction = "lda4umap",order = TRUE,min.cutoff = "q25") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Glutamine1"),reduction = "lda4umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Glutamine2"),reduction = "lda4umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Glutamine3"),reduction = "lda4umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Gls"),cols = c("lightgrey", "Blue"),reduction = "lda4umap")
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE,pt.size =0)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
VlnPlot(pbmc,features = c("Glutamine1","Glutamine2","Glutamine3","plasmacell"),sort = TRUE,pt.size =0)
DotPlot(pbmc,features =c("Glutamine1","Glutamine2","Glutamine3","plasmacell"))


Glutamine_1<-read.table("d:/GSE60927.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$plasmacell<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)

Glutamine_1<-read.table("d:/GOBP_GLUTAMINE_FAMILY_AMINO_ACID_BIOSYNTHETIC_PROCESS.v2023.2.Hs.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-capitalize(tolower(Glutamine_1_new))
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Glutamine1<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)

Glutamine_1<-read.table("d:/GOBP_GLUTAMINE_FAMILY_AMINO_ACID_CATABOLIC_PROCESS.v2023.2.Hs.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-capitalize(tolower(Glutamine_1_new))
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Glutamine2<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)

Glutamine_1<-read.table("d:/GOBP_GLUTAMINE_FAMILY_AMINO_ACID_METABOLIC_PROCESS.v2023.2.Hs.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-capitalize(tolower(Glutamine_1_new))
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Glutamine3<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)
remove(ref)
gc()
saveRDS(pbmc,"d:/CJJGSE109732-1.rds")
s.genes <- capitalize(tolower(cc.genes$s.genes))
g2m.genes <- capitalize(tolower(cc.genes$g2m.genes))

pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc))
remove(g2m.genes,s.genes)
gc()
saveRDS(pbmc,"d:/CJJGSE109732-1_cellcycle.rds")

library(sctransform)
pbmc <- SCTransform(pbmc, verbose = FALSE)
gc()
saveRDS(pbmc,"d:/CJJGSE109732-1_cellcycle_sctransform.rds")



Idents(pbmc)<-pbmc$newcluster
VlnPlot(pbmc,features = c("plasmacell"),sort = TRUE)
DotPlot(pbmc,features =c("plasmacell"))
FeaturePlot(pbmc,features =c("plasmacell"),max.cutoff = "q95") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
Idents(pbmc)<-pbmc$newcluster
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"),assay = "RNA")
library(viridis)

FeaturePlot(pbmc,features =c("Irf4"),order = TRUE,cols = c("lightgrey", "#B22222"))
FeaturePlot(pbmc,features =c("Kdm6b"),order = TRUE,cols = c("lightgrey", "#B22222"))
FeaturePlot(pbmc,features =c("Ly75"),order = TRUE,cols = c("lightgrey", "#B22222"))
FeaturePlot(pbmc,features =c("plasmacell"),max.cutoff = "q75",reduction = "lda2umap",order = TRUE,min.cutoff = "q25") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features =c("Irf4"),cols = c("lightgrey", "#B22222"))
FeaturePlot(pbmc,features =c("Kdm6b"),cols = c("lightgrey", "#B22222"))
FeaturePlot(pbmc,features =c("Ly75"),cols = c("lightgrey", "#B22222"))
FeaturePlot(pbmc,features =c("plasmacell"))+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
VlnPlot(pbmc,features = c("plasmacell"),sort = TRUE)
DotPlot(pbmc,features =c("plasmacell"))


# Idents(pbmc)<-pbmc$Cluster
# new.cluster.ids <- c("DZ_2", "Fraction3_2", "Fraction3_1", "Fraction3_3", "DZ_1", "DZ_3",
#                      "Fraction2_3", "Fraction2_2", "prePB","Fraction2_1","Fraction1_2","prePB","prePB","prePB","prePB")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, label = TRUE, pt.size = 0.5) + NoLegend()
# pbmc@meta.data$newnewcluster<-Idents(pbmc)
# pbmc<-RunLDA(pbmc,labels = pbmc$newnewcluster)
# pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:10)
# pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:10)
# DimPlot(pbmc,reduction = "ldaumap")
# DimPlot(pbmc,reduction = "ldatsne")
# VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
# DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
# FeaturePlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
# DimPlot(pbmc)
# library(CelliD)
# pathway<-read.table("D:/GSE109732_GeneExpression_DZ_Fraction1-3_GC-PB.txt",sep = "\t",header = TRUE)
# pathway<-pathway[!duplicated(pathway$Gene),]
# rownames(pathway)<-pathway$Gene
# pathway<-pathway[,-1]
# pathway<-as.list(pathway)
# pbmc<-RunMCA(pbmc)
# HGT_organoid_gs <- RunCellHGT(pbmc, pathways = pathway, n.features = 200,dims = 1:50)
# HGT_organoid_gs_prediction <- rownames(HGT_organoid_gs)[apply(HGT_organoid_gs, 2, which.max)]
# HGT_organoid_signif <- ifelse(apply(HGT_organoid_gs, 2, max)>2, yes = HGT_organoid_gs_prediction, "unassigned")
# pbmc$HGT_organoid_prediction <- HGT_organoid_signif
