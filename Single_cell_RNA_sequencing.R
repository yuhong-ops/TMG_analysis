
library("Seurat") 
library("dplyr") 
library("data.table") 
library("ggplot2") 
library("ggsci")
library("RColorBrewer") 
library("SingleR") ##单细胞类型自动注释 
library("celldex")  
##  多核并行
library(future)
plan(multicore, workers = 20)
options(future.globals.maxSize = 100000 * 1024^5)

#rm(list=ls())# 清除环境变量
#gc() ## 释放内存

setwd("D:/Documents/TMG_methylation_paper") 
# Meta_sc<-read.csv("thymoma_single_cell_seq/meta.csv") 
# rownames(Meta_sc)<-Meta_sc$NAME 
# 
# table(Meta_sc$organ__ontology_label) 
# 
# cluster_sc<-read.csv("thymoma_single_cell_seq/cluster.csv") 
# 
# Expr_sc<-fread("thymoma_single_cell_seq/expression.csv", sep = ',', header = TRUE) 
# rownames(Expr_sc)<-Expr_sc$V1 
# Expr_sc<-Expr_sc %>% dplyr::select(-V1) 
# 
# View(rownames(Expr_sc))  
# 
# 
# meta_clu<-cbind(Meta_sc,cluster_sc[,c(2,3)])  
# 
# ##初始化Seurat对象  
# sc_ob<-CreateSeuratObject(counts = Expr_sc, project = "thymoma_sc_project",  
#                           assay = "RNA",  
#                           names.field = 1,  
#                           names.delim = "_",  
#                           meta.data = Meta_sc)  
# sc_ob@meta.data$nCount_RNA  
# sc_ob@meta.data$nFeature_RNA  
# sc_ob@meta.data$percent.mt  
# 
# VlnPlot(sc_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)  
# 
# hist(colSums(sc_ob$RNA@data),  
#      breaks = 100, 
#      main = "Total expression before normalisation",  
#      xlab = "Sum of expression")  
# 
# ## scale data  
# all.genes <- rownames(sc_ob)  
# sc_ob <- ScaleData(sc_ob, features = all.genes)  
# 
# ## find variable feature 
# sc_ob <- FindVariableFeatures(sc_ob, selection.method = "vst", nfeatures = 2000) 
# 
# ## run PCA 
# sc_ob<- RunPCA(sc_ob, features = VariableFeatures(object = sc_ob)) 
# 
# ## 确定数据的维度 
# ElbowPlot(sc_ob) 
# 
# # 细胞聚类 
# sc_ob<-FindNeighbors(sc_ob,dims=1:10) 
# sc_ob<-FindClusters(sc_ob,resolution =0.5)
# ## 非线性降维
# sc_ob <- RunUMAP(sc_ob, dims = 1:10)
# dev.off()
# DimPlot(sc_ob, reduction = "umap")
# # 显示聚类标签
# DimPlot(sc_ob, reduction = "umap",label = TRUE)
# 
# ## 使用tSNE聚类
# sc_ob <- RunTSNE(sc_ob, dims = 1:10)
# DimPlot(sc_ob, reduction = "tsne", label=TRUE)
# 
# saveRDS(sc_ob,file="sc_ob.rds")



# 
# sc_ob<-readRDS("sc_ob.rds")
# 
# ## cell annoation
# #sc_ob_for_SingleR <- GetAssayData(sc_ob, slot="data")
# #hpca.se <- HumanPrimaryCellAtlasData()
# #sc_ob.hesc <- SingleR(test = sc_ob_for_SingleR , ref = hpca.se, labels = hpca.se$label.main)
# 
# 
# #seurat 和 singleR的table表
# #sc_ob@meta.data$labels <-sc_ob.hesc$labels
# #DimPlot(sc_ob, group.by = c("seurat_clusters", "labels"),reduction = "umap",label=TRUE)
# #DimPlot(sc_ob, group.by = c("labels"),reduction = "umap",label=TRUE)
# DimPlot(sc_ob, group.by = c("cell_type__ontology_label"),reduction = "umap",label=TRUE) + scale_color_npg()
# DimPlot(sc_ob, group.by = c("cell_type__custom"),reduction = "umap",label=TRUE) 
# FeaturePlot(sc_ob,features = c("NEFM"))

# ## subset thymus cells and reanalysis
# sc_thy = sc_ob[,sc_ob@meta.data$organ__ontology_label =="thymus"]
# all.genes.thy <- rownames(sc_thy)
# sc_thy <- ScaleData(sc_thy, features = all.genes.thy)
# sc_thy <- FindVariableFeatures(sc_thy, selection.method = "vst", nfeatures = 2000)
# sc_thy <- RunPCA(sc_thy, features = VariableFeatures(object = sc_thy))
# ElbowPlot(sc_thy)
# sc_thy<-FindNeighbors(sc_thy,dims=1:10)
# sc_thy<-FindClusters(sc_thy,resolution =0.5)
# sc_thy <- RunUMAP(sc_thy, dims = 1:10)
# saveRDS(sc_thy,file="sc_thy.rds")
sc_thy<-readRDS("sc_thy.rds")
# DimPlot(sc_thy, reduction = "umap")
# DimPlot(sc_thy, reduction = "umap",label = TRUE)
DimPlot(sc_thy, group.by = "cell_type__custom",reduction = "umap",label=TRUE)
# 
# ##
DimPlot(sc_thy, group.by = c("cell_type__ontology_label"),reduction = "umap",label=TRUE) + scale_color_npg()
# DimPlot(sc_thy, group.by = c("cell_type__custom"),reduction = "umap",label=TRUE) 


### DE analysis
sc_thy@meta.data$cell_type__custom2<-sc_thy@meta.data$cell_type__custom
sc_thy@meta.data$cell_type__custom2[sc_thy@meta.data$cell_type__custom %in% c("nmTEC","mTEC (I)","mTEC (II)")]<-"mTEC"


sc_thy@meta.data$cell_type__custom3<-sc_thy@meta.data$cell_type__custom
sc_thy@meta.data$cell_type__custom3[sc_thy@meta.data$cell_type__custom %in% c("nmTEC","mTEC (I)","mTEC (II)")]<-"mTEC"
sc_thy@meta.data$cell_type__custom3[!(sc_thy@meta.data$cell_type__custom %in% c("nmTEC","mTEC (I)","mTEC (II)"))]<-"non-mTEC"


DimPlot(sc_thy, group.by = c("cell_type__custom3"),reduction = "umap",label=TRUE) + scale_color_npg()

## DEG:nmTEC 与 其他细胞对比
#deg = FindMarkers(sc_thy,ident.1 = 'nmTEC',
#                    group.by="cell_type__custom2", min.pct = 0.01)
# #saveRDS(deg,file="deg.rds")
deg<-readRDS("deg.rds")

## DEG2:mTEC 与 其他细胞对比，ming.pct=0.01
## DEG3:mTEC 与 其他细胞对比，ming.pct=0.001

deg2 = FindMarkers(sc_thy,ident.1 = 'mTEC',
                  group.by="cell_type__custom2", ming.pct = 0.01)

# 
#saveRDS(deg2,file="deg3.rds")

deg2<-readRDS("deg2.rds")

deg3<-readRDS("deg3.rds")

#deg4 = FindMarkers(sc_thy,ident.1 = 'mTEC',
#                   group.by="cell_type__custom2", ming.pct = 0.0001)

?FindMarkers
#saveRDS(deg4,file="deg4.rds")
deg2<-deg4

deg3[row.names(deg3)=="NEFM",]
deg3[row.names(deg3)=="RYR3",]
deg3[row.names(deg3)=="CHRNA1",]
deg3[row.names(deg3)=="GABRA5",]

deg2[row.names(deg2)=="NEFM",]
deg[row.names(deg)=="NEFM",]
deg2[row.names(deg2)=="RYR3",]
deg[row.names(deg)=="RYR3",]
deg2[row.names(deg2)=="GABRA5",]
deg[row.names(deg)=="GABRA5",]
deg2[row.names(deg2)=="DNMT1",]
deg[row.names(deg)=="DNMT1",]
deg2[row.names(deg2)=="UHRF1",]
deg[row.names(deg)=="UHRF1",]
deg2[row.names(deg2)=="AIRE",]

deg2[row.names(deg2)=="HDAC1",]
deg2[row.names(deg2)=="RYR1",]
deg2[row.names(deg2)=="RYR3",]
deg2[row.names(deg2)=="CHRNA1",]
deg2[row.names(deg2)=="TET1",]
deg2[row.names(deg2)=="TET2",]
deg2[row.names(deg2)=="TET3",]

deg[row.names(deg2)=="HDAC1",]
deg[row.names(deg2)=="TET1",]
deg[row.names(deg2)=="TET2",]
deg[row.names(deg2)=="TET3",]
deg[row.names(deg)=="MTF2",]


deg2[row.names(deg2)=="DNMT3A",]
deg2[row.names(deg2)=="MECP2",]
deg2[row.names(deg2)=="MBD1",]
deg2[row.names(deg2)=="MBD2",]
deg2[row.names(deg2)=="MBD3",]
deg2[row.names(deg2)=="MBD4",]
deg2[row.names(deg2)=="AIRE",]
deg[row.names(deg)=="AIRE",]
deg2[row.names(deg2)=="MECP2",]
deg2[row.names(deg2)=="MTF2",]
deg2[row.names(deg2)=="OGT",]
deg2[row.names(deg2)=="GATA3",]

write.csv(deg3,"single cell RNA seq mTEC vs non-mTEC.csv")

DotPlot(sc_thy,features = c("NEFM","RYR3","GABRA5","CHRNA1","NEFL"),assay='RNA',group.by ="cell_type__custom2")
DotPlot(sc_thy,features = c("DNMT1","UHRF1","CTCF","AIRE"),assay='RNA',group.by ="cell_type__custom2")

DotPlot(sc_thy,features = unique(c("DNMT1","UHRF1","HDAC1",'HDAC2',"CTCF",'SUZ12','RBBP4',"NEFM","RYR3","GABRA5")),assay='RNA',
        
        group.by ="cell_type__custom3",dot.scale = 16,cols=c("#4E9AD4","#DB5434"))

DotPlot(sc_thy,features = c("CHRNA1","RYR3","UHRF1","AIRE"),assay='RNA',group.by ="cell_type__custom2")

##画图表达量的barplot

VlnPlot(sc_thy,features="DNMT1",group.by="cell_type__custom3")
VlnPlot(sc_thy,features="UHRF1",group.by="cell_type__custom3")
VlnPlot(sc_thy,features="HDAC1",group.by="cell_type__custom3")
VlnPlot(sc_thy,features="AIRE",group.by="cell_type__custom3")


DoHeatmap(sc_thy, features = c("DNMT1","UHRF1","HDAC1","AIRE"), group.by ="cell_type__custom3", cells = 1:500, size = 4,
          angle = 90) + NoLegend()

## 将有意义的甲基化基因表达量提取出来
mymatrix <- as.data.frame(sc_thy@assays$RNA@data[c("DNMT1","UHRF1","HDAC1"),])
mymatrix2<-t(mymatrix)%>%as.data.frame()

mymatrix2[,ncol(mymatrix2)+1]<-sc_thy$cell_type__custom3
colnames(mymatrix2)[ncol(mymatrix2)] <- "Group"

write.csv(mymatrix2,"methylation_gene_sc_data.csv")






mymatrix2 %>%
  group_by(Group) %>%
  summarise_at(vars(DNMT1,UHRF1), list(name = mean))

library(ggpubr)
library(ggprism)


library(ggplot2)
p1<- ggplot2::ggplot(mymatrix2,aes(x=Group,y=DNMT1,fill=Group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange'))


ggplot(mymatrix2,aes(Group,DNMT1,fill=Group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", width = 0.3)+
  theme_bw()

?stat_summary

p2<- ggplot2::ggplot(mymatrix2,aes(x=Group,y=UHRF1,fill=Group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange'))

p3<- ggplot2::ggplot(mymatrix2,aes(x=Group,y=HDAC1,fill=Group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange'))
p1
p2
p3


ggboxplot(mymatrix2, x = "Group", y = "DNMT1",
                color = "Group", palette = "npg",
                add = "none")


p2<-ggboxplot(mymatrix2, x = "Group", y = "UHRF1",
          color = "Group", palette = "npg",
          add = "none")

p2<-ggboxplot(mymatrix2, x = "Group", y = "DNMT1",
          color = "Group", palette = "npg",
          add = "none")

ggboxplot(mymatrix2, x = "Group", y = "HDAC1",
          color = "Group", palette = "npg",
          add = "none")



p2 + stat_compare_means(method = "t.test", method.args = list(alternative = "less"),
                        label =  "p.signif", 
                        label.x = 1.5, label.y = 10)






library(org.Hs.eg.db)
library(clusterProfiler)
k<-keys(org.Hs.eg.db,keytype='ENSEMBL') 
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,columns=c("ENTREZID",'SYMBOL'), keytype = 'ENSEMBL') 

 
Go<-function(a){ge<-en2ENSE %>% filter(SYMBOL %in% a) 
                 ge_path<-enrichGO(gene          = ge$ENTREZID, 
                            OrgDb         = org.Hs.eg.db, 
                            ont           = "ALL", 
                            pAdjustMethod = "fdr", 
                            pvalueCutoff  = 1, 
                            qvalueCutoff  = 1,
                            minGSSiz=1,
                            maxGSSiz=100,
                            readable      = TRUE)  
                 return(ge_path)
}

DE2_list<-deg2 %>% mutate(SYMBOL=rownames(deg2)) %>% left_join(en2ENSE,by="SYMBOL") 
 
DE2_gene_list_up<-rownames(deg2)[(deg2$p_val_adj<0.05)&(deg2$avg_log2FC>0)] 
DMP2_ge_up_path<-Go(DE2_gene_list_up)
View(DMP2_ge_up_path@result) 

DE2_gene_list_down<-rownames(deg2)[(deg2$p_val_adj<0.01)&(deg2$avg_log2FC<0)] 
DMP2_ge_down_path<-Go(DE2_gene_list_down)
View(DMP2_ge_down_path@result) 
write.csv(DMP2_ge_down_path@result,"DMP2_ge_down_path2.csv") 
DMP2_ge_down_path@result<-read.csv("DMP2_ge_down_path.csv")


DE_list<-deg %>% mutate(SYMBOL=rownames(deg)) %>% left_join(en2ENSE,by="SYMBOL")

DE_gene_list_down<-rownames(deg)[(deg$p_val_adj<0.05)&(deg$avg_log2FC<0)]
DMP_ge_down_path<-Go(DE_gene_list_down)
View(DMP_ge_down_path@result)



library(stringr)

meth_path<-DMP2_ge_down_path
meth_path@result<-DMP2_ge_down_path@result %>% filter((str_detect(DMP2_ge_down_path@result$Description,pattern = "methy") |
                                                        str_detect(DMP2_ge_down_path@result$Description,pattern = "epigene"))
                                                        & DMP2_ge_down_path@result$pvalue<0.05  ) 
View(meth_path@result)
p1<-barplot(meth_path,color="pvalue") + scale_fill_continuous(high = "#4E9AD4", 
                      low = "#DB5434", guide= guide_colorbar(reverse = TRUE))
acety_path<-DMP2_ge_down_path
acety_path@result<-DMP2_ge_down_path@result %>% filter((str_detect(DMP2_ge_down_path@result$Description,pattern = "acetyla") )
                                                      & DMP2_ge_down_path@result$pvalue<0.05  )
p2<-barplot(acety_path,color="pvalue") + scale_fill_continuous(high = "#4E9AD4", 
                                                          low = "#DB5434", guide= guide_colorbar(reverse = TRUE))


ggpubr::ggarrange(p1,p2,ncol=2)



##gsea analysis

gene_lis<-  DE2_list %>% filter(p_val_adj<0.05) %>% select(ENTREZID,avg_log2FC) %>% distinct(ENTREZID,.keep_all = T) %>% na.omit()
gene_list <- gene_lis[,2]
names(gene_list)<-gene_lis[,1]
gene_list = sort(gene_list, decreasing = TRUE)

GO_gse_entrez_sc <- gseGO(geneList     = gene_list,
                      ont          = "ALL", 
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "ENTREZID",
                      pvalueCutoff = 1,
                      minGSSize = 1,
                      maxGSSize = 500)   #实际为padj
GO_gse_sc <- DOSE::setReadable(GO_gse_entrez_sc, 
                           OrgDb=org.Hs.eg.db,
                           keyType='ENTREZID')#转化id 
View(GO_gse_sc@result)

deg2[row.names(deg2)=="CHRNA1",]
?gseGO

FindMarkers


### subset stromal cells and reanalysis

sc_stro = sc_thy[,sc_thy@meta.data$cell_type__ontology_label =="stromal cell"]
sc_mtec = sc_thy[,sc_thy@meta.data$cell_type__custom3 =="mTEC"]

all.genes.stro <- rownames(sc_stro)
sc_stro <- ScaleData(sc_stro, features = all.genes.stro)
sc_stro <- FindVariableFeatures(sc_stro, selection.method = "vst", nfeatures = 2000)
sc_stro<- RunPCA(sc_stro, features = VariableFeatures(object = sc_stro))
ElbowPlot(sc_stro)
sc_stro<-FindNeighbors(sc_stro,dims=1:7)
sc_stro<-FindClusters(sc_stro,resolution =0.5)
sc_stro <- RunUMAP(sc_stro, dims = 1:7)
DimPlot(sc_stro, reduction = "umap")
DimPlot(sc_stro, reduction = "umap",label = TRUE)
DimPlot(sc_stro, group.by = "cell_type__custom",reduction = "umap",label=TRUE) + scale_color_npg()
sc_stro@meta.data

DotPlot(sc_stro,features = c("NEFM","RYR3","GABRA5","CHRNA1","NEFL"),assay='RNA',group.by ="cell_type__custom")
DotPlot(sc_stro,features = c("DNMT1","UHRF1","CTCF"),assay='RNA',group.by ="cell_type__custom")
 
 
# Stromal cells  
#  endothelial cells (positive for PECAM1/CD31, VWF),   
#  normal fibroblasts (FN1, EGFL6),  
#  tumor-associated fibroblasts (TAFs; PDGFRA, ADH1B)  
#  thymic epithelial cells (TECs; KRT19, S100A14)
#   cTEC (CCL25, PSMB11)
#   mTEC (CCL19, KRT7) clusters
#     mTEC(I) KRT15,  IFI27; 
#     mTEC(II) CLDN4, KRT7; 
#     nmTECs  GABRA5, MAP2, NEFL, NEFM, SOX15, TF. 
#       nmTECs: KRT6 and GABRA5 as marker genes 
FeaturePlot(sc_thy,features = c("AIRE"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("NEFM"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("RYR3"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("GABRA5"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("CHRNA1"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("RYR1"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("TTN"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("KRT5"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("KRT6A"),pt.size=2,order=T) 
FeaturePlot(sc_thy,features = c("KRT8"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("KRT19"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("S100A14"),pt.size=2,order=T)

FeaturePlot(sc_thy,features = c("DNMT1"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("UHRF1"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("DNMT3A"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("DNMT3B"),pt.size=2,order=T)


FeaturePlot(sc_thy,features = c("TET1"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("TET2"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("TET3"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("DNMT1"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("UHRF1"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("DNMT3A"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("DNMT3B"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("AIRE"),pt.size=2,order=T)



FeaturePlot(sc_stro,features = c("TET1"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("TET2"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("TET3"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("CHRNA1"),pt.size=2,order=T)

FeaturePlot(sc_thy,features = c("CCL19"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("CK19"),pt.size=2,order=T)


FeaturePlot(sc_stro,features = c("NEFM"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("RYR3"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("GABRA5"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("KRT8"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("KRT19"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("S100A14"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("CCL19"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("CCL25"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("PSMB11"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("MTS5"),pt.size=2,order=T)


FeaturePlot(sc_stro,features = c("KRT7"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("KRT5"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("KRT6A"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("MUC1"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("MUC1"),pt.size=2,order=T)

FeaturePlot(sc_stro,features = c("DNMT1"),pt.size=2,order=T)


### 髓质标记物

FeaturePlot(sc_stro,features = c("CD40"),pt.size=2,order=T)
FeaturePlot(sc_stro,features = c("AIRE"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("CD40"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("AIRE"),pt.size=2,order=T)
FeaturePlot(sc_thy,features = c("AIRE"),pt.size=2,order=T)

FeaturePlot(sc_thy,features = c("CD4"),pt.size=2,order=T)



## 使用tSNE聚类
sc_ob <- RunTSNE(sc_ob, dims = 1:10)
DimPlot(sc_ob, reduction = "tsne", label=TRUE)

saveRDS(sc_ob,file="sc_ob.rds")
sc_ob<-readRDS("sc_ob.rds")

head(sc_ob$RNA@data[,1:5])

sc_ob[["percent.mt"]] <- PercentageFeatureSet(sc_ob, pattern = "^MT-")
head(sc_ob@meta.data, 5)


thymus_meta_clu <-meta_clu %>% filter(organ__ontology_label == "thymus")
thymus_cell_id<-thymus_meta_clu$NAME

## thymus expression count matrix
Expr_sc_th<-Expr_sc[,thymus_cell_id,with=F] 
rownames(Expr_sc_th)<-rownames(Expr_sc)
rm(Expr_sc) #删除原始count，节省空间

which(rownames(Expr_sc_th)=="CHRNA1")
which(rownames(Expr_sc_th)=="NEFM")
which(rownames(Expr_sc_th)=="RYR3")
which(rownames(Expr_sc_th)=="GABRA5")
which(rownames(Expr_sc_th)=="SPOCK1")
which(rownames(Expr_sc_th)=="MCOLN2")
which(rownames(Expr_sc_th)=="TTN")
which(rownames(Expr_sc_th)=="RYR1")

thymus_meta_clu$CHRNA1<-as.vector(t(Expr_sc_th[5185,]))
thymus_meta_clu$NEFM<-as.vector(t(Expr_sc_th[15036,]))
thymus_meta_clu$NEFM_scale<-rescale(as.vector(t(Expr_sc_th[15036,])))
thymus_meta_clu$RYR3<-as.vector(t(Expr_sc_th[25328,]))
thymus_meta_clu$GABRA5<-as.vector(t(Expr_sc_th[25257,]))
thymus_meta_clu$SPOCK1<-as.vector(t(Expr_sc_th[10537,]))
thymus_meta_clu$MCOLN2<-as.vector(t(Expr_sc_th[1395,]))
thymus_meta_clu$TTN<-as.vector(t(Expr_sc_th[5240,]))
thymus_meta_clu$RYR1<-as.vector(t(Expr_sc_th[32002,]))

## umap plot
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(CHRNA1)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=CHRNA1),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(NEFM)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=NEFM),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(RYR3)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=RYR3),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(GABRA5)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=GABRA5),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(SPOCK1)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=SPOCK1),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")

  
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(TTN)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=TTN),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_t<-thymus_meta_clu %>% arrange(RYR1)
  ggplot(data=thymus_meta_clu_t) + 
    geom_point(aes(x=X,y=Y,color=RYR1),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  

## stromal expression count matrix

  
  thymus_meta_clu_stromal <-thymus_meta_clu %>% 
    filter(cell_type__ontology_label == "stromal cell" & cell_type__custom != "Endothelial cell")
  thymus_cell_id_stromal<-thymus_meta_clu_stromal$NAME
  
  ## thymus expression count matrix
  Expr_sc_th_stro<-Expr_sc_th[,thymus_cell_id_stromal,with=F] 
  rownames(Expr_sc_th_stro)<-rownames(Expr_sc_th)
  
  
  ### thymus_meta_clu_stromal$CHRNA1<-as.vector(t(Expr_sc_th_stro[5185,]))
  
  ## umap plot

  ggplot(data=thymus_meta_clu_stromal) + 
    geom_point(aes(x=X,y=Y,color=cell_type__custom),size=1.2,alpha=1)+
    theme_classic(base_size=16)  + scale_color_npg()
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(CHRNA1)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=CHRNA1),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(NEFM)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=NEFM),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(RYR3)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=RYR3),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(GABRA5)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=GABRA5),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(SPOCK1)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=SPOCK1),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(TTN)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=TTN),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")
  
  
  thymus_meta_clu_stromal_t<-thymus_meta_clu_stromal %>% arrange(RYR1)
  ggplot(data=thymus_meta_clu_stromal_t) + 
    geom_point(aes(x=X,y=Y,color=RYR1),size=2,alpha=1)+
    theme_classic(base_size=16)  +scale_colour_gradient(low="lightgrey", high="blue")


##

sc_mtec = sc_thy[,sc_thy@meta.data$cell_type__custom3 =="mTEC"]
table(sc_mtec@meta.data$biosample_id)
sc_mtec

counts_mtec<-GetAssayData(object =sc_mtec, layer = "counts")

counts_mtx<-as.data.frame(counts_mtec[c("NEFM","RYR3","GABRA5","UHRF1","DNMT1","AIRE"),])
View(counts_mtx)

FeaturePlot(sc_thy,features = c("DNMT1"),pt.size=2,order=T)
FeaturePlot(sc_mtec,features = c("DNMT1"),pt.size=2,order=T)
FeaturePlot(sc_mtec,features = c("DNMT1"),pt.size=2,order=T)
FeaturePlot(sc_mtec,features = c("AIRE"),pt.size=2,order=T)

library(muscat)
sc_mtec$biosample_id

sc_mtec_sc<-as.SingleCellExperiment(sc_mtec)

pb <- aggregateData(sc_mtec_sc,
                    fun = "sum",
                    by = c("biosample_id"))
assay(pb)[c("NEFM","RYR3","GABRA5","UHRF1","DNMT1"),]



### DE3 上调和下调基因列表

DE_number<-data.frame(Group=c('Up','Down'),
                      DE_gene_number=c(length(DE3_gene_list_up),length(DE3_gene_list_down)))
p <- ggplot(DE_number, aes(x=Group, y=DE_gene_number,fill=Group)) +
  geom_bar(stat="identity",  colour="black", position=position_dodge(),width = 0.6)+
  theme_classic() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))
p + scale_fill_npg()
