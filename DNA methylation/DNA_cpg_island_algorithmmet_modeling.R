library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(reshape)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(IMA2)
#library(IMA)
library(ggplot2)
library(ggsci)

#source("./script/TCGA_methylation.R") 
#source("./script/GSE94769_methylation.R") 

# extract ensemble and entrze gene ID mapping files 
k<-keys(org.Hs.eg.db,keytype='ENSEMBL') 
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k, 
                               columns=c("ENTREZID",'SYMBOL'),  
                               keytype = 'ENSEMBL')  

#setwd("C:/Users/yuhon/Documents/TMG_methylation_paper") #laptop  

data(Islands.UCSC) 
data(Locations)  
data(Other)  
?Other
 
table(Other$UCSC_RefGene_Group)  
#methylation island data from ill450hg19 package  

isl_UC<- data.frame(Islands.UCSC) %>% 
  mutate(me_id = row.names(.)) 

#View(data.frame(Other)) 
me_ge_info<- data.frame(Other) %>% 
  dplyr::select("UCSC_RefGene_Name","UCSC_RefGene_Accession","UCSC_RefGene_Group") %>%   
  mutate(me_id = row.names(.))    

isl_info<-data.frame(Locations) %>%   
  mutate(me_id = row.names(.)) %>%    
  inner_join(isl_UC,by="me_id") %>%  
  inner_join(me_ge_info,by="me_id")   

#View(data.frame(Other)) 
 
# extract methylation island methylation loci data,extract the cpg island with at least 3 methylation loci
## 导入所有甲基化位点的normalization之后的beta值
myNorm_met<-read.csv("./analy_data/tcga_met_normalization_TCGA_database_20220813.csv",check.names = F) 
row.names(myNorm_met)<-myNorm_met[,1] 
myNorm_met<-myNorm_met[,-1] %>% as.matrix() 
myNorm_met<-na.omit(myNorm_met)

# 导入每个甲基化位点的详细注释信息
me_da_t<-myNorm_met %>% data.frame(check.names=F) %>% mutate(me_id=row.names(.))   
me_da_anno<-isl_info %>% right_join(me_da_t, by="me_id")
me_da_anno<-me_da_anno %>% mutate(Relation_to_Island_2=Relation_to_Island)
me_da_anno$Relation_to_Island_2[me_da_anno$Relation_to_Island_2 %in% c("N_Shore","S_Shore")]<-"Shore"
me_da_anno$Relation_to_Island_2[me_da_anno$Relation_to_Island_2 %in% c("N_Shelf","S_Shelf")]<-"Shelf" 
me_da_anno$Relation_to_Island_2[me_da_anno$Relation_to_Island_2 == "Island"]<-"CpG island"  
me_da_anno$Relation_to_Island_2[me_da_anno$Relation_to_Island_2 == "OpenSea"]<-"Open sea"  

## 柱状图：位于基因的shore，shelf，island，opensea区域的甲基化位点数量
library(ggsci)
p1<-ggplot(data=me_da_anno,aes(x=Relation_to_Island_2))+
  geom_bar(aes(fill=Relation_to_Island_2))+ theme_classic()+ scale_fill_npg()
p1

# GSE 数据集
me_da_g<-myNorm_met_gse %>% data.frame(check.names=F) %>% mutate(me_id=row.names(.))
 me_da_isanno<-isl_info %>% right_join(me_da_g, by="me_id")
 me_da_isanno<-me_da_isanno %>% mutate(Relation_to_Island_2=Relation_to_Island)
 me_da_isanno$Relation_to_Island_2[me_da_isanno$Relation_to_Island_2 %in% c("N_Shore","S_Shore")]<-"Shore"
 me_da_isanno$Relation_to_Island_2[me_da_isanno$Relation_to_Island_2 %in% c("N_Shelf","S_Shelf")]<-"Shelf"
 me_da_isanno$Relation_to_Island_2[me_da_isanno$Relation_to_Island_2 == "Island"]<-"CpG island"
 me_da_isanno$Relation_to_Island_2[me_da_isanno$Relation_to_Island_2 == "OpenSea"]<-"Open sea"

 #维恩图，两个队列的重复
 length(me_da_isanno$me_id)
 length(me_da_anno$me_id)
 intersect(me_da_isanno$me_id,me_da_anno$me_id)
 length(me_da_anno$me_id)-length(intersect(me_da_isanno$me_id,me_da_anno$me_id))
 length(me_da_isanno$me_id)-length(intersect(me_da_isanno$me_id,me_da_anno$me_id))


library(limma)

# ## GSE dataset, 画图，所有甲基化位点所在区域
# library(ggsci)
# p2<-ggplot(data=me_da_isanno,aes(x=Relation_to_Island_2))+
#   geom_bar(aes(fill=Relation_to_Island_2))+ theme_classic()+ scale_fill_npg()
# 
# met_region<-rbind(data.frame(Group="Primary cohort",Region=me_da_anno$Relation_to_Island_2),
#                   data.frame(Group="Verification cohort",Region=me_da_isanno$Relation_to_Island_2))
# 
# ggplot(data=met_region,aes(x=Region))+
#   geom_bar(aes(fill=Group),position=position_dodge(width=0.9))+ theme_classic()+ scale_fill_npg()
# 
# table(me_da_isanno$Relation_to_Island) 
# sum(table(me_da_isanno$Relation_to_Island))


##  筛选出属于CpG岛的甲基化位点及注释
me_da_isl<-isl_info %>% right_join(me_da_t, by="me_id") %>%   
  filter(Islands_Name!="" & Relation_to_Island =="Island")

  ##  筛选出每个CpG岛甲基化位点的数量至少为3个的CpG岛
me_num_pisl<-me_da_isl %>% group_by(Islands_Name) %>% summarise(count = length(Islands_Name)) %>%   
  filter(count>=3) %>% arrange(Islands_Name)

  ### 按CpG岛的名字排序
me_da_isl<-me_da_isl %>% filter(Islands_Name %in% me_num_pisl$Islands_Name) %>% arrange(Islands_Name) 

###与启动子相关的CpG岛的ID（TSS1500，TSS200， 5’UTR）
me_num_pisl_id<-me_da_isl %>% group_by(Islands_Name) %>% summarise(count = length(Islands_Name),
                                                                UCSC_RefGene_Group= paste(UCSC_RefGene_Group,collapse= ";")) %>%   
  arrange(Islands_Name) %>%filter(str_detect(.$UCSC_RefGene_Group, pattern = "TSS1500")|
                                                         str_detect(.$UCSC_RefGene_Group, pattern = "TSS200") |
                                                         str_detect(.$UCSC_RefGene_Group, pattern = "5'UTR")) %>%
  pull(Islands_Name)


#### 每个sample 每个CpG岛 的beta值的均值和中位数（TCGA），运行慢！！！ 
me_da_isl_mean_beta<-me_da_isl%>%group_by(Islands_Name)%>%summarise_at(vars(colnames(me_da_isl)[10:130]),mean)  
example<-me_da_isl_mean_beta[,2]%>% data.frame()
head(example)
ggplot(example,aes(x = TCGA.3G.AB0O)) + 
  geom_density(color = "black", fill = "grey") +   
  theme_classic()
me_da_isl_median_beta<-me_da_isl%>%group_by(Islands_Name)%>%summarise_at(vars(colnames(me_da_isl)[10:130]),median)  

#me_da_isl_MG<-me_da_isl %>% dplyr::select(1:10,meta_data$patient_id[meta_data$MG=="YES"]) %>%  
#  melt(,id.vars=colnames(me_da_isl)[1:10]) %>%group_by(Islands_Name)%>%summarise_at(vars(value),median)  

#me_da_isl_NO<-me_da_isl %>% dplyr::select(1:10,meta_data$patient_id[meta_data$MG=="NO"]) %>%  
#  melt(,id.vars=colnames(me_da_isl)[1:10]) %>%group_by(Islands_Name)%>%summarise_at(vars(value),median)  

me_da_isl_MG_mean<-me_da_isl %>% dplyr::select(1:10,meta_data$patient_id[meta_data$MG=="YES"]) %>%   
  melt(,id.vars=colnames(me_da_isl)[1:10]) %>%group_by(Islands_Name)%>%summarise_at(vars(value),mean)  

me_da_isl_NO_mean<-me_da_isl %>% dplyr::select(1:10,meta_data$patient_id[meta_data$MG=="NO"]) %>%  
  melt(,id.vars=colnames(me_da_isl)[1:10]) %>%group_by(Islands_Name)%>%summarise_at(vars(value),mean) 


##计算每个CpG岛的差异表达分析（运行慢！！，后面直接导入）

# ## TCGA 
# me_list<-list() 
# Sum_table<-data.frame() 
# j<-1 
# 
# for (i in c(1:length(me_num_pisl$Islands_Name))) { 
#   me_da_t<-me_da_isl[(j:((me_num_pisl$count[i])+j-1)),] 
#   
#   me_da_t2<-me_da_t %>% data.table::melt(colnames(.)[1:10],colnames(.)[11:length(colnames(.))], 
#                              variable.name="patient_id",value.name="me_beta") %>% arrange(pos) %>%  
#     left_join(met_meta, by=c("patient_id") )%>%  
#     mutate(me_isl_id=paste0(me_id,"-",patient_id)) 
#   
#   expr_isl_t<-me_da_t2 %>% dplyr::select(me_isl_id,me_beta) %>%  
#     spread(key=me_isl_id,value=me_beta) %>%  
#     dplyr::select(me_da_t2$me_isl_id) 
#   
#   met_loci<-factor(me_da_t2$me_id) 
#   MG<-factor(me_da_t2$MG,levels=c("NO","YES"))
#   design<-model.matrix(~met_loci+MG)
# 
#   fit <- lmFit(expr_isl_t, design)
#   fit <- eBayes(fit)
#   topTable<-topTable(fit, coef="MGYES") %>% 
#     mutate(chr=me_da_t$chr[1], 
#            Islands_Name= me_da_t$Islands_Name[1], 
#            me_id = paste(me_da_t$me_id,collapse= ";"), 
#            UCSC_RefGene_Name= paste(me_da_t$UCSC_RefGene_Name,collapse= ";"),
#            me_count = me_num_pisl$count[i])
#   
#   logFC<-log2(mean(me_da_t2$me_beta[me_da_t2$MG=="YES"])/mean(me_da_t2$me_beta[me_da_t2$MG=="NO"]))
#   topTable$logFC<-logFC
#   
#   Sum_table<-rbind(Sum_table,topTable)
#   j<-j+me_num_pisl$count[i]
# }
# 
# read CPG island tcga DE analysis result 
# Sum_table_tcga<-Sum_table %>% dplyr::select(Islands_Name,chr,me_id,me_count,UCSC_RefGene_Name,
#                                   everything())
# length(Sum_table_tcga$Islands_Name)
# write.csv(Sum_table_tcga,"CPG_island_differential_linear_model_TCGA_new.csv")


# ### #### ##################################  below only NFLM
# #chr8:24770908-24772547


## 手动选取某几个甲基化位点作为CpG island 进行DE 分析
## CpG_is_DE的函数的输入是

CpG_is_DE<-function(x){
  me_list<-list()
  Sum_table<-data.frame()
  me_da_t<-me_da_anno %>% filter(me_id %in% c("cg21552300","cg18003659"))
  
  me_da_t2<-me_da_t %>% data.table::melt(colnames(.)[1:9],colnames(.)[10:(length(colnames(.))-1)],
                                         variable.name="patient_id",value.name="me_beta") %>% arrange(pos) %>%
    left_join(met_meta, by=c("patient_id") )%>%
    mutate(me_isl_id=paste0(me_id,"-",patient_id))
  
  expr_isl_t<-me_da_t2 %>% dplyr::select(me_isl_id,me_beta) %>%
    spread(key=me_isl_id,value=me_beta) %>%
    dplyr::select(me_da_t2$me_isl_id)
  
  met_loci<-factor(me_da_t2$me_id)
  
  MG<-factor(me_da_t2$MG,levels=c("NO","YES"))
  design<-model.matrix(~met_loci+MG)
  design<-model.matrix(~MG)
  
  logFC<-log2(mean(me_da_t2$me_beta[me_da_t2$MG=="YES"])/mean(me_da_t2$me_beta[me_da_t2$MG=="NO"]))
  
  fit <- lmFit(expr_isl_t, design)
  fit <- eBayes(fit)
  topTable<-topTable(fit, coef="MGYES") %>%
    mutate(chr=me_da_t$chr[1],
           Islands_Name= me_da_t$Islands_Name[1],
           me_id = paste(me_da_t$me_id,collapse= ";"),
           UCSC_RefGene_Name= paste(me_da_t$UCSC_RefGene_Name,collapse= ";"),
           me_count = 15)
  topTable$logFC<-logFC
}
CpG_is_DE


me_list<-list()
Sum_table<-data.frame()
me_da_t<-me_da_anno %>% filter(me_id %in% c("cg01083589","cg02451888", "cg03607993"))

me_da_t2<-me_da_t %>% data.table::melt(colnames(.)[1:9],colnames(.)[10:(length(colnames(.))-1)],
                                       variable.name="patient_id",value.name="me_beta") %>% arrange(pos) %>%
  left_join(met_meta, by=c("patient_id") )%>%
  mutate(me_isl_id=paste0(me_id,"-",patient_id))

expr_isl_t<-me_da_t2 %>% dplyr::select(me_isl_id,me_beta) %>%
  spread(key=me_isl_id,value=me_beta) %>%
  dplyr::select(me_da_t2$me_isl_id)

met_loci<-factor(me_da_t2$me_id)

MG<-factor(me_da_t2$MG,levels=c("NO","YES"))
design<-model.matrix(~met_loci+MG)
View(design)
design<-model.matrix(~MG)

logFC<-log2(mean(me_da_t2$me_beta[me_da_t2$MG=="YES"])/mean(me_da_t2$me_beta[me_da_t2$MG=="NO"]))

fit <- lmFit(expr_isl_t, design)
fit <- eBayes(fit)
topTable<-topTable(fit, coef="MGYES") %>%
  mutate(chr=me_da_t$chr[1],
         Islands_Name= me_da_t$Islands_Name[1],
         me_id = paste(me_da_t$me_id,collapse= ";"),
         UCSC_RefGene_Name= paste(me_da_t$UCSC_RefGene_Name,collapse= ";"),
         me_count = 15)
topTable$logFC<-logFC

## 导入CpG岛差异表达分析结果
CPG_tcga<-read.csv("CPG_island_differential_linear_model_TCGA_new.csv",header=T)
CPG_tcga$adj.P.Val<-p.adjust(CPG_tcga$P.Value,method = "fdr")
CPG_tcga[CPG_tcga$adj.P.Val<0.05,]

CPG_tcga <- CPG_tcga %>% dplyr::arrange(P.Value)

## 将差异表达分析的CpG岛关联的基因的列表提取出来，去掉重复项目
CPG_gene_list<-CPG_tcga$UCSC_RefGene_Name %>% str_split(pattern = ";")
CPG_gene_list_c<-vector()
 for (i in c(1:length(CPG_gene_list))){
   CPG_gene_t<-unique(CPG_gene_list[[i]])
   CPG_gene_list_c[i]<-str_c(CPG_gene_t,collapse = ";")
 }
CPG_tcga$UCSC_RefGene_Name <-CPG_gene_list_c %>% str_remove(pattern="^;") %>% str_remove(pattern=";$")



#### TCGA cohort 位于基因组不同位置（island, shore, opensea and shelf）的CpG位点的甲基化程度的密度图   
#TCGA CPG island, shore, opensea and shelf density plot (beta value) of MG and NMG patients
table(isl_info$Relation_to_Island)
Island_met_loci<-isl_info %>% filter(Relation_to_Island == "Island") %>% pull(me_id)
OpenSea_met_loci<-isl_info %>% filter(Relation_to_Island == "OpenSea") %>% pull(me_id)
Shelf_met_loci<-isl_info %>% filter(Relation_to_Island %in% c("S_Shelf","N_Shelf")) %>% pull(me_id)
Shore_met_loci<-isl_info %>% filter(Relation_to_Island %in% c("S_Shore","N_Shore")) %>% pull(me_id)

p<-ggplot(CPG_tcga, aes(x = AveExpr)) + 
  geom_density(color = "black", fill = "grey") +   
   theme_classic()
p+ scale_color_npg()
 
myNorm_met_t<-data.table(myNorm_met)
Met_tcga_norm_mean<-data.frame(AveExpr=rowMeans(myNorm_met))
Met_tcga_norm_mean_NMG<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                       dplyr::select(meta_data$patient_id[meta_data$MG == "NO"])))
Met_tcga_norm_mean_MG_Island<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                              dplyr::select(meta_data$patient_id[meta_data$MG == "YES"]) %>%
                                                              filter(row.names(myNorm_met) %in% Island_met_loci)))
Met_tcga_norm_mean_MG_opensea<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                      dplyr::select(meta_data$patient_id[meta_data$MG == "YES"]) %>%
                                                      filter(row.names(myNorm_met) %in% OpenSea_met_loci)))
Met_tcga_norm_mean_MG_Shelf<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                              dplyr::select(meta_data$patient_id[meta_data$MG == "YES"]) %>%
                                                              filter(row.names(myNorm_met) %in% Shelf_met_loci)))
Met_tcga_norm_mean_MG_Shore<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                              dplyr::select(meta_data$patient_id[meta_data$MG == "YES"]) %>%
                                                              filter(row.names(myNorm_met) %in% Shore_met_loci)))
Met_tcga_norm_mean_NMG_Island<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                             dplyr::select(meta_data$patient_id[meta_data$MG == "NO"]) %>%
                                                            filter(row.names(myNorm_met) %in% Island_met_loci)))
Met_tcga_norm_mean_NMG_opensea<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                             dplyr::select(meta_data$patient_id[meta_data$MG == "NO"]) %>%
                                                              filter(row.names(myNorm_met) %in% OpenSea_met_loci)))
Met_tcga_norm_mean_NMG_Shelf<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                              dplyr::select(meta_data$patient_id[meta_data$MG == "NO"]) %>%
                                                              filter(row.names(myNorm_met) %in% Shelf_met_loci)))
Met_tcga_norm_mean_NMG_Shore<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
                                                              dplyr::select(meta_data$patient_id[meta_data$MG == "NO"]) %>%
                                                              filter(row.names(myNorm_met) %in% Shore_met_loci)))
p1<-ggplot(Met_tcga_norm_mean_MG_Island, aes(x = AveExpr))  + geom_density(color = "black", fill = "grey") + 
   ggtitle("Island MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p2<-ggplot(Met_tcga_norm_mean_NMG_Island, aes(x = AveExpr))  + geom_density(color = "black", fill = "grey") + 
   ggtitle("Island MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p3<-ggplot(Met_tcga_norm_mean_MG_Shore, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
   ggtitle("Shore MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p4<-ggplot(Met_tcga_norm_mean_NMG_Shore, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
   ggtitle("Shore MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p5<-ggplot(Met_tcga_norm_mean_MG_Shelf, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
   ggtitle("Shelf MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p6<-ggplot(Met_tcga_norm_mean_NMG_Shelf, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
   ggtitle("Shelf MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p7<-ggplot(Met_tcga_norm_mean_MG_opensea, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
   ggtitle("OpenSea MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
p8<-ggplot(Met_tcga_norm_mean_NMG_opensea, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
   ggtitle("OpenSea MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4)

## GSE队列 位于基因组不同位置（island, shore, opensea and shelf）的CpG位点的甲基化程度的密度图
# myNorm_met<-myNorm_met_gse
# meta_data_g<-meta_data_thymoma_kajiura
# Met_tcga_norm_mean_MG_Island<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                              dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "yes"]) %>%
#                                                              filter(row.names(myNorm_met) %in% Island_met_loci)))
# Met_tcga_norm_mean_MG_opensea<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                               dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "yes"]) %>%
#                                                               filter(row.names(myNorm_met) %in% OpenSea_met_loci)))
# Met_tcga_norm_mean_MG_Shelf<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                             dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "yes"]) %>%
#                                                             filter(row.names(myNorm_met) %in% Shelf_met_loci)))
# Met_tcga_norm_mean_MG_Shore<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                             dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "yes"]) %>%
#                                                             filter(row.names(myNorm_met) %in% Shore_met_loci)))
# Met_tcga_norm_mean_NMG_Island<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                               dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "no"]) %>%
#                                                               filter(row.names(myNorm_met) %in% Island_met_loci)))
# Met_tcga_norm_mean_NMG_opensea<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                                dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "no"]) %>%
#                                                                filter(row.names(myNorm_met) %in% OpenSea_met_loci)))
# Met_tcga_norm_mean_NMG_Shelf<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                              dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "no"]) %>%
#                                                              filter(row.names(myNorm_met) %in% Shelf_met_loci)))
# Met_tcga_norm_mean_NMG_Shore<-data.frame(AveExpr=rowMeans(data.table(myNorm_met) %>% 
#                                                              dplyr::select(meta_data$Sample[meta_data_g$Myasthenia_gravis == "no"]) %>%
#                                                              filter(row.names(myNorm_met) %in% Shore_met_loci)))
# p1<-ggplot(Met_tcga_norm_mean_MG_Island, aes(x = AveExpr))  + geom_density(color = "black", fill = "grey") + 
#    ggtitle("Island MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p2<-ggplot(Met_tcga_norm_mean_NMG_Island, aes(x = AveExpr))  + geom_density(color = "black", fill = "grey") + 
#    ggtitle("Island MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p3<-ggplot(Met_tcga_norm_mean_MG_Shore, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
#    ggtitle("Shore MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p4<-ggplot(Met_tcga_norm_mean_NMG_Shore, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
#    ggtitle("Shore MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p5<-ggplot(Met_tcga_norm_mean_MG_Shelf, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
#    ggtitle("Shelf MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p6<-ggplot(Met_tcga_norm_mean_NMG_Shelf, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
#    ggtitle("Shelf MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p7<-ggplot(Met_tcga_norm_mean_MG_opensea, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
#    ggtitle("OpenSea MG+")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# p8<-ggplot(Met_tcga_norm_mean_NMG_opensea, aes(x = AveExpr))  +  geom_density(color = "black", fill = "grey") + 
#    ggtitle("OpenSea MG-")  + theme_classic()+ theme(plot.title = element_text(vjust = -6,hjust=0.5,size=11))
# ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4)

# 
# #### 与某一个基因相关联的差异甲基化的CpG岛 significant TCGA cpg island of specific genes
 CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))
 CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "AIRE"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "GABRA5"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFL"))
CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "GRIK4"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "CACNA1A"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "CHRNA1"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "TTN"))
# CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR1"))
# me_da_isl_MG_mean[me_da_isl_MG$Islands_Name=="chr15:33602816-33604003",]
# me_da_isl_NO_mean[me_da_isl_NO$Islands_Name=="chr15:33602816-33604003",]
 CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))
 CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "GABRA5"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFL"))
CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))
CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR1"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "GRIK4"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "CACNA1A"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "CHRNA1"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "TTN"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR1"))
# 

## CPG 岛的甲基化上调下调的数量的柱状图 （TCGA） 
CPG_tcga_sig<-CPG_tcga[CPG_tcga$adj.P.Val<0.05,]
##write.csv(CPG_tcga_sig,"CPG_island_TCGA.csv")
cpg_tcga_down<-CPG_tcga[CPG_tcga$adj.P.Val<0.05 & CPG_tcga$logFC <0, ]
cpg_tcga_up<-CPG_tcga[CPG_tcga$adj.P.Val<0.05 & CPG_tcga$logFC >0, ]
DE_number<-data.frame(Group=c('Down','Up'),
                       DE_met_number=c(length(cpg_tcga_down$logFC),length(cpg_tcga_up$logFC)))
p <- ggplot(DE_number, aes(x=Group, y=DE_met_number,fill=Group)) +
   geom_bar(stat="identity",  colour="black", position=position_dodge(),width = 0.6)+
     theme_classic() + 
   theme(axis.line.x = element_line(colour = "black"),
         axis.line.y = element_line(colour = "black"))
p

# ## tcga 甲基化通路分析 gene pathway analysis
#  CPG_tcga_ge<-CPG_tcga_les0.5 %>% filter(adj.P.Val <0.05) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
#  CPG_tcga_ge_up<-CPG_tcga_les0.5 %>% filter(adj.P.Val <0.05 & logFC>0) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
#  CPG_tcga_ge_down<-CPG_tcga_les0.5 %>% filter(adj.P.Val <0.05 & logFC<0) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
#  
#  CPG_tcga_ge<-CPG_tcga_les0.5 %>% filter(adj.P.Val <0.05) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
#  CPG_tcga_ge_up<-CPG_tcga_les0.5 %>% filter(adj.P.Val <0.05 & logFC>0) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
#  CPG_tcga_ge_down<-CPG_tcga_les0.5 %>% filter(adj.P.Val <0.05 & logFC<0) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
# 
#  cpg_geid<-en2ENSE %>% filter(SYMBOL %in% unique(unlist(CPG_tcga_ge))) %>% pull(ENTREZID) %>% unique(.)
#  cpg_geid_up<-en2ENSE %>% filter(SYMBOL %in% unique(unlist(CPG_tcga_ge_up))) %>% pull(ENTREZID) %>% unique(.)
#  cpg_geid_down<-en2ENSE %>% filter(SYMBOL %in% unique(unlist(CPG_tcga_ge_down))) %>% pull(ENTREZID) %>% unique(.)
#  
 #### 自动运行时运行到这里就行了
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 cpg_tmg_cc_down <- enrichGO(gene          = cpg_geid_down,
                          universe      = result_tmg$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "none",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         minGSSize = 2,
                         maxGSSize = 150,
                         readable      = TRUE)
View(cpg_tmg_cc_down@result)
 neur_tcga_me_path<-cpg_tmg_cc_down
 neur_tcga_me_path@result<-cpg_tmg_cc_down@result %>% filter(str_detect(cpg_tmg_cc_down@result$Description,pattern = "neur") | 
                                                  str_detect(cpg_tmg_cc_down@result$Description,pattern = "glia"))
 ## only neural related pathways
 neur_tcga_me_path@result<-cpg_tmg_cc_down@result %>% filter(str_detect(cpg_tmg_cc_down@result$Description,pattern = "neur") | 
                                                               str_detect(cpg_tmg_cc_down@result$Description,pattern = "glia")|
                                                               str_detect(cpg_tmg_cc_down@result$Description,pattern = "nerv")|
                                                               str_detect(cpg_tmg_cc_down@result$Description,pattern = "axo")|
                                                               str_detect(cpg_tmg_cc_down@result$Description,pattern = "spin")|
                                                               str_detect(cpg_tmg_cc_down@result$Description,pattern = "skeletal muscle"))
 
 View(neur_tcga_me_path@result)
# write.csv(neur_tcga_me_path@result,"CPG_island_TCGA_down_neur_GO.csv")
 dotplot(neur_tcga_me_path,font.size=14)+
#   scale_color_material("#E6736F",reverse = T)
   scale_color_gradient(low="#436D5A",high="#8CD4B0")
 
barplot(neur_tcga_me_path)+
  scale_color_material("purple",reverse = T)
 View(cpg_tmg_cc_down@result)
# 
 cpg_tmg_cc_up <- enrichGO(gene          = cpg_geid_up,
                        universe      = result_tmg$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                       ont           = "ALL",
                        pAdjustMethod = "none",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1,
                        minGSSize = 2,
                       maxGSSize = 150,
                        readable      = TRUE)
 
 View(cpg_tmg_cc_up@result)
 dotplot(cpg_tmg_cc_up,showCategory = 37,font.size=14)+
   scale_color_gradient(low="#436D5A",high="#8CD4B0")
# write.csv(cpg_tmg_cc_up,"CPG_island_TCGA_up_GO.csv")

#########################差异甲基化的CpG岛和差异表达的mRNA基因关联分析
##选出只有一个关联基因的差异甲基化的CpG岛
 
CPG_tcga_single_gene_sig<-CPG_tcga %>% 
 filter(UCSC_RefGene_Name != "" & (!str_detect(CPG_tcga$UCSC_RefGene_Name,pattern =";"))) %>% 
   filter(adj.P.Val<0.05)
CPG_tcga_single_gene_sig<-CPG_tcga %>% filter(UCSC_RefGene_Name != "") %>% filter(adj.P.Val<0.05)
 length(CPG_tcga_single_gene_sig$Islands_Name)
res_tmg_t_sig<-res_tmg_t %>% filter(padj<0.05)
length(res_tmg_t_sig$baseMean)

## 筛选出关联基因的mRNA差异表达的差异甲基化的CPG islands
inter_cpg<-CPG_tcga_single_gene_sig %>% filter(UCSC_RefGene_Name %in% res_tmg_t_sig$gene_name)

## instersect gene names
inter<-intersect(CPG_tcga_single_gene_sig$UCSC_RefGene_Name,res_tmg_t_sig$gene_name)
cpg_inter<-CPG_tcga_single_gene_sig %>% filter(UCSC_RefGene_Name %in% inter) %>% dplyr::select(logFC,UCSC_RefGene_Name )
cpg_inter$logFC<-cpg_inter$logFC
mrna_sig<-res_tmg_t_sig %>% filter(gene_name %in% inter) %>% dplyr::select(log2FoldChange,gene_name )
inter_df<-cpg_inter %>% left_join(mrna_sig,by=c("UCSC_RefGene_Name"="gene_name"))
  #geom_smooth(aes(x=logFC,y=log2FoldChange),method = "lm")




### 后面先不管

if(0) {
## selected pathways
##ion channel pathway
channel_tcga_me_path<-go_up_tmg_cc
channel_tcga_me_path@result<-go_up_tmg_cc@result %>% filter(str_detect(go_up_tmg_cc@result$Description,pattern = "channel") |
                                                           str_detect(go_up_tmg_cc@result$Description,pattern = "transport")) %>%  arrange(pvalue)
dotplot(channel_tcga_me_path, showCategory = 15,font.size=14) +
  scale_color_material("red",reverse = T)

write.csv(channel_tcga_me_path@result,"mRNA_cpg_island_co_up_channel.csv")

## neuro related pathway
neuro_tcga_me_path<-go_up_tmg_cc
neuro_tcga_me_path@result<-go_up_tmg_cc@result %>%
  filter(str_detect(go_up_tmg_cc@result$Description,pattern = "syna") |
                                                              str_detect(go_up_tmg_cc@result$Description,pattern = "neu")|
                                                              str_detect(go_up_tmg_cc@result$Description,pattern = "nerv")|
                                                            str_detect(go_up_tmg_cc@result$Description,pattern = "sarc")) %>%
  filter(Count>1 | Description=="postsynaptic cytoskeleton" )


dotplot(neuro_tcga_me_path, showCategory = 15,font.size=14) +
  scale_color_material("indigo",reverse = T)
write.csv(neuro_tcga_me_path@result,"mRNA_cpg_island_co_up_neuro.csv")
}

##纳入的神经相关的基因关联的CpG island通路的热图
# ### heatmap of inlcuded genes
# 
# library(ComplexHeatmap)
# 
# 
# abc<-neuro_tcga_me_path2@result %>% arrange(geneID)
# 
# neuro_list<-unique(unlist(str_split(abc$geneID,pattern = "/")))
# 
# 
# axon_t_df<-data.frame(geneID=neuro_list)
# 
# mt<-matrix(data=NA,nrow=length(neuro_list),
#            ncol=length(abc$Description)) %>% data.frame()
# row.names(mt)<-neuro_list
# colnames(mt)<-neuro_tcga_me_path2@result$Description
# 
# for (i in c(1:length(neuro_tcga_me_path2@result$Description))) {
#   for (j in c(1:length(neuro_list)))
#   { 
#     if ( neuro_list[j] %in% str_split(neuro_tcga_me_path2@result$geneID[i], pattern = "/")[[1]]) 
#     {mt[j,i]<-1}
#     else {mt[j,i]<-0}
#   }
# }
# 
# axon_t_df<- cbind(axon_t_df,mt) %>% data.frame() %>% arrange()
#   dplyr::mutate(rowSums = rowSums(.[,4:length(colnames(.))])) %>% arrange(desc(rowSums)) %>% dplyr::select(-rowSums)
# row.names(axon_t_df)<-axon_t_df$geneID
# 
# 
# axon_t_df_t<-t(axon_t_df[,2:length(colnames(axon_t_df))]) %>% data.frame() %>% arrange(GABRA5,HSPA2,RYR3,
#                                                                                        VWC2,PTPRS,SPOCK1,
#                                                                                        MCOLN2,IFT172,CCDC88C,NEFM)
# 
# 
# 
# Heatmap(axon_t_df_t,cluster_rows = F,cluster_columns = F,
#         col=(c("white",2:10)), border = T,rect_gp= gpar(col = "black", lwd = 1),
#         row_names_rot = 0,
#         column_names_rot=30,
#         column_names_side="bottom",
#         row_names_side="left",
#         row_names_gp = gpar(fontsize = 10),
#         column_names_gp = gpar(fontsize = 6))
# 
# 
# 
# 
# ha2 = columnAnnotation(MG_vs_No_MG = anno_barplot((as.matrix(axon_t_df[,c(2)])),gp = gpar(fill = 3)))
# 
# 
# Heatmap(t(axon_t_df[,4:length(colnames(axon_t_df))]),cluster_rows = F,cluster_columns = F,
#         col=(c("white",2:10)), border = T,rect_gp= gpar(col = "black", lwd = 1),   
#         top_annotation = ha2,
#         row_names_rot = 0,
#         column_names_rot=60,
#         column_names_side="bottom",
#         row_names_side="left",
#         row_names_gp = gpar(fontsize = 8),
#         column_names_gp = gpar(fontsize = 7))
# 
# 
# 
# 
if(0) {
go_do_tmg_cc <- enrichGO(gene          = cpg_exp_geid_down,
                         universe      = result_tmg$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "none",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
barplot(go_do_tmg_cc,showCategory = 45)
dotplot(go_do_tmg_cc,showCategory = 12)

neuro_tcga_me_path<-go_do_tmg_cc
neuro_tcga_me_path@result<-go_do_tmg_cc@result %>%
  filter(str_detect(go_do_tmg_cc@result$Description,pattern = "syna") |
           str_detect(go_do_tmg_cc@result$Description,pattern = "neu")|
           str_detect(go_do_tmg_cc@result$Description,pattern = "nerv")|
           str_detect(go_do_tmg_cc@result$Description,pattern = "sarc")) %>%
  filter(Count>1 | Description=="postsynaptic cytoskeleton" )

View(go_do_tmg_cc@result)
dotplot(neuro_tcga_me_path, showCategory = 15,font.size=14) +
  scale_color_material("indigo",reverse = T)

cpg_exp_ge_do<-intersect(unique(unlist(CPG_tcga_ge)),result_tmg_down$gene_name)



## GSE dataset
me_da_islg<-isl_info %>% right_join(me_da_g, by="me_id") %>% 
   filter(Islands_Name!="" & Relation_to_Island =="Island")

me_num_pislg<-me_da_islg %>% group_by(Islands_Name) %>% summarise(count = length(Islands_Name)) %>% 
   filter(count>=3) %>% arrange(Islands_Name)
 me_num_pislg
 
 me_da_islg<-me_da_islg %>% filter(Islands_Name %in% me_num_pislg$Islands_Name) %>% arrange(Islands_Name)
 me_da_islg$Islands_Name
}
# 
# 
# 
# me_da_islg[me_da_islg$Islands_Name=="chr8:24770908-24772547",]
# 
# me_num_pislg[me_num_pislg$Islands_Name=="chr8:24770908-24772547",]
# 
# 
# 
# 
# 
# ## GSE
# 
# Sum_table<-data.frame()
# j<-1
# len<-length(me_num_pislg$Islands_Name)
# 
# for (i in c(1:len)) {
#   me_da_t<-me_da_islg[(j:((me_num_pislg$count[i])+j-1)),]
#   
#   me_da_t2<-me_da_t %>% data.table::melt(colnames(.)[1:10],colnames(.)[11:length(colnames(.))],
#                              variable.name="patient_id",value.name="me_beta") %>% arrange(pos) %>% 
#     left_join(meta_data_thymoma_kajiura, by=c("patient_id"="Sample") )%>% 
#     mutate(me_isl_id=paste0(me_id,"-",patient_id))
#   
#   expr_isl_t<-me_da_t2 %>% dplyr::select(me_isl_id,me_beta) %>% 
#     spread(key=me_isl_id,value=me_beta) %>% 
#     dplyr::select(me_da_t2$me_isl_id)
#   
#   met_loci<-factor(me_da_t2$me_id)
#   MG<-factor(me_da_t2$Myasthenia_gravis,levels=c("yes","no"))
#   design<-model.matrix(~met_loci+MG)
#   
#   fit <- lmFit(expr_isl_t, design)
#   fit <- eBayes(fit)
#   topTable<-topTable(fit, coef="MGno") %>% 
#     mutate(chr=me_da_t$chr[1], 
#            Islands_Name= me_da_t$Islands_Name[1], 
#            me_id = paste(me_da_t$me_id,collapse= ";"), 
#            UCSC_RefGene_Name= paste(me_da_t$UCSC_RefGene_Name,collapse= ";"),
#            me_count = me_num_pislg$count[i])
#   
#   logFC<-log2(mean(me_da_t2$me_beta[me_da_t2$Myasthenia_gravis=="yes"])/mean(me_da_t2$me_beta[me_da_t2$Myasthenia_gravis=="no"]))
#   topTable$logFC<-logFC
#   
#   Sum_table<-rbind(Sum_table,topTable)
#   j<-j+me_num_pislg$count[i]
# }
# 
# Sum_table_gse<-Sum_table %>% dplyr::select(Islands_Name,chr,me_id,me_count,UCSC_RefGene_Name,
#                                     everything())
# write.csv(Sum_table_gse,"CPG_island_differential_linear_model_gse_new.csv")
# 
 CPG_gse<-read.csv("CPG_island_differential_linear_model_gse_new.csv",header=T)
# 
# 
 CPG_gse$adj.P.Val<-p.adjust(CPG_gse$P.Value,method = "fdr")
# 
 CPG_gse <- CPG_gse %>% dplyr::arrange(P.Value)
# 
 CPG_gene_list<-CPG_gse$UCSC_RefGene_Name %>% str_split(pattern = ";")
 CPG_gene_list_c<-vector()
 for (i in c(1:length(CPG_gene_list))){
   CPG_gene_t<-unique(CPG_gene_list[[i]])
   CPG_gene_list_c[i]<-str_c(CPG_gene_t,collapse = ";")
 }
 
 CPG_gse$UCSC_RefGene_Name <-CPG_gene_list_c %>% str_remove(pattern="^;") %>% str_remove(pattern=";$")
# CPG_gse
# 
# 
 CPG_gse_sig<-CPG_gse[CPG_gse$adj.P.Val<0.05, ]
# write.csv(CPG_gse_sig,"CPG_island_gse.csv")
# 
# plot(CPG_gse$logFC)
# plot(CPG_tcga$logFC)
# 
# ## scatter plot of gse and tcga CPG island fold change for tcga mRNA and CPG island co-significant genes
# 
 inter_gse_cpg<-CPG_gse %>% filter(Islands_Name %in% inter_cpg$Islands_Name)
 inter_gse_cpg$logFC<-inter_gse_cpg$logFC
# 
 cpg_gse_inter_fc<-inter_cpg %>% dplyr::select(Islands_Name,tcga_logFC=logFC,UCSC_RefGene_Name) %>% 
   left_join(inter_gse_cpg[,c("Islands_Name","logFC")],by="Islands_Name") %>% dplyr::rename(gse_logFC=logFC)
# 
# write.csv(cpg_gse_inter_fc,"cpg_gse_tcga_inter_fc.csv")
# 
 unique(cpg_gse_inter_fc$UCSC_RefGene_Name)
 cpg_gse_inter_fc[cpg_gse_inter_fc$tcga_logFC >0 & cpg_gse_inter_fc$gse_logFC >0,]
 cpg_gse_inter_fc[cpg_gse_inter_fc$tcga_logFC <0 & cpg_gse_inter_fc$gse_logFC <0,]
# 
# 
 a<-inter_df$UCSC_RefGene_Name[inter_df$logFC>0 & inter_df$log2FoldChange<0]
 b<-inter_df$UCSC_RefGene_Name[inter_df$logFC<0 & inter_df$log2FoldChange>0]
# 
 c<-cpg_gse_inter_fc %>% filter(UCSC_RefGene_Name %in% b) %>% filter(gse_logFC<0) %>%  pull(UCSC_RefGene_Name)
# 
 d<-inter_df$UCSC_RefGene_Name[inter_df$logFC>0 & inter_df$log2FoldChange>0]
 e<-inter_df$UCSC_RefGene_Name[inter_df$logFC<0 & inter_df$log2FoldChange<0]
 f<-cpg_gse_inter_fc %>% filter(UCSC_RefGene_Name %in% d) %>% filter(gse_logFC>0) %>%  pull(UCSC_RefGene_Name)
 g<-cpg_gse_inter_fc %>% filter(UCSC_RefGene_Name %in% e) %>% filter(gse_logFC<0) %>%  pull(UCSC_RefGene_Name)
# 
# 
# 
# 
# 
 library(plyr)
 ggplot(data = cpg_gse_inter_fc) + #
   geom_point(aes(x=gse_logFC,y=tcga_logFC),color="black") + 
   geom_hline(yintercept = 0) + 
   geom_vline(xintercept = 0) + 
   theme_classic(base_size=16) + 
   geom_text_repel(data=subset(cpg_gse_inter_fc, UCSC_RefGene_Name %in% c("NEFM","RYR3","GABRA5")), 
                   aes(x=gse_logFC,y=tcga_logFC,label=UCSC_RefGene_Name))+  
   geom_point(data=subset(cpg_gse_inter_fc, UCSC_RefGene_Name %in% c(a,c)),
              color="red",aes(x=gse_logFC,y=tcga_logFC),size=1) + 
   geom_smooth(aes(x=gse_logFC,y=tcga_logFC),color="black", linetype="dashed",method = "lm") +
   geom_point(data=subset(cpg_gse_inter_fc, UCSC_RefGene_Name %in% c(f,g)),
              color="blue",aes(x=gse_logFC,y=tcga_logFC),size=1)+ 
   geom_point(data=subset(cpg_gse_inter_fc, UCSC_RefGene_Name %in% c("NEFM","RYR3","GABRA5")),
              color="red",aes(x=gse_logFC,y=tcga_logFC),size=3) 
# 
# 
# ?geom_smooth
# 
# ## scatter plot of gse and tcga CPG island fold change for tcga CPG island significant genes
 CPG_tcga_single_gene_sig_add_gse<-CPG_tcga_single_gene_sig %>% dplyr::select(Islands_Name,tcga_logFC=logFC,UCSC_RefGene_Name )%>%
   left_join(CPG_gse[,c("Islands_Name","logFC")],by="Islands_Name") %>% dplyr::rename(gse_logFC=logFC)
# 
 ggplot(data = CPG_tcga_single_gene_sig_add_gse) + 
   geom_point(aes(x=gse_logFC,y=tcga_logFC),color="black",size=0.9) + 
   geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
   theme_classic(base_size=16) + 
   geom_smooth(aes(x=gse_logFC,y=tcga_logFC),method = "lm") 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ############# GSE and TCGA CPG co sig analysis
# 
# CPG_tcga_single_gene_sig  ## tcga isl sig data
# 
# CPG_gse_single_gene_sig<-CPG_gse %>% filter(UCSC_RefGene_Name != "" & (!str_detect(CPG_gse$UCSC_RefGene_Name,pattern =";"))) %>% filter(adj.P.Val<0.1) ## gse isl sig data
# 
# ## gse and tcga co sig island name
# CPG_gse_tcga_inter_list<-intersect(CPG_tcga_single_gene_sig$Islands_Name,CPG_gse_single_gene_sig$Islands_Name)
# ## gse and tcga co sig related island name
# CPG_gse_tcga_inter_list_gene_id<-CPG_gse_single_gene_sig %>% filter(Islands_Name %in% CPG_gse_tcga_inter_list) %>% pull(UCSC_RefGene_Name)
#    ## gse and tcga co sig related island related gene pathway analysis 
# 
# 
# CPG_gse_tcga_inter_gene_id_ent<-en2ENSE %>% filter(SYMBOL %in% CPG_gse_tcga_inter_list_gene_id) %>% pull(ENTREZID)
# cpg_exp_geid_down<-en2ENSE %>% filter(SYMBOL %in% cpg_exp_ge_down) %>% pull(ENTREZID)
# 
# 
# 
# CPG_gse_tcga_inter_gene_id_down<-intersect(CPG_tcga_single_gene_sig$UCSC_RefGene_Name[CPG_tcga_single_gene_sig$logFC>0],
#                                                CPG_gse_single_gene_sig$UCSC_RefGene_Name[CPG_gse_single_gene_sig$logFC>0])
# CPG_gse_tcga_inter_gene_id_ent_down<-en2ENSE %>% filter(SYMBOL %in% CPG_gse_tcga_inter_gene_id_down) %>% pull(ENTREZID)
# 
# CPG_gse_tcga_inter_gene_id_up<-intersect(CPG_tcga_single_gene_sig$UCSC_RefGene_Name[CPG_tcga_single_gene_sig$logFC<0],
#                                              CPG_gse_single_gene_sig$UCSC_RefGene_Name[CPG_gse_single_gene_sig$logFC<0])
# 
# CPG_gse_tcga_inter_gene_id_ent_up<-en2ENSE %>% filter(SYMBOL %in% CPG_gse_tcga_inter_gene_id_up) %>% pull(ENTREZID)
# 
# 
# go_tmg <- enrichGO(gene          = CPG_gse_tcga_inter_gene_id_ent,
#                          OrgDb         = org.Hs.eg.db,
#                          ont           = "ALL",
#                          pAdjustMethod = "none",
#                          minGSSize = 3,
#                          pvalueCutoff  = 0.05,
#                          readable      = TRUE)
# go_tmg@result %>% filter(go_tmg@result$Count>1)
# 
# go_tmg_down <- enrichGO(gene          = CPG_gse_tcga_inter_gene_id_ent_down,
#                    OrgDb         = org.Hs.eg.db,
#                    ont           = "ALL",
#                    pAdjustMethod = "none",
#                    minGSSize = 1,
#                    pvalueCutoff  = 0.05,
#                    readable      = TRUE)
# go_tmg_down@result
# 
# 
# ## selected pathways
# 
# 
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "GABRA5"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFL"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))
# CPG_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "GRIK4"))
# 
# 
# 
# 
# 
# ## gse gene pathway analysis
# 
# CPG_gse_ge<-CPG_gse %>% filter(P.Value <0.0005) %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
# 
# cpg_gse_geid<-en2ENSE %>% filter(SYMBOL %in% unique(unlist(CPG_gse_ge))) %>% pull(ENTREZID) %>% unique(.)
# 
# 
# cpg_tmg_cc_gse <- enrichGO(gene          = cpg_gse_geid,
#                        universe      = result_tmg$ENTREZID,
#                        OrgDb         = org.Hs.eg.db,
#                        ont           = "CC",
#                        pAdjustMethod = "none",
#                        pvalueCutoff  = 0.05,
#                        readable      = TRUE)
# 
# barplot(cpg_tmg_cc_gse,showCategory = 37)
# ##
# 
# 
# Sum_table_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))
# Sum_table_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))
# 
# 
# Sum_table_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))
# Sum_table_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))
# 
# 
# 
# 
# 
# 
# ####
# 
# ################################################################################################################
# 
# ## TCGA CPG island tSNE plot 
# ## TCGA CPG island tSNE analysis
# library(Rtsne)
# library(lumi)
# 
# me_da_isl_mean_beta
# 
# myNorm_M<-beta2m(me_da_isl_mean_beta[,2:122])
# rownames(myNorm_M)<-me_da_isl_mean_beta$Islands_Name
# 
# View(myNorm_M)
# 
# myNorm_M_v<-rowVars(as.matrix(myNorm_M))
# myNorm_M_t<-myNorm_M %>%  as.data.frame() %>% arrange(desc(myNorm_M_v))
# 
# 
# 
# set.seed(321) # 设置随机数种子
# dim(myNorm_M_t)
# tsne_out = Rtsne(
#   t(myNorm_M_t[1:1000,]),
#   dims = 2,
#   pca = T,
#   max_iter = 5000,
#   theta = 0.1,
#   perplexity = 10,
#   verbose = F
# ) 
# tsne_result = as.data.frame(tsne_out$Y)
# colnames(tsne_result) = c("tSNE1","tSNE2")
# row.names(tsne_result)<-met_meta$patient_id
# 
# ggplot(tsne_result,aes(tSNE1,tSNE2,color=met_meta$MG)) +
#   geom_point(size=2) +
#   stat_ellipse(level = 0.5) +
#   scale_color_npg() + 
#   theme_classic()
# 
# 
# ggplot(tsne_result,aes(tSNE1,tSNE2,color=met_meta$histo_type)) + 
#   geom_point(size=2) +
#   stat_ellipse(level = 0.5) +
#   scale_color_npg() + 
#   theme_classic()
# 
# 
# 
# ###
# 
# 
# 
# 