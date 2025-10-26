library(stringr)
library("ggpmisc")

## import the CPG island 差异甲基化分析结果
CPG_tcga<-read.csv("CPG_island_differential_linear_model_TCGA_new.csv",header=T)

#Calcualte the RPKM value of all expressed genes, and calculate the mean RPKM of NSCs
ensembl_list <-ExpDataCPM %>%
  mutate(ensemblID = gsub('[.][0-9]+$','',.$ensemblID)) %>% pull(1)

#ensembl_length<-getGeneLengthAndGCContent(ensembl_list, "hsa")
#write.csv(ensembl_length,'ensemble_gene_length.csv')
ensembl_length<-read.csv('D:/Documents/TMG_methylation_paper/reference_data/ensemble_gene_length.csv')

cpmMatrixFiltered_unlog <- data.frame(gene_id=as.character(RNA_expr_THYM$gene_id), Count2CPM(RNA_expr_THYM[,-1]),check.names = F)
cpm_matrix_Ense<-left_join(ensembl_length[,c(1,2)],cpmMatrixFiltered_unlog,by=c('X'='gene_id'))

## FPKM
#change count of fragments to RFKM: RFKM ( log transformed)= CPM*1000/gene_length  + 1
ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,
                         log2(cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2]+1),check.names = F)
#FPKM not log transformed
#ExpData_FPKM<-data.frame(ensemblID=cpm_matrix_Ense$X,
#                       cpm_matrix_Ense[,c(3:(length(colnames(cpm_matrix_Ense))-1))]*1000/cpm_matrix_Ense[,2],check.names = F)

k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,
                               columns=c("ENTREZID",'SYMBOL'),
                               keytype = 'ENSEMBL')
ExpData_FPKM <- ExpData_FPKM %>% left_join(en2ENSE[,c(1,3)],by=c("ensemblID"="ENSEMBL"))


### 单个基因：CpG岛甲基化水平与关联基因的mRNA的关系

### NEFM
 # CpG平均值
me_da_isl_median_beta

CpG_mean_NEFM<-me_da_isl[me_da_isl$Islands_Name =="chr8:24770908-24772547",] %>% 
  dplyr::select(c(11:length(colnames(me_da_isl)))) %>% apply(.,2,median)
CpG_mean_NEFM<-data.frame(CpG=CpG_mean_NEFM,ID=names(CpG_mean_NEFM)) 


 # mRNA水平
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "NEFM") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
#FPKM<-ExpData_FPKM %>% filter(SYMBOL == "DNMT1") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
FPKM$ID<-row.names(FPKM)
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"

 ## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_NEFM,by="ID") %>% na.omit() %>% left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))


p1<-ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("NEFM","-","chr8:24770908-24772547"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)


## 单个基因：CpG岛甲基化水平与关联基因的mRNA散点图，并且添加趋势线
ggplot(CpG_mRNA, aes(CpG,mRNA,color=MG)) + 
  geom_point(aes(color=MG),size =1) +  
  geom_smooth(aes(color = MG),method = "lm",se = F, show.legend = T,size=0.5)+  
  ggtitle("NEFM")+ 
  theme_classic()+ 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_lancet()


CpG_mRNA_co<-CpG_mRNA %>% filter(MG=="NO")
summary(lm(CpG_mRNA_co$mRNA~CpG_mRNA_co$CpG))
CpG_mRNA_ca<-CpG_mRNA %>% filter(MG=="YES")
summary(lm(CpG_mRNA_ca$mRNA~CpG_mRNA_ca$CpG))
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))


#chr21:45705428-45706044
#chr21:45713509-45713813
View(me_da_isl[me_da_isl$Islands_Name =="chr8:24770908-24772547",])
CpG_mean_NEFM<-CpG_locus_mRNA("NEFM","cg09234518")
CpG_mRNA<-CpG_mean_NEFM

ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("NEFM"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)

#write.csv(CpG_mRNA_co,"./data_middle/CpG_mRNA_co.csv")
#write.csv(CpG_mRNA_ca,"./data_middle/CpG_mRNA_ca.csv")

## RYR3
me_da_isl[str_detect(me_da_isl$UCSC_RefGene_Name,pattern ="RYR3"),]
CpG_mean_RYR3<-me_da_isl[me_da_isl$Islands_Name =="chr15:33602816-33604003",] %>% 
  dplyr::select(c(11:length(colnames(me_da_isl)))) %>% apply(.,2,median)
CpG_mean_RYR3<-data.frame(CpG=CpG_mean_RYR3,ID=names(CpG_mean_RYR3))
# mRNA水平
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "RYR3") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame() 
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "DNMT1") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame() 
FPKM$ID<-row.names(FPKM) 
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"
## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_RYR3,by="ID") %>% na.omit() %>% left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))
## 画图
ggplot(CpG_mRNA, aes(CpG,mRNA,color=MG)) + 
  geom_point(aes(color=MG),size =1) +  
  geom_smooth(aes(color = MG),method = "lm",se = F, show.legend = T,size=0.5)+  
  ggtitle("RYR3")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_color_lancet() 


CpG_mRNA_co<-CpG_mRNA %>% filter(MG=="NO")
summary(lm(CpG_mRNA_co$mRNA~CpG_mRNA_co$CpG))
CpG_mRNA_ca<-CpG_mRNA %>% filter(MG=="YES")
summary(lm(CpG_mRNA_ca$mRNA~CpG_mRNA_ca$CpG))
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

## RYR3单个位点
View(me_da_isl[me_da_isl$Islands_Name =="chr15:33602816-33604003",])
CpG_mean_NEFM<-CpG_locus_mRNA("RYR3","cg25405123")
CpG_mRNA<-CpG_mean_NEFM

ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("NEFM"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)

#write.csv(CpG_mRNA_co,"./data_middle/CpG_mRNA_RYR3_co.csv")
#write.csv(CpG_mRNA_ca,"./data_middle/CpG_mRNA_RYR3_ca.csv")

## GABARB
CpG_mean_GABRA5<-me_da_isl[str_detect(me_da_isl$UCSC_RefGene_Name,pattern ="GABRA5"),]
CpG_mean_GABRA5<-me_da_isl[me_da_isl$Islands_Name =="chr15:27112030-27113479",] %>% 
  dplyr::select(c(11:length(colnames(me_da_isl)))) %>% apply(.,2,mean)
CpG_mean_GABRA5<-data.frame(CpG=CpG_mean_GABRA5,ID=names(CpG_mean_GABRA5))   
# mRNA水平
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "GABRA5") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "DNMT1") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
FPKM$ID<-row.names(FPKM)
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"
## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_GABRA5,by="ID") %>% na.omit() %>% left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))
## 画图
ggplot(CpG_mRNA, aes(CpG,mRNA,color=MG)) + 
  geom_point(aes(color=MG),size =1) +  
  geom_smooth(aes(color = MG),method = "lm",se = F, show.legend = T,size=0.5)+  
  ggtitle("GABRA5")+ 
  theme_classic()+ 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_color_lancet() 

p1<-ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("GABRA5","-",""))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)


CpG_mRNA_co<-CpG_mRNA %>% filter(MG=="NO") 
summary(lm(CpG_mRNA_co$mRNA~CpG_mRNA_co$CpG)) 
CpG_mRNA_ca<-CpG_mRNA %>% filter(MG=="YES")
summary(lm(CpG_mRNA_ca$mRNA~CpG_mRNA_ca$CpG))

#write.csv(CpG_mRNA_co,"./data_middle/CpG_mRNA_GABRA5_co.csv")
#write.csv(CpG_mRNA_ca,"./data_middle/CpG_mRNA_GABRA5_ca.csv")


## CHRNA1 基因与单个CpG loci 关系，一共七个CpG loci
me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),]

 ## tss200
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg21552300")
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg18003659")
 ## body
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg13910973")
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg20719675")
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg22925639")
 ## ts1500
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg00748431")
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg05649009")
CpG_mean_CHRNA1<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="CHRNA1;"),] %>% filter(me_id=="cg24719366")

## 写一个funciton，输出为correlation plot 和相关性检验
CpG_locus_mRNA<-function(gene_sym,cp_id) {
CpG_mean_gene <- me_da_anno %>% filter(me_id==cp_id) %>% 
  dplyr::select(c(11:length(colnames(me_da_isl)))) %>% colMeans() 
CpG_mean_gene<-data.frame(CpG=CpG_mean_gene,ID=names(CpG_mean_gene)) 
   # mRNA水平
FPKM<-ExpData_FPKM %>% filter(SYMBOL == gene_sym) %>% 
  dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
FPKM$ID<-row.names(FPKM)
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"
   ## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_gene,by="ID") %>% na.omit() %>% 
  left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))
   ## 画图
return(CpG_mRNA)
}


summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))
CpG_mRNA_coca<-CpG_mRNA%>% filter(MG=="YES") 
summary(lm(CpG_mRNA_coca$mRNA~CpG_mRNA_coca$CpG))
CpG_mRNA_coca<-CpG_mRNA%>% filter(MG=="NO") 
summary(lm(CpG_mRNA_coca$mRNA~CpG_mRNA_coca$CpG))

correlation_plot_cpg_locus<-function(gene_sym,cp_id) {
  CpG_mean_gene <- me_da_anno %>% filter(me_id==cp_id) %>% 
    dplyr::select(c(11:length(colnames(me_da_isl)))) %>% colMeans() 
  CpG_mean_gene<-data.frame(CpG=CpG_mean_gene,ID=names(CpG_mean_gene)) 
  # mRNA水平
  FPKM<-ExpData_FPKM %>% filter(SYMBOL == gene_sym) %>% 
    dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
  FPKM$ID<-row.names(FPKM)
  colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"
  ## mRNA + CpG
  CpG_mRNA<-FPKM %>% left_join(CpG_mean_gene,by="ID") %>% na.omit() %>% 
    left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))
  ## 画图
  p1<-ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
    geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
    geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
    geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
    ggtitle(paste0(gene_sym,"-",cp_id))+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = y ~ x),
                    mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                    parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
    stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)
  p2<-summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))
  return(p1)
}


CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg18003659")
a1<-correlation_plot_cpg_locus("CHRNA1","cg18003659")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg21552300")
a2<-correlation_plot_cpg_locus("CHRNA1","cg21552300")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg13910973")
a3<-correlation_plot_cpg_locus("CHRNA1","cg13910973")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg20719675")
a4<-correlation_plot_cpg_locus("CHRNA1","cg20719675")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg22925639")
a5<-correlation_plot_cpg_locus("CHRNA1","cg22925639")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg00748431")
a6<-correlation_plot_cpg_locus("CHRNA1","cg00748431")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg05649009")
a7<-correlation_plot_cpg_locus("CHRNA1","cg05649009")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

CpG_mRNA<-CpG_locus_mRNA("CHRNA1","cg24719366")
a8<-correlation_plot_cpg_locus("CHRNA1","cg24719366")
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

library(gridExtra)
grid.arrange(a1,a2,a3,a4,a5,a6,a7,a8,ncol=2)

## AIRE
  AIRE_cp_ID<-na.omit(me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="AIRE"),])
View(AIRE_cp_ID)

#chr21:45705428-45706044
#chr21:45713509-45713813

CpG_mean_AIRE1<-me_da_isl[me_da_isl$Islands_Name =="chr21:45713509-45713813",] %>% 
  dplyr::select(c(11:length(colnames(me_da_isl)))) %>% apply(.,2,median)
View(me_da_isl[me_da_isl$Islands_Name =="chr21:45705428-45706044",])
CpG_mean_AIRE1<-data.frame(CpG=CpG_mean_AIRE1,ID=names(CpG_mean_AIRE1))


# mRNA水平
FPKM<-ExpDataCPM %>% filter(ensemblID == "ENSG00000160224") %>% dplyr::select(-c("ensemblID","GeneSymbol")) %>% t() %>% as.data.frame() 
FPKM$ID<-row.names(FPKM) 
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"
## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_AIRE1,by="ID") %>% na.omit() %>% left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))
## 画图
ggplot(CpG_mRNA, aes(CpG,mRNA,color=MG)) + 
  geom_point(aes(color=MG),size =1) +  
  geom_smooth(aes(color = MG),method = "lm",se = F, show.legend = T,size=0.5)+  
  ggtitle("AIRE1")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_color_lancet() 

p1<-ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("AIRE","-","chr21:45705428-45706044"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)

p2<-ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("AIRE","-","chr21:45713509-45713813"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)

library(gridExtra)
grid.arrange(p1,p2,ncol=2)

##AIRE单个甲基化位点与mRNA关系

#chr21:45705428-45706044
#chr21:45713509-45713813
View(me_da_isl[me_da_isl$Islands_Name =="chr21:45705428-45706044",])
CpG_mean_AIRE1<-CpG_locus_mRNA("AIRE","cg09510531")
CpG_mRNA<-FPKM %>% left_join(CpG_mean_AIRE1,by="ID") %>% na.omit() 

ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("AIRE"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)

p2<-summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))
return(p1)
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG))

## TTN
CpG_mean_TTN<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="TTN") & str_detect(me_da_anno$chr,pattern ="chr2"),] 
View(CpG_mean_TTN)
#View(CpG_mean_TTN)
## cg05185019,cg10859358
CpG_mean_TTN<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="TTN") & str_detect(me_da_anno$chr,pattern ="chr2"),] %>% filter(me_id=="cg05185019")
CpG_mean_TTN<-me_da_anno[str_detect(me_da_anno$UCSC_RefGene_Name,pattern ="TTN") & str_detect(me_da_anno$chr,pattern ="chr2"),] %>% filter(me_id=="cg10859358")

CpG_mean_TTN<-CpG_mean_TTN %>% dplyr::select(c(11:length(colnames(me_da_isl)))) %>% colMeans() 
CpG_mean_TTN<-data.frame(CpG=CpG_mean_TTN,ID=names(CpG_mean_TTN)) 

# mRNA水平
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "TTN") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame()
FPKM$ID<-row.names(FPKM)
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"
## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_TTN,by="ID") %>% na.omit() %>% left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))

## 画图
ggplot(CpG_mRNA, aes(CpG,mRNA,color=MG)) + 
  geom_point(aes(color=MG),size =1) +  
  geom_smooth(aes(color = MG),method = "lm",se = F, show.legend = T,size=0.5)+  
  ggtitle("TTN-cg05185019") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_lancet() 


### RYR1基因
# CpG平均值
CpG_mean_RYR<-me_da_isl[str_detect(me_da_isl$UCSC_RefGene_Name,pattern ="RYR1"),]
View(CpG_mean_RYR)

CpG_mean_RYR1<-me_da_isl[me_da_isl$Islands_Name =="chr19:38924239-38924519",] %>% 
  dplyr::select(c(11:length(colnames(me_da_isl)))) %>% colMeans()
CpG_mean_RYR1<-data.frame(CpG=CpG_mean_RYR1,ID=names(CpG_mean_RYR1))
# mRNA水平
FPKM<-ExpData_FPKM %>% filter(SYMBOL == "RYR1") %>% dplyr::select(-c("ensemblID","SYMBOL")) %>% t() %>% as.data.frame() 
FPKM$ID<-row.names(FPKM)
colnames(FPKM)[colnames(FPKM)=="V1"]<-"mRNA"

## mRNA + CpG
CpG_mRNA<-FPKM %>% left_join(CpG_mean_RYR1,by="ID") %>% na.omit() %>% left_join(meta_data[,c("patient_id","MG")],by=c("ID"="patient_id"))

## 画图
ggplot(CpG_mRNA, aes(CpG,mRNA,color=MG)) +  
  geom_point(aes(color=MG),size =1) +  
  geom_smooth(aes(color = MG),method = "lm",se = F, show.legend = T,size=0.5)+  
  ggtitle("RYR1")+ 
  theme_classic()+ 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_color_lancet() 
summary(lm(CpG_mRNA$mRNA~CpG_mRNA$CpG)) 

ggplot(CpG_mRNA, aes(CpG,mRNA)) + 
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="black")+  
  geom_point(aes(color=MG),size =1) +  scale_color_manual(values=c("#7DBE9D","#E6736F"))+
  geom_smooth(CpG_mRNA, mapping=aes(x=CpG,y=mRNA,color=MG),method = "lm",se = F, show.legend = T,size=0.5,linetype="dashed")+
  ggtitle(paste0("RYR1-chr19:38924239-38924519"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)


### 所有基因的甲基化与mRNA表达的关联
CPG_tcga<-read.csv("CPG_island_differential_linear_model_TCGA_new.csv",header=T)
CPG_tcga_ge_id<-CPG_tcga %>% dplyr::select(UCSC_RefGene_Name) %>% pull() %>% str_split(pattern = ";")
CPG_tcga_ge_id<-lapply(CPG_tcga_ge_id,unique)
CPG_tcga_ge_id_i<-CPG_tcga_ge_id

## 将list 转换为 data.frame
CPG_tcga_ge<-as.data.frame(t(sapply(CPG_tcga_ge_id,"[",i=1:max(sapply(CPG_tcga_ge_id_i,length)))))
CPG_tcga_ge<-apply(CPG_tcga_ge,2,function(x) ifelse(x=="",NA,x))

for (i in c(1:length(CPG_tcga_ge_id))) {
  CPG_tcga_ge_id[[i]]<-sort(CPG_tcga_ge[i,])
}

CPG_tcga_ge_id

##只选取cpg island 关联只有一个基因的 
CPG_tcga_ge<-as.data.frame(t(sapply(CPG_tcga_ge_id,"[",i=1:max(sapply(CPG_tcga_ge_id,length)))))
CPG_tcga_gene<-cbind(CPG_tcga,CPG_tcga_ge)[which(sapply(CPG_tcga_ge_id,length)==1),] %>% dplyr::select(1:5,7:13) 
CPG_tcga_gene_p_order<-CPG_tcga_gene[order(CPG_tcga_gene$adj.P.Val),] 

## cpg island 差异化甲基化分析，上调和下调的
cpg_gene_down<-CPG_tcga_gene_p_order%>% filter(adj.P.Val<0.05,logFC<0) %>% pull(V2)
cpg_gene_up<-CPG_tcga_gene_p_order%>% filter(adj.P.Val<0.05,logFC>0) %>% pull(V2)
cpg_gene_id_up<-en2ENSE %>% filter(SYMBOL %in% cpg_gene_up) %>% pull(ENTREZID) %>% unique()
cpg_gene_id_down<-en2ENSE %>% filter(SYMBOL %in% cpg_gene_down) %>% pull(ENTREZID) %>% unique()

go_up_cpg_gene <- enrichGO(gene          = cpg_gene_id_up,
                         universe      = result_tmg$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "none",
                         minGSSize = 1,
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
View(go_up_cpg_gene@result)

exp(0.1)/(1+exp(0.1))
exp(0.9)/(1+exp(0.9))

go_down_cpg_gene <- enrichGO(gene          = cpg_gene_id_down,
                           universe      = result_tmg$ENTREZID,
                           OrgDb         = org.Hs.eg.db,
                           ont           = "ALL",
                           pAdjustMethod = "none",
                           minGSSize = 1,
                           pvalueCutoff  = 0.05,
                           readable      = TRUE)
View(go_down_cpg_gene@result)
CPG_is_gene<-CPG_tcga_gene %>% dplyr::select(Islands_Name,V2)

### RNA DE analysis和甲基化同时有差异的基因及数值
res_tmg_t_sig<-res_tmg_t %>% filter(padj<0.05)
inter<-intersect(c(cpg_gene_down,cpg_gene_up),res_tmg_t_sig$gene_name)
cpg_inter<-CPG_tcga_gene_p_order %>% filter(V2 %in% inter) %>% dplyr::select(logFC,UCSC_RefGene_Name=V2,adj.P.Val )
mrna_sig<-res_tmg_t_sig %>% filter(gene_name %in% inter) %>% dplyr::select(log2FoldChange,gene_name )
inter_df<-cpg_inter %>% left_join(mrna_sig,by=c("UCSC_RefGene_Name"="gene_name"))

##  cpg岛甲基化 mRNA co-significant散点图
ggplot(data = inter_df)+
  geom_point(aes(x=logFC,y=log2FoldChange)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_vline(xintercept = 0,linetype="dashed") +
  theme_classic(base_size=14) +
  xlab("logFC CPG Island methylation beta value") +
  ylab("logFC mRNA expression")+
  geom_point(data=subset(inter_df, (logFC<0 & log2FoldChange>0)),
             color="#E64B35FF",aes(x=logFC,y=log2FoldChange),size=3)+
  geom_point(data=subset(inter_df, (logFC>0 & log2FoldChange<0)),
             color="#3C5488FF",aes(x=logFC,y=log2FoldChange),size=3)

# +  geom_text_repel(data=subset(inter_df, ((logFC<0 & log2FoldChange>0)|(logFC>0 & log2FoldChange<0))),
 #                 aes(x=logFC,y=log2FoldChange,label=UCSC_RefGene_Name))

## 差异表达的mRNA和差异甲基化的CpG island呈相反趋势的基因列表
cpg_down_exp_up<-inter_df %>% filter(logFC<0,log2FoldChange>0)%>% pull(UCSC_RefGene_Name) %>% unique 
cpg_up_exp_down<-inter_df %>% filter(logFC>0,log2FoldChange<0) %>% pull(UCSC_RefGene_Name) %>% unique 


cpg_up_exp_up<-inter_df %>% filter(logFC>0,log2FoldChange>0)%>% pull(UCSC_RefGene_Name) %>% unique 
cpg_down_exp_down<-inter_df %>% filter(logFC<0,log2FoldChange<0) %>% pull(UCSC_RefGene_Name) %>% unique 


cpg_down_geid_up<-en2ENSE %>% filter(SYMBOL %in% cpg_down_exp_up) %>% pull(ENTREZID)  %>% unique
cpg_up_geid_down<-en2ENSE %>% filter(SYMBOL %in% cpg_up_exp_down) %>% pull(ENTREZID)  %>% unique

## 通路分析
cpg_down_geid_up_go <- enrichGO(gene          = cpg_down_geid_up,
                         universe      = result_tmg$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "none",
                         minGSSize = 1,
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
View(cpg_down_geid_up_go@result)

cpg_up_geid_down_go <- enrichGO(gene          = cpg_up_geid_down,
                                universe      = result_tmg$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "ALL",
                                pAdjustMethod = "none",
                                minGSSize = 1,
                                pvalueCutoff  = 0.05,
                                readable      = TRUE)
View(cpg_up_geid_down_go@result)

##ion channel pathway
#channel_tcga_me_path@result<-go_up_tmg_cc@result %>% filter(str_detect(go_up_tmg_cc@result$Description,pattern = "channel") |
##                                                           str_detect(go_up_tmg_cc@result$Description,pattern = "transport"))
#neuro related pathway  
#neuro_tcga_me_path@result<-go_up_tmg_cc@result %>%
#  filter(str_detect(go_up_tmg_cc@result$Description,pattern = "syna") |
#           str_detect(go_up_tmg_cc@result$Description,pattern = "neu")|
#           str_detect(go_up_tmg_cc@result$Description,pattern = "nerv")|
#           str_detect(go_up_tmg_cc@result$Description,pattern = "sarc")) %>%
#  filter(Count>1 | Description=="postsynaptic cytoskeleton" )



### 启动子与非启动子区域的CpG island 甲基化与关联基因mRNA表达的相关性比较

 ### 每个病人每个island的mean beta
me_da_isl_mean_beta<-me_da_isl%>%group_by(Islands_Name)%>%summarise_at(vars(colnames(me_da_isl)[10:130]),mean)  
 # me_da_isl_median_beta<-me_da_isl%>%group_by(Islands_Name)%>%summarise_at(vars(colnames(me_da_isl)[10:130]),median) 
 
 ## 所有病人每个island mean beta
me_da_isl_mean_beta$mean<-rowMeans(me_da_isl_mean_beta[,2:length(me_da_isl_mean_beta)])
me_da_isl_mean_beta_all<-me_da_isl_mean_beta %>% dplyr::select(Islands_Name,mean)
 #me_da_isl_mean_beta_all<-me_da_isl_mean_beta %>% dplyr::select(Islands_Name,"TCGA-XU-A936")

 ## 所有基因
CPG_is_gene<-CPG_is_gene %>% left_join(me_da_isl_mean_beta_all, by="Islands_Name")
colnames(CPG_is_gene)[c(2,3)]<-c("Gene","Beta")

 ## 所有患者每个基因平均表达表
ExpData_FPKM$mean<-rowMeans(ExpData_FPKM[,2:(length(ExpData_FPKM)-1)])
ExpData_FPKM_all<-ExpData_FPKM%>% dplyr::select(ensemblID,SYMBOL,mean)

 #ExpData_FPKM_all<-ExpData_FPKM%>% dplyr::select(ensemblID,SYMBOL,"TCGA-XU-A936")
 ## me_num_pisl_id:与启动子相关的island (计算过程见met_isl_model)
CPG_is_gene_mRNA<-CPG_is_gene %>% left_join(ExpData_FPKM_all, by = c("Gene" = "SYMBOL")) %>% 
  filter(Islands_Name %in% me_num_pisl_id)

CPG_is_gene_mRNA_not<-CPG_is_gene %>% left_join(ExpData_FPKM_all, by = c("Gene" = "SYMBOL")) %>% 
  filter(!(Islands_Name %in% me_num_pisl_id))


 ## 画图
promoter_CpG_plot<-ggplot(CPG_is_gene_mRNA, aes(Beta,mean)) + 
  geom_point(size =0.5,color="#468BCA")+  
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="#DB706C")+ 
  coord_cartesian(ylim = c(0,11))+ 
  labs(x="Methylation (Beta value)",y="mRNA (FPKM)") + 
  theme_classic() 
nonpromoter_CpG_plot<-ggplot(CPG_is_gene_mRNA_not, aes(Beta,mean)) + 
  geom_point(size =0.5,color="#468BCA")+   
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="#DB706C")+ 
  coord_cartesian(ylim = c(0,11)) + 
  labs(x="Methylation (Beta value)",y="mRNA (FPKM)") + 
  theme_classic()  
library(gridExtra)
grid.arrange(promoter_CpG_plot,nonpromoter_CpG_plot,ncol=2,widths=c(1,1))

 ## 线性回归统计分析
summary(lm(CPG_is_gene_mRNA$mean~CPG_is_gene_mRNA$Beta))
summary(lm(CPG_is_gene_mRNA_not$mean~CPG_is_gene_mRNA_not$Beta))

### DNMT1与CpG island 甲基化水平相关性图
me_da_isl_mean_beta

CPG_is_gene_is_all<-CPG_is_gene[,c(1,2)] %>% left_join(me_da_isl_mean_beta, by = c("Islands_Name")) %>% 
  filter(Islands_Name %in% me_num_pisl_id)


CPG_is_gene_is_lg_form<-data.table::melt(CPG_is_gene_is_all[,c(1:(length(colnames(CPG_is_gene_is_all))-1))], id.vars = 1:2, variable.name = "sample")

DNMT1_mRNA<-ExpData_FPKM %>% filter(SYMBOL=="DNMT1") %>% select(2:(length(colnames(ExpData_FPKM))-2)) %>% t() %>% 
  data.frame()
DNMT1_mRNA<-ExpData_FPKM %>% filter(SYMBOL=="DNMT3A") %>% select(2:(length(colnames(ExpData_FPKM))-2)) %>% t() %>% 
  data.frame()
colnames(DNMT1_mRNA)<-"DNMT1_mRNA"
DNMT1_mRNA$sample<-rownames(DNMT1_mRNA)
CPG_is_mRNA_all_lg<-CPG_is_gene_is_lg_form %>% left_join(DNMT1_mRNA,by="sample")


ggplot(CPG_is_mRNA_all_lg, aes(DNMT1_mRNA,value)) + 
  geom_point(size =0.001,color="#468BCA")+
  ggtitle(paste0("DNMT1 mRNA-CpG island methylation"))+ 
  theme_classic() 
summary(lm(CPG_is_mRNA_all_lg$DNMT1_mRNA~CPG_is_mRNA_all_lg$value))

plot(CPG_is_mRNA_all_lg$DNMT1_mRNA,CPG_is_mRNA_all_lg$value)


+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH")+
  stat_poly_eq(aes(label = paste(..eq.label..)), formula = y ~ x, parse = T,family = "SH",label.x = 0.95,label.y = 1.0)



###  单个CPG位点甲基化与mRNA表达量
  dplyr::select(1:9,132) %>% filter(!(.$UCSC_RefGene_Name=="")) %>% 
me_da_anno_ge<-me_da_anno %>% mutate(beta_mean=apply(.[,c(10:130)],1,mean))%>% 
  filter((!str_detect(.$UCSC_RefGene_Name,pattern=";")) )
me_da_anno_ge$Relation_to_Island[me_da_anno_ge$Relation_to_Island=="N_Shore"]<-"Shore"
me_da_anno_ge$Relation_to_Island[me_da_anno_ge$Relation_to_Island=="S_Shore"]<-"Shore"
me_da_anno_ge$Relation_to_Island[me_da_anno_ge$Relation_to_Island=="S_Shelf"]<-"Shelf"
me_da_anno_ge$Relation_to_Island[me_da_anno_ge$Relation_to_Island=="N_Shelf"]<-"Shelf"

me_ge_mRNA<-me_da_anno_ge %>% left_join(ExpData_FPKM_all, by = c("UCSC_RefGene_Name" = "SYMBOL"))

table(me_da_anno_ge$UCSC_RefGene_Group)

plot_me_mRNA<-function(x,y) {
promoter_me_is_plot<-ggplot(me_ge_mRNA[me_ge_mRNA$Relation_to_Island==x,], aes(beta_mean,mean)) + 
  geom_point(size =0.002,color="#468BCA")+  
  geom_smooth(method = "lm",se = F, show.legend = T,size=0.5,color="#DB706C")+
  coord_cartesian(ylim = c(0,11))+
  labs(x="Methylation (beta value)",y="mRNA (FPKM)") + ggtitle(y)+
  theme_classic()
promoter_me_is_plot
}

p1<-plot_me_mRNA("Island","Island")
me_ge_mRNA_isl<-me_ge_mRNA[me_ge_mRNA$Relation_to_Island=="Island",]
summary(lm(me_ge_mRNA_isl$mean~me_ge_mRNA_isl$beta_mean))
View(me_ge_mRNA_isl)
p2<-plot_me_mRNA("Shore","Shore")
me_ge_mRNA_sho<-me_ge_mRNA[me_ge_mRNA$Relation_to_Island=="Shore",]
summary(lm(me_ge_mRNA_sho$mean~me_ge_mRNA_sho$beta_mean))

p3<-plot_me_mRNA("Shelf","Shelf")
me_ge_mRNA_she<-me_ge_mRNA[me_ge_mRNA$Relation_to_Island=="Shelf",]
summary(lm(me_ge_mRNA_she$mean~me_ge_mRNA_she$beta_mean))

p4<-plot_me_mRNA("OpenSea","OpenSea")
me_ge_mRNA_ope<-me_ge_mRNA[me_ge_mRNA$Relation_to_Island=="OpenSea",]
summary(lm(me_ge_mRNA_ope$mean~me_ge_mRNA_ope$beta_mean))

library(gridExtra)
grid.arrange(p1,p2,p3,p4,ncol=4)

###启动子与非启动子相关的CpG island 与MG相关性 tSNE plot
## TCGA CPG island tSNE plot
## TCGA CPG island tSNE analysis
library(Rtsne)
library(lumi)

me_da_isl_mean_beta
me_da_isl_mean_beta_is<-me_da_isl_mean_beta %>% filter(Islands_Name %in% CPG_is_gene_mRNA_not$Islands_Name)

myNorm_M<-beta2m(me_da_isl_mean_beta_is[,2:122])
rownames(myNorm_M)<-me_da_isl_mean_beta_is$Islands_Name

myNorm_M_v<-rowVars(as.matrix(myNorm_M))
myNorm_M_t<-myNorm_M %>%  as.data.frame() %>% arrange(desc(myNorm_M_v))



set.seed(321) # 设置随机数种子
dim(myNorm_M_t)
tsne_out = Rtsne(
  t(myNorm_M_t[1:1000,]),
  dims = 2,
  pca = T,
  max_iter = 5000,
  theta = 0.1,
  perplexity = 10,
  verbose = F
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
row.names(tsne_result)<-met_meta$patient_id

ggplot(tsne_result,aes(tSNE1,tSNE2,color=met_meta$MG)) +
  geom_point(size=2) +
  stat_ellipse(level = 0.5) +
  scale_color_npg() +
  theme_classic()


ggplot(tsne_result,aes(tSNE1,tSNE2,color=met_meta$histo_type)) +
  geom_point(size=2) +
  stat_ellipse(level = 0.5) +
  scale_color_npg() +
  theme_classic()






