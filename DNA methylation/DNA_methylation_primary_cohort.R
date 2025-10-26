library("ChAMP")
library(tidyverse)
library("ggplot2")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(DMRcatedata)


## annotation
data(Islands.UCSC)
data(Locations)
data(Other)

#methylation

setwd("D:/Documents/TMG_methylation_paper")
Methy_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/legacy/TCGA.THYM.sampleMap_HumanMethylation450/HumanMethylation450")
colnames(Methy_THYM)[2:length(Methy_THYM)]<-substr(colnames(Methy_THYM)[2:length(Methy_THYM)],1,12)

##clinical data
clinical_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/harmonized/Clinical/Clinical_Supplement/7cca5722-26cf-4ac8-a4a6-b803459f1861/nationwidechildrens.org_clinical_patient_thym_2.txt")
#clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-3G-AB14"]<-"NO"
#clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-X7-A8DF"]<-"NO"
#View(clinical_THYM)
meta_data<-clinical_THYM %>% dplyr::select(patient_id=bcr_patient_barcode,
                                    gender,height,weight,
                                    histo_type = histological_type, 
                                    MG = history_myasthenia_gravis,  
                                    section_MG = section_myasthenia_gravis,  
                                    age_patho_diagnosis =   
                                      age_at_initial_pathologic_diagnosis,  
                                    masaoka_stage) %>%  
  dplyr::filter(MG == "YES"|MG =="NO" ) %>% 
  arrange(patient_id) 
met_meta<-meta_data %>% dplyr::filter(meta_data$patient_id %in% (colnames(Methy_THYM)[-1]))  
 
met_meta$histo_type[met_meta$histo_type=="Thymoma; Type B2|Thymoma; Type B3"]<-'Thymoma; Type B3' 
met_meta$histo_type[met_meta$histo_type=="Thymoma; Type A|Thymoma; Type AB"]<-'Thymoma; Type AB' 
met_meta$histo_type[met_meta$histo_type=="Thymoma; Type B1|Thymoma; Type B2"]<-'Thymoma; Type B2' 

### table of TCGA thymoma clinical characteristics of included patients  
table(met_meta$gender) 
table(met_meta$MG)  
table(met_meta$histo_type) 
table(met_meta$masaoka_stage)  
sort(as.numeric(met_meta$age_patho_diagnosis[-34]))  


## omit the methylation loci that has NAs 
met_THYM<- Methy_THYM %>% dplyr::select("sample",met_meta$patient_id) 
length(met_THYM$sample) 
met_THYM_m <- met_THYM %>% dplyr::select(-sample) %>% as.matrix() 
row.names(met_THYM_m)<-met_THYM$sample 
met_THYM_m<-na.omit(met_THYM_m) 

##filter 
# myload_met <- champ.filter(beta=met_THYM_m,pd=met_meta,fixOutlier = T,filterBeads = F,filterDetP = F,autoimpute = F) 

##quality control 
#champ.QC(beta=myload_met$beta,pheno =myload_met$pd$patient_id) 

## normalization 
#myNorm_met <- champ.norm(beta = myload_met$beta,arraytype="450K", cores=8) 
#write.csv(myNorm_met,"tcga_met_normalization_TCGA_database_20220813.csv") 
 
## normalized methylation data 
myNorm_met<-read.csv("./analy_data/tcga_met_normalization_TCGA_database_20220813.csv",check.names = F) 
row.names(myNorm_met)<-myNorm_met[,1] 
myNorm_met<-myNorm_met[,-1] %>% as.matrix() 
myNorm_met<-na.omit(myNorm_met) 


# Differential methylation position (DMP) analysis 

## champ DMP analysis 
myDMP<-champ.DMP(beta = myNorm_met,pheno=met_meta$MG, 
          #cov1= met_meta$histo_type, 
          adjust.method = "BH", adjPVal = 1) 
myDMP1<-myDMP$NO_to_YES
myDMP1 %>% filter(gene == "TTN") 
myDMP1 %>% filter(gene == "CHRNA1") 
myDMP1 %>% filter(gene == "NEFM") 
myDMP1 %>% filter(gene == "RYR1") 
myDMP1  %>% filter(gene == "RYR2") 
myDMP1 %>% filter(gene == "RYR3") 
myDMP1  %>% filter(gene == "CHRND") 
myDMP1 %>% filter(gene == "CHRNG") 
myDMP1 %>% filter(gene == "CHRNE") 
myDMP1  %>% filter(gene == "CHRNB1") 
myDMP1 %>% filter(gene == "BACH2") 
myDMP1 %>% filter(gene == "KCNA4") 




# 
# ## limma DMP analysis
# library(limma)
# 
# compare.group <- unique(met_meta$MG)[1:2]
# 
# p <- met_meta$MG[which(met_meta$MG %in% compare.group)]
# beta <- myNorm_met[,which(met_meta$MG %in% compare.group)]
# design <- model.matrix( ~ 0+ factor(p)+factor(met_meta$histo_type))
# colnames(design) <- c("control", "case","typeAB","typeB1","typeB2","typeB3","typeC")
# contrast.matrix <- makeContrasts(contrasts=paste(colnames(design)[2:1],collapse="-"), levels = design)
# print(contrast.matrix)
# 
# fit <- lmFit(beta, design)
# fit2 <- contrasts.fit(fit,contrast.matrix)
# fit3 <- eBayes(fit2)
# myDMP2 <- topTable(fit3,coef=1,number=nrow(beta),adjust.method="BH")
# myDMP2 <- myDMP2 %>% mutate(Met_ID=rownames(.)) %>% dplyr::select(Met_ID,everything()) %>%
#   left_join(m450anno[,c(1:6,8:10)],by="Met_ID") %>%
#   dplyr::select(Met_ID, chr, pos, strand,Islands_Name,
#          Relation_to_Island,UCSC_RefGene_Name, everything())
# 
# myDMP2 %>% filter(adj.P.Val<0.05)


######################################################
citation("ChAMP")
### DMR analysis
 # get M value
 M_myNorm_met<-logit2(myNorm_met)


 # dmrcate CPG annotation
group<-factor(met_meta$MG)
design<-model.matrix(~group)

myannotation <- cpg.annotate(datatype = "array",
                             object = M_myNorm_met, design = design, coef = ncol(design), fdr = 1,pcutoff = 1,
                             analysis.type = "differential", arraytype="450K",
                             annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                             what = "M")

# dmrcate DMR analysis and plot NEFM

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2,pcutoff = 1)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19") ## 有的网络连接不上
results.ranges[which(results.ranges$overlapping.genes=="AIRE"),] # 5488

# which(results.ranges$overlapping.genes=="NEFM") # 5488
# 
# results.ranges
# 
# groups <- c(yes="magenta", no="forestgreen")
# cols <- groups[meta_data_thymoma_kajiura$Myasthenia_gravis]
# plot<-DMR.plot(ranges=results.ranges, dmr=5488, CpGs=myNorm_met_gse, what="Beta",
#                arraytype = "450K", genome="hg19",phen.col=cols)
# 
# ## significant DMP genes and RNA expression genes
# 
# library(reshape2)
# library(org.Hs.eg.db)
# k<-keys(org.Hs.eg.db,keytype='ENSEMBL') 
# en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,columns=c("ENTREZID",'SYMBOL'), keytype = 'ENSEMBL') 
# 
# 
# DMP_gene_list<-unique(myDMP$gene[myDMP$adj.P.Val<0.05],)
# DMP_down<-myDMP1 %>% filter(logFC<0 & adj.P.Val<0.05) %>% pull(gene) %>% unique(.)
# DMP_up<-myDMP1 %>% filter(logFC>0 & adj.P.Val<0.05) %>% pull(gene)  %>% unique(.)
# 
# DMP_ge_down<-en2ENSE %>% filter(SYMBOL %in% DMP_down)
# DMP_ge_up<-en2ENSE %>% filter(SYMBOL %in% DMP_up)
# 
# DMP_ge_down_path<-enrichGO(gene          = DMP_ge_down$ENTREZID,
#                         universe      = na.omit(result_tmg$ENTREZID),
#                         OrgDb         = org.Hs.eg.db,
#                         ont           = "ALL",
#                         pAdjustMethod = "none",
#                         minGSSize = 1,
#                         maxGSSize = 150,
#                         pvalueCutoff  = 0.05,
#                         #qvalueCutoff  = 0.05,
#                         readable      = TRUE)
# View(DMP_ge_down_path@result)
# 
# DMP_ge_up_path<-enrichGO(gene          = DMP_ge_up$ENTREZID,
#                            universe      = na.omit(result_tmg$ENTREZID),
#                            OrgDb         = org.Hs.eg.db,
#                            ont           = "ALL",
#                            pAdjustMethod = "fdr",
#                            minGSSize = 1,
#                            maxGSSize = 150,
#                            pvalueCutoff  = 0.05,
#                            qvalueCutoff  = 0.05,
#                            readable      = TRUE)
# View(DMP_ge_up_path@result)
# 
# RNA_gene_list<-unique(result_tmg$gene_name[result_tmg$padj<0.05])
# RNA_ge_up<-result_tmg  %>% filter(padj<0.05 & log2FoldChange>0)  %>% pull(gene_name)
# RNA_ge_down<-result_tmg %>% filter(padj<0.05 & log2FoldChange<0) %>% pull(gene_name)
# 
# both_gene<-intersect(DMP_gene_list,RNA_gene_list)
# both_DMP_d<-intersect(DMP_down,RNA_ge_up)
# both_DMP_u<-intersect(DMP_up,RNA_ge_down) # return null
# 
# DMP_ge_ent<-en2ENSE %>% filter(SYMBOL %in% DMP_gene_list) %>% dplyr::select(ENTREZID) %>% pull()
# both_gene_entr<-en2ENSE %>% filter(SYMBOL %in% both_DMP_d) %>% dplyr::select(ENTREZID) %>% pull()
# both_DMP_d_entr<-en2ENSE %>% filter(SYMBOL %in% both_gene) %>% dplyr::select(ENTREZID) %>% pull()
# 
# both_DMP_d_cc<-enrichGO(gene          = both_gene_entr,
#          universe      = na.omit(result_tmg$ENTREZID),
#          OrgDb         = org.Hs.eg.db,
#          ont           = "ALL",
#          pAdjustMethod = "fdr",
#          minGSSize = 1,
#          maxGSSize = 150,
#          pvalueCutoff  = 0.05,
#          #qvalueCutoff  = 0.05,
#          readable      = TRUE)
# 
# both_DMP_d_bp<-enrichGO(gene          = both_gene_entr,
#                         universe      = na.omit(result_tmg$ENTREZID),
#                         OrgDb         = org.Hs.eg.db,
#                         ont           = "BP",
#                         pAdjustMethod = "fdr",
#                         minGSSize = 1,
#                         maxGSSize = 150,
#                         pvalueCutoff  = 0.05,
#                         qvalueCutoff  = 0.05,
#                         readable      = TRUE)
# 
# View(both_DMP_d_cc@result)
# View(both_DMP_d_bp@result)
# barplot(both_gene,showCategory = 47)
# intersect(DMP_gene_list,result_tmg_up$gene_name)
# 
# #################################################################################################################
# ## tSNE analysis and plot
# library(Rtsne)
# library(lumi)
# 
# myNorm_M<-beta2m(myNorm_met)
# 
# myNorm_M_v<-rowVars(myNorm_M)
# myNorm_M_t<-myNorm_M %>%  as.data.frame() %>% arrange(desc(myNorm_M_v))
# 
# set.seed(321) # 设置随机数种子 
# dim(myNorm_M_t) 
#  
# tsne_out = Rtsne( 
#   t(myNorm_M_t[1:10000,]), dims = 2, 
#   pca = T, max_iter = 5000, theta = 0.1,  
#   perplexity = 10, verbose = F)  
# 
# tsne_result = as.data.frame(tsne_out$Y)  
# colnames(tsne_result) = c("tSNE1","tSNE2")  
# row.names(tsne_result)<-met_meta$patient_id  
# 
# p<-ggplot(tsne_result,aes(tSNE1,tSNE2,color=met_meta$MG)) +
#   geom_point(size=1.5)+ scale_color_npg() + 
#   stat_ellipse(level = 0.5) 
# p + theme_classic() 
# 
# 
#  
# ggplot(tsne_result,aes(tSNE1,tSNE2,color=met_meta$histo_type)) +  
#   geom_point(size=1.5)  + scale_color_npg() +  theme_classic() +  
#   stat_ellipse(level = 0.5)
# 
# write.csv(tsne_result,"tsne_tcga_met.csv")
