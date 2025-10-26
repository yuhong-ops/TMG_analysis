library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)
library(ggfortify)
library(ggrepel)
library(ggplot2)
library(ggalt)
library(dplyr)
library(ggsci)
library(org.Hs.eg.db)

setwd("D:/Documents/TMG_methylation_paper") #laptop
#setwd("C:/Users/Administrator/Documents") # desktop

source('./script/functions.R')  
R.utils::setOption("clusterProfiler.download.method",'auto')  

#TCGA thymus RNA  
RNA_expr_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/legacy/Gene_expression/Gene_expression_quantification_UCSC_htseq/TCGA-THYM.htseq_counts.tsv")

RNA_expr_THYM<-RNA_expr_THYM[-c(60484:60488),]  
   #log2(count+1) tranform to raw count  
RNA_expr_THYM[,2:122]<-2^RNA_expr_THYM[,2:122]-1  
 
RNA_expr_THYM$Ensembl_ID<-substr(RNA_expr_THYM$Ensembl_ID,1,15)  

colnames(RNA_expr_THYM)[1]<-"gene_id" 
RNA_expr_THYM <- RNA_expr_THYM %>% dplyr::select(order(colnames(RNA_expr_THYM))) 

 #remove normal tissue 
RNA_expr_THYM<-RNA_expr_THYM %>% dplyr::select(-c("TCGA-X7-A8D7-11A","TCGA-X7-A8D6-11A"))
colnames(RNA_expr_THYM)[2:length(colnames(RNA_expr_THYM))]<-substr(colnames(RNA_expr_THYM)[2:length(colnames(RNA_expr_THYM))],1,12)


#clinical

clinical_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/harmonized/Clinical/Clinical_Supplement/7cca5722-26cf-4ac8-a4a6-b803459f1861/nationwidechildrens.org_clinical_patient_thym_2.txt",
)

#update myasthenia gravis information

##clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-3G-AB14"]<-"NO"
#clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-X7-A8DF"]<-"NO"
#clinical_THYM$history_myasthenia_gravis

meta_data<-clinical_THYM %>% dplyr::select(patient_id=bcr_patient_barcode,
                                    gender,height,weight,
                                    histo_type = histological_type,
                                    MG = history_myasthenia_gravis,
                                    section_MG = section_myasthenia_gravis,
                                    age_patho_diagnosis = 
                                      age_at_initial_pathologic_diagnosis) %>%
                              filter(MG == "YES"|MG =="NO" ) %>%
                              arrange(patient_id)



meta_data_RNA<-meta_data %>% filter(meta_data$patient_id %in% (colnames(RNA_expr_THYM)[-1]))



meta_data_RNA$histo_type[meta_data_RNA$histo_type=="Thymoma; Type B2|Thymoma; Type B3"]<-'Thymoma; Type B3'
meta_data_RNA$histo_type[meta_data_RNA$histo_type=="Thymoma; Type A|Thymoma; Type AB"]<-'Thymoma; Type AB'
meta_data_RNA$histo_type[meta_data_RNA$histo_type=="Thymoma; Type B1|Thymoma; Type B2"]<-'Thymoma; Type B2'

table(meta_data_RNA$MG)

#TCGA_paper_clinical_data

#clinical_THYM_paper<-read_csv(file="./TCGApaper/TCGA database/TCGApaper/TCGA_paper.csv")


#看TCGA文章里面和数据库里面标本的临床资料是否一???

#clinical_THYM_paper$patient_short_barcode %in% clinical_THYM$bcr_patient_barcode

#clinical_THYM_paper$patient_short_barcode %in% clinical_THYM$bcr_patient_barcode

#setdiff(clinical_THYM_paper$patient_short_barcode,meta_data$patient_id)
#"TCGA-3G-AB14" 1 MG NO dead 2 MG unknown alive  
#"TCGA-5K-AAAP" 1 MG unknown 2 MG unknown 一致
#"TCGA-X7-A8DF" 1 MG NO      2 MG not available
# 根据TCGA paper内容更新

#setdiff(meta_data$patient_id,clinical_THYM_paper$patient_short_barcode)

#intersect(clinical_THYM_paper$patient_short_barcode,

#intersect(clinical_THYM_paper$patient_short_barcode,
#          meta_data$patient_id)

#intersect(clinical_THYM_paper$patient_short_barcode[clinical_THYM_paper$history_myasthenia_gravis=="YES"],
#          meta_data$patient_id[meta_data$MG == "YES"])

#intersect(clinical_THYM_paper$patient_short_barcode[clinical_THYM_paper$history_myasthenia_gravis=="NO"],
#          meta_data$patient_id[meta_data$MG == "NO"])



## RNA that has data
RNA_expr_THYM_meta<- RNA_expr_THYM %>% dplyr::select("gene_id",meta_data_RNA$patient_id)

## geneName import
geneNames<-read.csv("./TCGApaper/geneNames.csv")%>% 
  mutate(gene_id = substr(.$gene_id,1,15) ) # remove ensemble gene_id version

#mitogenes
MitoGenes <- geneNames[geneNames$seqnames=="chrM",]
RiboGenes <- geneNames[geneNames$gene_type %in% c("rRNA", "Mt_rRNA"),]
ProtGenes <- geneNames[geneNames$gene_type=="protein_coding",]

# Now we remove the chr names from the geneNames because there are duplicate
# gene names in chrX and chrY 
geneNames %<>% dplyr::select(gene_id, gene_name, gene_type) %>%
  unique 


## 2 Quality control
## 2.1 Transcript biotypes

# Counts for mitochondrial and rRNA genes
MitoCountSum <- colSums(RNA_expr_THYM_meta %>% filter(gene_id %in% MitoGenes$gene_id) %>% dplyr::select(-gene_id))
MitoCountSum <- data.frame(patient_id=names(MitoCountSum), mito_counts=round(MitoCountSum), stringsAsFactors=FALSE)
RiboCountSum <- colSums(RNA_expr_THYM_meta %>% filter(gene_id %in% RiboGenes$gene_id) %>% dplyr::select(-gene_id))
RiboCountSum <- data.frame(patient_id=names(RiboCountSum), ribo_counts=round(RiboCountSum), stringsAsFactors=FALSE)
ProtCountSum <- colSums(RNA_expr_THYM_meta %>% filter(gene_id %in% (ProtGenes %>% filter(seqnames!="chrM"))$gene_id) %>% dplyr::select(-gene_id))
ProtCountSum <- data.frame(patient_id=names(ProtCountSum), prot_counts=round(ProtCountSum), stringsAsFactors=FALSE)
TotCountSum  <- colSums(RNA_expr_THYM_meta %>% dplyr::select(-gene_id))
TotCountSum <- data.frame(patient_id=names(TotCountSum), counts=TotCountSum, stringsAsFactors=FALSE)

ReadFractions <- inner_join(MitoCountSum, RiboCountSum, by="patient_id") %>%
  inner_join(ProtCountSum, by="patient_id") %>%
  inner_join(TotCountSum, by="patient_id") %>%
  mutate(mito_fraction=mito_counts/counts, ribo_fraction=ribo_counts/counts, prot_fraction=prot_counts/counts,
         other_fraction=1-(mito_fraction+ribo_fraction+prot_fraction), lib_size=as.integer(round(counts))) %>%
  dplyr::select(patient_id, mito_fraction, ribo_fraction, prot_fraction, other_fraction, lib_size=counts)

meta_data_RNA %<>% left_join(ReadFractions,by="patient_id")

ReadFractions %>% dplyr::select(-lib_size) %>% gather("Gene_biotype", "Proportion", -patient_id) %>%
  arrange(Proportion) %>% mutate(Gene_biotype=factor(Gene_biotype, levels=c("ribo_fraction","mito_fraction","other_fraction","prot_fraction"))) %>%
  ggplot(aes(x=patient_id, y=Proportion, fill=Gene_biotype)) +
  geom_bar(position="stack", stat="identity", colour="white") +
  labs(y="Proportion of reads", x="Sample") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())



##2.2 Highly expressed genes
## highly expressed 

# Counts excluding mitochondrial genes and mitochondrial reads fraction over
# total lib. size
minProp <- 0.01
MitoCountFiltered <- RNA_expr_THYM_meta %>% filter(!gene_id %in% MitoGenes$gene_id)
CountFreqs <- prop.table(as.matrix(MitoCountFiltered[,-1]), 2)
rownames(CountFreqs) <- MitoCountFiltered$gene_id

#CountFreqs[CountFreqs < minProp] <- NA

HighExpGenes <- CountFreqs[apply(CountFreqs, 1, function(x) any(x >= minProp)),] %>%
  as.data.frame %>% mutate(gene_id=rownames(.)) %>%
  left_join(geneNames) %>% 
  gather("patient_id", "Proportion", -gene_id, -gene_name, -gene_type)

s_order <- HighExpGenes %>% group_by(patient_id) %>% summarise(sum=sum(Proportion)) %>%
  arrange(sum) %>% pull(patient_id) 

g_order <- HighExpGenes %>% group_by(gene_name) %>% summarise(sum=sum(Proportion)) %>%
  arrange(sum) %>% pull(gene_name) 

HighExpGenes$patient_id <- factor(HighExpGenes$patient_id, levels=s_order)
HighExpGenes$gene_name <- factor(HighExpGenes$gene_name, levels=g_order)

ggplot(HighExpGenes, aes(patient_id, Proportion, fill=gene_name)) +
  geom_bar(stat = "identity", colour="white") +
  labs(y="Proportion of reads", x="Patient") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())


##Now that we have checked for highly abundant transcripts, we can transform the count data to log2(CPM) and do some more filtering on the non-mitochondrial genes. Genes are filtered out if:

##  no Gene Symbol (or duplicated)
##  expression is below the median expression for more than 20% of the samples

MitoCountFiltered<-RNA_expr_THYM_meta %>% filter(!gene_id %in% MitoGenes$gene_id)

countMatrixFiltered <- MitoCountFiltered

GenesKept <- rep(0, 3)
names(GenesKept) <- c("Total nuclear transcripts", "Trancripts above threshold", "With Gene Symbol")

# Create log2 CPM matrix after removal of mitochondria-encoded genes
cpmMatrixFiltered <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), log2(Count2CPM(countMatrixFiltered[,-1])+1),check.names=F)
GenesKept[1] <- nrow(cpmMatrixFiltered)

# Add gene symbols
ExpDataCPM <- left_join(cpmMatrixFiltered, geneNames %>% dplyr::select(gene_id, gene_name)) %>%
  dplyr::select(ensemblID=gene_id, GeneSymbol=gene_name, everything())

# Filter low-expressed genes similar to https://f1000research.com/articles/5-1438
# (although a bit stricter: keep genes that have at least 10 reads in 25% of
# the samples)
min_lib_size <- min(meta_data_RNA$lib_size)
min_cpm <- round(10/(min_lib_size/1e+6),1)
keep <- rowSums(cpmMatrixFiltered[,-1] > min_cpm) >= round(nrow(meta_data_RNA)*.25)
#keep %>% table

# Filter genes below noise level
ExpDataCPM <- ExpDataCPM [keep,]
GenesKept[2] <- nrow(ExpDataCPM)

# Filter out genes with no Gene Symbol or with duplicated Gene Symbol (highest
# expressed option of duplicates is kept)
ExpDataCPM %<>% filter(GeneSymbol != "", !is.na(GeneSymbol)) %>%
  mutate(Sum = do.call(pmax, select_if(., is.numeric))) %>%
  arrange(desc(Sum)) %>% 
  distinct(GeneSymbol, .keep_all=TRUE) %>%
  dplyr::select(-Sum) %>%
  arrange(ensemblID)

GenesKept[3] <- nrow(ExpDataCPM)
prot_coding_n <- c(geneNames %>% filter(gene_id %in% ExpDataCPM$ensemblID) %>%
                     pull(gene_type) %>% table %>% prop.table)["protein_coding"]


GenesKept %>% kable(caption="Number of genes passing filters") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

# Reorder samples in CPM matrix
colnames(ExpDataCPM)[3:length(colnames(ExpDataCPM))] == meta_data_RNA$patient_id


if ( ! all(colnames(ExpDataCPM)[-c(1,2)] == meta_data_RNA$patient_id) ) {
  warning("Reordering CPM matrix...")
  ExpDataCPM <- ExpDataCPM[,c(1,2,(match(meta_data_RNA$patient_id, colnames(ExpDataCPM[,-c(1,2)]))+2))]
} else {
  message("No need to reorder CPM matrix")
}

# Filter low-expressed genes from counts matrix and reorder if necessary
Counts <- data.frame(countMatrixFiltered,check.names=F) %>%
  filter(gene_id %in% ExpDataCPM$ensemblID) %>%
  mutate_if(is.numeric, function(x) as.integer(round(x)))

if ( ! all(colnames(Counts)[-1] == meta_data_RNA$patient_id) ) {
  warning("Reordering count matrix...")
  Counts <- Counts[,c(1,(match(meta_data_RNA$patient_id, colnames(Counts[,-1]))+1))]
} else {
  message("No need to reorder count matrix")
}



##heatmap of sample correlations


library(ComplexHeatmap)



sc <- ExpDataCPM %>%
  filter(ensemblID %in% ProtGenes$gene_id) %>%
  dplyr::select(-ensemblID, -GeneSymbol) %>%
  as.matrix %>%
  cor
medians <- matrixStats::rowMedians(sc, na.rm=TRUE)
names(medians) <- colnames(sc)

# Create temporal Metadata data.frame
dfMeta <- meta_data_RNA %>% dplyr::select(gender, histo_type, MG)
rownames(dfMeta) <- meta_data_RNA$patient_id

cols_sex <- c(MALE="darkturquoise", FEMALE="palevioletred")

column_ha = HeatmapAnnotation(MG=dfMeta$MG,
                              Histology= dfMeta$histo_type,
                              Sex=dfMeta$gender,
                              na_col="white",
                              col=list(Sex=cols_sex))

dend1 = cluster_within_group(sc, dfMeta$MG)
diag(sc) <- NA
Heatmap(sc, column_title="Sample correlation in gene expression", col=viridis::viridis(10), top_annotation=column_ha,
        show_row_names=FALSE, show_column_names=F, show_row_dend=FALSE, column_dend_height=unit(1,"cm"),
        cluster_columns = dend1,cluster_rows = dend1)

?viridis
# 
####################################################################################################################
library("DESeq2")
# DE analysis

Mt <- meta_data_RNA %>% filter(MG == "NO" | MG == "YES")
Mt %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

Mt$histo_type<-factor(Mt$histo_type)
Mt$MG<-factor(Mt$MG)

Ct <- Counts[,colnames(Counts) %in% Mt$patient_id] %>% na.omit(.)
rownames(Ct) <- Counts$gene_id

Md <- as.formula(paste0("~histo_type + MG"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_tmg_t<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T, pAdjustMethod = "bonferroni") %>%
  data.frame() %>% mutate(gene_id=rownames(.)) %>% dplyr::arrange(pvalue) %>%
  left_join(geneNames,by="gene_id")

res_tmg_t_fdr<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T, pAdjustMethod = "fdr") %>%
  data.frame() %>% mutate(gene_id=rownames(.)) %>% dplyr::arrange(pvalue) %>%
  left_join(geneNames,by="gene_id")

#res_tmg_t<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T, pAdjustMethod = "fdr") %>%
#  data.frame() %>% mutate(gene_id=rownames(.)) %>% dplyr::arrange(pvalue) %>%
#  left_join(geneNames,by="gene_id")



res_tmg_t[res_tmg_t$gene_name=='CHRNA1',]
res_tmg_t[res_tmg_t$gene_name=='AIRE',]
res_tmg_t[res_tmg_t$gene_name=='RYR3',]
res_tmg_t[res_tmg_t$gene_name=='GABRA5',]

res_tmg_t[res_tmg_t$gene_name=='TTN',]
res_tmg_t[res_tmg_t$gene_name=='NEFM',]

res_tmg_t[res_tmg_t$gene_name=='DNMT1',]
res_tmg_t[res_tmg_t$gene_name=='UHRF1',]
res_tmg_t[res_tmg_t$gene_name=='HDAC1',]
res_tmg_t[res_tmg_t$gene_name=='CTCF',]
res_tmg_t[res_tmg_t$gene_name=='SMCHD1',]
res_tmg_t[res_tmg_t$gene_name=='MPHOSPH8',]
res_tmg_t[res_tmg_t$gene_name %in% 
            c('DNMT1','HDAC1','CTCF','UHRF1','KMT2A','CTCF',
              'BAZ2A','SMCHD1','ARID4B','MPHOSPH8','ATF7IP',
              'TASOR','PPHLN1','RESF1','GATA3','TOX','RBM33',
              'BAZ2A','HDAC2','JARID2','RBBP4','SNRPD2','MTF2',
              'SUZ12','BOD1L1','KANSL1'),]

res_tmg_t[res_tmg_t$gene_name=='DNMT3A',]
res_tmg_t[res_tmg_t$gene_name=='DNMT3B',]
res_tmg_t[res_tmg_t$gene_name=='TET1',]
res_tmg_t[res_tmg_t$gene_name=='TET2',]
res_tmg_t[res_tmg_t$gene_name=='TET3',]
res_tmg_t[res_tmg_t$gene_name=='MECP2',]
res_tmg_t[res_tmg_t$gene_name=='CHRNG',]





res_tmg_t %>% filter(res_tmg_t$padj<0.05 & res_tmg_t$log2FoldChange>0)
res_tmg_t %>% filter(res_tmg_t$padj<0.05 & res_tmg_t$log2FoldChange<0)
res_tmg_t
res_tmg_t_sig<- res_tmg_t %>% filter(res_tmg_t$padj<0.05)
res_tmg_t_sig<- res_tmg_t %>% filter(res_tmg_t$padj<0.2)
write.csv(res_tmg_t,"TAMG_vs_nonMG_TCGA2.csv")


# 
# 
# length(res_tmg_t %>% filter(res_tmg_t$padj<0.05 & res_tmg_t$log2FoldChange>0) %>% pull(baseMean))
# length(res_tmg_t %>% filter(res_tmg_t$padj<0.05 & res_tmg_t$log2FoldChange<0) %>% pull(baseMean))
# 
# ### methylation related genes
res_tmg_t[res_tmg_t$gene_name=='DNMT1',]
# res_tmg_t[res_tmg_t$gene_name=='UHRF1',]
res_tmg_t[res_tmg_t$gene_name=='DNMT3A',]
# res_tmg_t[res_tmg_t$gene_name=='DNMT3B',]
res_tmg_t[res_tmg_t$gene_name=='DNMT3C',]
res_tmg_t[res_tmg_t$gene_name=='DNMT3L',]
res_tmg_t[res_tmg_t$gene_name=='TET1',]
res_tmg_t[res_tmg_t$gene_name=='TET2',]
res_tmg_t[res_tmg_t$gene_name=='TET3',]
res_tmg_t[res_tmg_t$gene_name=='AIRE',]
res_tmg_t[res_tmg_t$gene_name=='CHRNA1',]


# 
# 
 colData(dds)
# 
# #################################################################################################################
 ## PCA analyusis
# 
 library(limma)

 
 Mt$histo_type<-factor(Mt$histo_type)
 Mt$MG<-factor(Mt$MG)
 
 design<-model.matrix(~Mt$MG)
 v <- voom(Ct,design,normalize="quantile") 
 
 fit <- lmFit(v, design)
 fit <- eBayes(fit)
 topTable<-topTable(fit, coef="Mt$MGYES",number=length(row.names(Ct)))
 
 topTable%>% filter(row.names(topTable)=="ENSG00000138435")
 
 
 rld<-vst(dds)   #rld_2<-rlog(dds) slow
 plotPCA(rld,intgroup = "MG") +stat_ellipse(level = 0.9)
 
 
library(RColorBrewer)
 pal <- brewer.pal(8,"Dark2") 
 plotMDS(rld@assays@data@listData[[1]], top=500, gene.selection="common", 
         col=pal[Mt$MG],pch=19) 
 
 ## tSNE analysis
 library(Rtsne)
 rld_t<-rld@assays@data@listData[[1]]
 rld_t_v<-rowVars(rld_t)
 rld_t<-rld_t %>% cbind(rld_t_v) %>% as.data.frame() %>% arrange(desc(rld_t_v)) %>% dplyr::select(-rld_t_v)
# 
# 
 set.seed(321) # 设置随机数种子
# 
 tsne_out = Rtsne(
   t(rld_t[1:1000,]),
   dims = 2,
   pca = T,
   max_iter = 4000,
   theta = 0.1,
   perplexity = 20,
   verbose = F
 ) 
 

 tsne_result = as.data.frame(tsne_out$Y)
 colnames(tsne_result) = c("tSNE1","tSNE2")
 row.names(tsne_result)<-Mt$patient_id
 
 #write.csv(tsne_result,"tsne_tcga.csv")
 
 p<-ggplot(tsne_result,aes(tSNE1,tSNE2,color=Mt$MG)) +
   geom_point(size=2) +stat_ellipse(level = 0.5) +theme_classic()
 p+scale_color_npg()
 
 
 p<-ggplot(tsne_result,aes(tSNE1,tSNE2,color=Mt$histo_type)) +
   geom_point(size=2) + stat_ellipse(level = 0.5) +theme_classic()
 p+scale_color_npg()
# 
# 
# 
# #######################################################################################################################
# #important genes foldchange
# 
# lg2fc_chrna1<-res_tmg_t[res_tmg_t$gene_name=='CHRNA1',]$log2FoldChange
# 
# res_tmg_t[res_tmg_t$gene_name=='CHRNA1',] 
# res_tmg_t[res_tmg_t$gene_name=='CHRNG',] 
# res_tmg_t[res_tmg_t$gene_name=='GABRA5',] 
# res_tmg_t[res_tmg_t$gene_name=='MAP2',] 
# 
# res_tmg_t[res_tmg_t$gene_name=='NEFM',]# "similiarity with CHRNA1"
# res_tmg_t[res_tmg_t$gene_name=='NEFL',] 
# 
# res_tmg_t[res_tmg_t$gene_name=='RYR1',] 
# res_tmg_t[res_tmg_t$gene_name=='RYR2',] 
# res_tmg_t[res_tmg_t$gene_name=='RYR3',] 
# res_tmg_t[res_tmg_t$gene_name=='GRIK4',] 
# citation("DESeq2")
# 
# res_tmg_t[res_tmg_t$gene_name=='CHRNE',] 
# Counts[Counts$gene_id=="ENSG00000138435",]## CHRNA1 foldchange
# res_tmg_t[res_tmg_t$gene_name=='NEFM',] 
# 
# 
# Counts[Counts$gene_id=="ENSG00000198838",] ### RYR3
# 
# 
# chrna1_counts_mg<-Counts %>% filter(gene_id=="ENSG00000138435") %>% dplyr::select(meta_data_RNA$patient_id[meta_data_RNA$MG=="YES"]) 
# chrna1_counts_nmg<-Counts %>% filter(gene_id=="ENSG00000138435") %>% dplyr::select(meta_data_RNA$patient_id[meta_data_RNA$MG=="NO"]) 
# median(t(chrna1_counts_mg))/median(t(chrna1_counts_nmg)) 
# 
# meta_data$patient_id[meta_data$MG=="YES"] 
# meta_data$patient_id[meta_data$MG=="NO"] 
# 
# 
# 
# ##########################################################################################################################
## pathway analysis

library(clusterProfiler)


# extract ensemble and entrze gene ID mapping files
k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,
                               columns=c("ENTREZID",'SYMBOL'),
                               keytype = 'ENSEMBL')

result_tmg <- res_tmg_t %>% left_join(en2ENSE[,1:2],by=(c('gene_id'='ENSEMBL'))) %>%
  dplyr::select(gene_name,gene_id,ENTREZID,gene_type,log2FoldChange,pvalue,padj)

result_tmg[result_tmg$gene_name=="AIRE",]

result_tmg_all<-result_tmg  %>% filter(padj<0.05)

result_tmg_up<-result_tmg  %>% filter(padj<0.05 & log2FoldChange>0)
result_tmg_down<-result_tmg %>% filter(pvalue<0.05 & log2FoldChange<0)




#### DE number plot

DE_number<-data.frame(Group=c('Up','Down'),
                      DE_gene_number=c(length(result_tmg_up$gene_id),length(result_tmg_down$gene_id)))
p <- ggplot(DE_number, aes(x=Group, y=DE_gene_number,fill=Group)) +
  geom_bar(stat="identity",  colour="black", position=position_dodge(),width = 0.6)+
  theme_classic() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))
p + scale_fill_npg()


## KEGG enrichment analysis
KEGG_up_tmg<-enrichKEGG(result_tmg_up$ENTREZID[1:100], organism = "hsa", keyType = "kegg",
           pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1,
           universe = na.omit(result_tmg$ENTREZID))
barplot(KEGG_up_tmg,showCategory = 28)


KEGGM_down_tmg<-enrichMKEGG(result_tmg_down$ENTREZID[1:100], organism = "hsa", minGSSize=1,
                          pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1,
            universe = na.omit(result_tmg$ENTREZID))
barplot(KEGGM_down_tmg,showCategory = 8)

##  GO analysis
go_up_tmg_bp <- enrichGO(gene          = result_tmg_up$ENTREZID,
                universe      = result_tmg$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "none",
                minGSSize = 3,
                maxGSSize = 150,
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1,
                readable      = T)
view(go_up_tmg_bp@result)


barplot(go_up_tmg_bp,showCategory = 20)
heatplot(go_up_tmg_bp)

p<-dotplot(go_up_tmg_bp,showCategory = 15,font.size=14)

library("reshape2")
p+scale_color_material("red",reverse = T)


View(go_up_tmg_bp@result)

go_down_tmg_bp <- enrichGO(gene          = result_tmg_down$ENTREZID,
                           universe      = result_tmg$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "none",
                         minGSSize = 2,
                         maxGSSize = 150,
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 1,
                         readable      = T)


barplot(go_down_tmg_bp,showCategory = 15)
p<-dotplot(go_down_tmg_bp,showCategory = 15,font.size=14)
p+scale_color_material("indigo",reverse = T)

View(go_down_tmg_bp@result)

write.csv(go_down_tmg_bp@result,"TAMG_vs_nonMG_TCGA_downregulated_GO.csv")


#GSEA analysis

gene_lis<-  result_tmg[,c(3,5)] %>% distinct(ENTREZID,.keep_all = T) %>% na.omit()
gene_list <- gene_lis[,2]
names(gene_list)<-gene_lis[,1]
gene_list = sort(gene_list, decreasing = TRUE)

GO_kk_entrez <- gseGO(geneList     = gene_list,
                      ont          = "ALL", 
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "ENTREZID",
                      pvalueCutoff = 1)   #实际为padj
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb=org.Hs.eg.db,
                           keyType='ENTREZID')#转化id 
View(GO_kk@result)




###提取目标基因的表达量，导入graphpad进行画图

antigen_cpm<-ExpDataCPM %>% filter(GeneSymbol %in% c("CHRNA1","TTN","RYR1","NEFL","NEFM","RYR3","GABRA5","AIRE"))
antigen_cpm2<-antigen_cpm[,c(3:length(colnames(antigen_cpm)))]%>% t() %>% data.frame()
colnames(antigen_cpm2)<-antigen_cpm$GeneSymbol
antigen_cpm2<-antigen_cpm2%>% mutate(patient_id=row.names(.)) %>% left_join(meta_data_RNA[,c(1,6)],by="patient_id")

write.csv(antigen_cpm2,"自身抗原基因mRNA表达量.csv")

# #########
# 
# positive regulation of neuron projection development
#ligand-gated channel activity
# postsynaptic specialization membrane
# presynaptic membrane
# positive regulation of synapse assembly
#cerebellum development
#regulation of postsynaptic membrane potential；
#dendrite membrane；voltage-gated potassium channel complex；
#postsynaptic specialization organization；
#synaptic transmission, GABAergic
#calcium-ion regulated exocytosis
#spinal cord development；hippocampus development
#main axon
#eurofilament bundle assembly

View(go_up_tmg_bp@result)
# 
# ######
# 
# 
## heatmap: axon list related genes
# 
go_up_tmg_bp
library(ComplexHeatmap)
# 
go_up_tmg_bp_figure<-go_up_tmg_bp

go_up_tmg_bp_figure@result<-go_up_tmg_bp@result %>% 
  filter((pvalue<0.05 & (str_detect(.$Description,pattern="neuro" ) |
                           str_detect(.$Description,pattern="synap" )))|
           str_detect(.$Description,pattern="neurofilament" )|
           str_detect(.$Description,pattern="gamma-aminobutyric" )|
           str_detect(.$Description,pattern="axon" ) |
           str_detect(.$Description,pattern="ligand-gated channel activity"))


go_up_tmg_bp_figure@result<-go_up_tmg_bp@result %>% 
  filter((pvalue<0.05 & (str_detect(.$Description,pattern="positive regulation of neuron projection development" ) |
                           str_detect(.$Description,pattern="ligand-gated channel activity" )))|
           str_detect(.$Description,pattern="postsynaptic specialization membrane" )|
           str_detect(.$Description,pattern="presynaptic membrane" )|
           str_detect(.$Description,pattern="positive regulation of synapse assembly" ) |
           str_detect(.$Description,pattern="cerebellum development") |
           str_detect(.$Description,pattern="regulation of postsynaptic membrane potential")|
           str_detect(.$Description,pattern="dendrite membrane")|
           str_detect(.$Description,pattern="voltage-gated potassium channel complex")|
           str_detect(.$Description,pattern="postsynaptic specialization organization")|
           str_detect(.$Description,pattern="synaptic transmission, GABAergic")|
           str_detect(.$Description,pattern="calcium-ion regulated exocytosis")|
           str_detect(.$Description,pattern="spinal cord development")|
           str_detect(.$Description,pattern="hippocampus development")|
           str_detect(.$Description,pattern="main axon")|
           str_detect(.$Description,pattern="neurofilament bundle assembly"))


#ligand-gated channel activity

axon_list<-unique(unlist(str_split(go_up_tmg_bp_figure@result$geneID,pattern = "/")))



axon_t_df<-data.frame(gene_name=axon_list) %>% 
  left_join(res_tmg_t[,c("log2FoldChange","gene_name")],by="gene_name") 


mt<-matrix(data=NA,nrow=length(axon_list),
           ncol=length(go_up_tmg_bp_figure@result$Description)) %>% data.frame()
row.names(mt)<-axon_list
colnames(mt)<-go_up_tmg_bp_figure@result$Description

for (i in c(1:length(go_up_tmg_bp_figure@result$Description))) {
  for (j in c(1:length(axon_list)))
  { 
    if ( axon_list[j] %in% str_split(go_up_tmg_bp_figure@result$geneID[i],pattern = "/")[[1]]) 
    {mt[j,i]<-1}
    else {mt[j,i]<-0}
  }
}


deg2_fig<-deg3
deg2_fig$gene_name<-rownames(deg3)

axon_t_df<-data.frame(gene_name=axon_list) %>% 
  left_join(deg2_fig[,c("avg_log2FC","gene_name")],by="gene_name") %>% 
  left_join(res_tmg_t[,c("log2FoldChange","gene_name")],by="gene_name") %>% 
  rename(logFC_sc_RNA="avg_log2FC",
         logFC_bulk_RNA="log2FoldChange")

#axon_t_df[is.na (axon_t_df)] <- 0

axon_t_df<-cbind(axon_t_df,mt)

ha = columnAnnotation(logFC_sc_RNA = anno_barplot((as.matrix(axon_t_df[,c(2)])),gp = gpar(fill = 3)),
                      logFC_bulk_RNA=anno_barplot((as.matrix(axon_t_df[,c(3)])),gp = gpar(fill = 2)))


Heatmap(t(axon_t_df[,4:length(colnames(axon_t_df))]),cluster_rows = F,cluster_columns = F,
        col=(c("white","#DA4C35")), border = T,rect_gp= gpar(col = "black", lwd = 1),
        row_names_rot = 0,
        top_annotation = ha,
        column_names_rot=60,
        column_names_side="bottom",
        row_names_side="left",
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 6))

ComplexHeatmap::Heatmap(as.matrix(axon_t_df[,3:length(colnames(axon_t_df))]))


Heatmap(t(axon_t_df[,3:length(colnames(axon_t_df))]),col=(c("black","#DA4C35")), border = T,rect_gp= gpar(col = "white", lwd = 2))
# 
# 

