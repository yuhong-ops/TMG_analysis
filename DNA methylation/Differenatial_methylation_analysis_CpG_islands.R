### CPG_tcga_gene

### 选出与唯一基因prmoter 相关的CpG island, 然后进行多重p值矫正
CPG_tcga_gene_is<-CPG_tcga_gene %>% filter (Islands_Name %in% me_num_pisl_id)

CPG_tcga_gene %>% filter(Islands_Name == "chr6:31695894-31698245")

CPG_tcga_gene_is$adj.P.Val<-p.adjust(CPG_tcga_gene_is$P.Value,method = "fdr")

CPG_tcga_gene_is <- CPG_tcga_gene_is %>% dplyr::arrange(P.Value) %>% 
  left_join(en2ENSE,by=(c('V2'='SYMBOL')),multiple="first")
length(CPG_tcga_gene_is$Islands_Name)

write.csv(CPG_tcga_gene_is,"CPG_tcga_gene_is2.csv")

CPG_tcga_gene_is  

### 选出有significant的与gene promoter相关的CpG
CPG_tcga_gene_is_sig<-CPG_tcga_gene_is %>% filter (adj.P.Val<0.05)

## 画图，基因prmoter 相关的CpG island在TAMG上调和下调的柱状图

##write.csv(CPG_tcga_sig,"CPG_island_TCGA.csv")

write.csv(CPG_tcga_sig,"CPG_island_TCGA.csv")
cpg_tcga_down_pro<-CPG_tcga_gene_is_sig[CPG_tcga_gene_is_sig$logFC <0, ]
length(cpg_tcga_down_pro$Islands_Name)
cpg_tcga_up_pro<-CPG_tcga_gene_is_sig[CPG_tcga_gene_is_sig$logFC >0, ]
length(cpg_tcga_up_pro$Islands_Name)
DE_number<-data.frame(Group=c('Down','Up'),
                      DE_met_number=c(length(cpg_tcga_down_pro$logFC),length(cpg_tcga_up_pro$logFC)))
p <- ggplot(DE_number, aes(x=Group, y=DE_met_number,fill=Group)) +
  geom_bar(stat="identity",  colour="black", position=position_dodge(),width = 0.6)+
  theme_classic() + 
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))
p
##
##
CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "MTHFR"))
CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "AIRE"))
CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR1"))

CPG_gse %>% filter(str_detect(.$Islands_Name, pattern= "chr21:45705428-45706044"))
CPG_gse %>% filter(str_detect(.$Islands_Name, pattern= "chr21:45713509-45713813"))

CPG_tcga_gene_is %>% filter(str_detect(.$V2, pattern= "NEFM"))
CPG_tcga_gene_is %>% filter(str_detect(.$V2, pattern= "RYR3"))
CPG_tcga_gene_is %>% filter(str_detect(.$V2, pattern= "RYR1"))
CPG_tcga_gene_is %>% filter(str_detect(.$V2, pattern= "CHRNA1"))
CPG_tcga_gene_is %>% filter(str_detect(.$V2, pattern= "AIRE"))
CPG_tcga_gene_is %>% filter(str_detect(.$V2, pattern= "MTHFR"))

cpg_tmg_cc_down <- enrichGO(gene          = CPG_tcga_gene_is_sig$ENTREZID[CPG_tcga_gene_is_sig$logFC<0],
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
cpg_tmg_cc_down_neu<-cpg_tmg_cc_down

cpg_tmg_cc_down_neu@result<-cpg_tmg_cc_down@result %>%
  filter(str_detect(cpg_tmg_cc_down@result$Description,pattern = "syna") |
           str_detect(cpg_tmg_cc_down@result$Description,pattern = "neu")|
           str_detect(cpg_tmg_cc_down@result$Description,pattern = "nerv")|
           str_detect(cpg_tmg_cc_down@result$Description,pattern = "sarc")|
           str_detect(cpg_tmg_cc_down@result$Description,pattern = "axon")|
           str_detect(cpg_tmg_cc_down@result$Description,pattern = "skeletal")) %>%
  filter(Count>1 | Description=="postsynaptic cytoskeleton" )

barplot(cpg_tmg_cc_down_neu,showCategory = 12)

cpg_tmg_cc_up <- enrichGO(gene          = CPG_tcga_gene_is_sig$ENTREZID[CPG_tcga_gene_is_sig$logFC>0][1:300],
                            universe      = result_tmg$ENTREZID,
                            OrgDb         = org.Hs.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "none",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 1,
                            minGSSize = 2,
                            maxGSSize = 150,
                            readable      = TRUE)
View(cpg_tmg_cc_up@result)

barplot(cpg_tmg_cc_down)

cpg_tmg_cc_up_lymph<-cpg_tmg_cc_up

cpg_tmg_cc_up_lymph@result<-cpg_tmg_cc_up@result %>%
  filter(str_detect(cpg_tmg_cc_up@result$Description,pattern = "T cell") |
           str_detect(cpg_tmg_cc_up@result$Description,pattern = "DNA binding")|
           str_detect(cpg_tmg_cc_up@result$Description,pattern = "thymic")|
           str_detect(cpg_tmg_cc_up@result$Description,pattern = "methyl")|
           str_detect(cpg_tmg_cc_up@result$Description,pattern = "promoter")|
           str_detect(cpg_tmg_cc_up@result$Description,pattern = "lymph")) %>%
  filter(Count>1 | Description=="postsynaptic cytoskeleton" )

barplot(cpg_tmg_cc_up_lymph,showCategory = 12)



####



### RNA DE analysis和甲基化同时有差异的基因及数值
library(disprose)
CPG_tcga_gene_is_sig

#res_tmg_t_sig<-res_tmg_t_fdr %>% filter(padj<0.05)
res_tmg_t_sig
inter<-intersect(c(cpg_tcga_down_pro$V2,cpg_tcga_up_pro$V2),res_tmg_t_sig$gene_name)
cpg_inter<-CPG_tcga_gene_p_order %>% filter(V2 %in% inter) %>% dplyr::select( Islands_Name,logFC,UCSC_RefGene_Name=V2 )

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
             color="#3C5488FF",aes(x=logFC,y=log2FoldChange),size=3)+
  geom_text_repel(data=subset(inter_df, ((logFC<0 & log2FoldChange>0)|(logFC>0 & log2FoldChange<0))),
                  aes(x=logFC,y=log2FoldChange,label=UCSC_RefGene_Name))

## 差异表达的mRNA和差异甲基化的CpG island呈相反趋势的基因列表
cpg_down_exp_up<-inter_df %>% filter(logFC<0,log2FoldChange>0)%>% pull(UCSC_RefGene_Name) %>% unique 
cpg_up_exp_down<-inter_df %>% filter(logFC>0,log2FoldChange<0) %>% pull(UCSC_RefGene_Name) %>% unique

unique(inter_df$UCSC_RefGene_Name)

cpg_down_exp_down<-inter_df %>% filter(logFC<0,log2FoldChange<0)%>% pull(UCSC_RefGene_Name) %>% unique 
cpg_up_exp_up<-inter_df %>% filter(logFC>0,log2FoldChange>0) %>% pull(UCSC_RefGene_Name) %>% unique


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

cpg_down_geid_up_go_select<-cpg_down_geid_up_go
cpg_down_geid_up_go_select@result<-cpg_down_geid_up_go@result %>% filter(ID %in% c("GO:0005743","GO:1902495","GO:0015276","GO:0099572","GO:0043235","GO:0034704",
                                                                                   "GO:0033693","GO:0019722","GO:0016528","GO:0021953","GO:0010810","GO:0070233",
                                                                                   "GO:0031594","GO:0032281","GO:0016917")) %>% 
  arrange(p.adjust)

barplot(cpg_down_geid_up_go_select,showCategory = 16)


cpg_up_geid_down_go <- enrichGO(gene          = cpg_up_geid_down,
                                universe      = result_tmg$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "ALL",
                                pAdjustMethod = "none",
                                minGSSize = 1,
                                pvalueCutoff  = 0.05,
                                readable      = TRUE)
View(cpg_up_geid_down_go@result)

####

#####
inter_df

# ## scatter plot of gse and tcga CPG island fold change for tcga CPG island significant genes
CPG_tcga_single_gene_sig_add_gse<-CPG_tcga_gene_is_sig %>% dplyr::select(Islands_Name,tcga_logFC=logFC,Gene_name=V2 )%>%
  left_join(CPG_gse[,c("Islands_Name","logFC")],by="Islands_Name") %>% dplyr::rename(gse_logFC=logFC)
CPG_tcga_gene_is_sig %>% filter(V2 %in% "NEFM")

chr8:24770908-24772547
CPG_gse_sig$Islands_Name %>% filter(Islands_Name %in% "chr8:24770908-24772547")

CPG_tcga_single_gene_sig_add_gse_mRNA_sig<-cpg_inter %>% left_join(CPG_tcga_single_gene_sig_add_gse,by="Islands_Name")
# 
ggplot(data = CPG_tcga_single_gene_sig_add_gse) + 
  geom_point(aes(x=gse_logFC,y=tcga_logFC),color="black",size=0.9) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  theme_classic(base_size=16) + 
  geom_smooth(aes(x=gse_logFC,y=tcga_logFC),method = "lm") 
summary(lm(CPG_tcga_single_gene_sig_add_gse$tcga_logFC~CPG_tcga_single_gene_sig_add_gse$gse_logFC))

ggplot(data = CPG_tcga_single_gene_sig_add_gse_mRNA_sig) + 
  geom_point(aes(x=gse_logFC,y=tcga_logFC),color="black",size=0.9) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  theme_classic(base_size=16) + 
  geom_smooth(aes(x=gse_logFC,y=tcga_logFC),method = "lm") +
  geom_text_repel(data=subset(CPG_tcga_single_gene_sig_add_gse_mRNA_sig, UCSC_RefGene_Name %in% c("NEFM","RYR3","GABRA5")), 
                  aes(x=gse_logFC,y=tcga_logFC,label=UCSC_RefGene_Name))+  
  geom_point(data=subset(CPG_tcga_single_gene_sig_add_gse_mRNA_sig, UCSC_RefGene_Name %in% c("NEFM","RYR3","GABRA5")),
             color="red",aes(x=gse_logFC,y=tcga_logFC),size=3) 

summary(lm(CPG_tcga_single_gene_sig_add_gse_mRNA_sig$tcga_logFC~CPG_tcga_single_gene_sig_add_gse_mRNA_sig$gse_logFC))

