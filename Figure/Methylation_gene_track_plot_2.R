library(Gviz)

# gse dataset NEFM plot
# indicate which genome is being used
gen <- "hg19"

#### RYR3

# the index of the DMR that we will plot 
which(results.ranges$overlapping.genes=="RYR3") # 5488
dmrIndex <- 6915
results.ranges[which(results.ranges$overlapping.genes=="RYR3"),]

# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))


# CpG islands (download from Wu etal  http://www.haowulab.org/software/makeCGI/index.html )
islandHMM <- read.csv("reference_data/model-based-cpg-islands-hg19.txt",
                      sep="\t", stringsAsFactors=FALSE, header=T)


islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))

## ideogram track and genome axis track
iTrack <- IdeogramTrack(genome = gen, chromosome = "chr15", name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")


##build gene region track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
rTrack <- GeneRegionTrack(txdb, chromosome = "chr15", 
                          start = 33602816,  end = 34122816)

##Ensure that the methylation data is ordered by chromosome and base position.


CPG_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))

myDMP_tcga_RYR3<-myDMP1 %>% dplyr::filter(gene == "RYR3") %>% 
  dplyr:: filter ((MAPINFO < 33604003) & (MAPINFO > 33602816))

myDMP_gse_RYR3<-myDMP_gse$no_to_yes %>% dplyr::filter(gene == "RYR3") %>% 
  dplyr:: filter ((MAPINFO < 33604003) & (MAPINFO > 33602816))

myDMP_gse_RYR3 <- myDMP_gse_RYR3[order(myDMP_gse_RYR3$CHR,myDMP_gse_RYR3$MAPINFO),]


myDMP_tcga_RYR3 <- myDMP_tcga_RYR3[order(myDMP_tcga_RYR3$CHR,myDMP_tcga_RYR3$MAPINFO),]

#  methylation data
RYR3_beta_gse <- myNorm_met_gse[match(rownames(myDMP_gse_RYR3),rownames(myNorm_met_gse)),]
RYR3_beta_tcga <- myNorm_met[match(rownames(myDMP_tcga_RYR3),rownames(myNorm_met)),]


# create genomic ranges object from methylaton data
cpgData <- GRanges(seqnames=Rle(rep("chr15",7)),  ##改一下
                   ranges=IRanges(start=myDMP_gse_RYR3$MAPINFO, end=myDMP_gse_RYR3$MAPINFO),
                   strand=Rle(rep("*",nrow(myDMP_gse_RYR3))),
                   betas=RYR3_beta_gse)


cpgData_tcga <- GRanges(seqnames=Rle(rep("chr15",8)),  ##改一下
                        ranges=IRanges(start=myDMP_tcga_RYR3$MAPINFO, end=myDMP_tcga_RYR3$MAPINFO),
                        strand=Rle(rep("*",nrow(myDMP_tcga_RYR3))),
                        betas=RYR3_beta_tcga)

# extract data on CpGs in DMR
#cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])


# methylation data track
methTrack_gse <- DataTrack(range=cpgData, groups=meta_data_thymoma_kajiura$Myasthenia_gravis,genome = gen,
                           chromosome="chr15", ylim=c(0,0.6), col=pal,
                           type=c("a","p"), name="DNA Meth.\n(beta value)",
                           background.panel="white", legend=TRUE, cex.title=0.8,
                           cex.axis=0.8, cex.legend=0.8)

methTrack_tcga <- DataTrack(range=cpgData_tcga, groups=met_meta$MG,genome = gen,
                            chromosome="chr15", ylim=c(0,0.8), col=pal,
                            type=c("a","p"), name="DNA Meth.\n(beta value)",
                            background.panel="white", legend=TRUE, cex.title=0.8,
                            cex.axis=0.8, cex.legend=0.8)

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome="chr15",fill="darkgreen")


#Set up the track list and indicate the relative sizes of the different tracks. Finally, draw the plot using the plotTracks function (Figure 11).

tracks <- list(iTrack, gTrack, rTrack, islandTrack, methTrack_tcga, methTrack_gse)
tracks <- list(iTrack, gTrack, rTrack, islandTrack, methTrack_tcga)

sizes <- c(2,2,4,2,12,12) # set up the relative sizes of the tracks
plotTracks(tracks, from=33602700, to=34122816, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))

plotTracks(tracks, from=33602700, to=33607700, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))


## GABRA5
# the index of the DMR that we will plot 
which(results.ranges$overlapping.genes=="GABRA5") # 5488
dmrIndex <- 2719
results.ranges[which(results.ranges$overlapping.genes=="GABRA5"),]

# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <-  27184000 #as.numeric(start(results.ranges[dmrIndex]))
end <-  27189000 #as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))


# CpG islands (download from Wu etal  http://www.haowulab.org/software/makeCGI/index.html )
islandHMM <- read.csv("reference_data/model-based-cpg-islands-hg19.txt",
                      sep="\t", stringsAsFactors=FALSE, header=T)
head(islandHMM)


islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData

## ideogram track and genome axis track
iTrack <- IdeogramTrack(genome = gen, chromosome = "chr15", name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")


##build gene region track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
rTrack <- GeneRegionTrack(txdb, chromosome = "chr15", 
                          start = 27184000,  end = 27189000)

##Ensure that the methylation data is ordered by chromosome and base position.

myDMP_tcga_GABRA5<-myDMP1 %>% dplyr::filter(gene == "GABRA5") %>% 
  dplyr:: filter ((MAPINFO < 27189000) & (MAPINFO > 27184000))

myDMP_gse_GABRA5<-myDMP_gse$no_to_yes %>% dplyr::filter(gene == "GABRA5") %>% 
  dplyr:: filter ((MAPINFO < 27189000) & (MAPINFO > 27184000))

myDMP_gse_GABRA5 <- myDMP_gse_GABRA5[order(myDMP_gse_GABRA5$CHR,myDMP_gse_GABRA5$MAPINFO),]


myDMP_tcga_GABRA5 <- myDMP_tcga_GABRA5[order(myDMP_tcga_GABRA5$CHR,myDMP_tcga_GABRA5$MAPINFO),]

#  methylation data
GABRA5_beta_gse <- myNorm_met_gse[match(rownames(myDMP_gse_GABRA5),rownames(myNorm_met_gse)),]
GABRA5_beta_tcga <- myNorm_met[match(rownames(myDMP_tcga_GABRA5),rownames(myNorm_met)),]


# create genomic ranges object from methylaton data
cpgData <- GRanges(seqnames=Rle(rep("chr15",4)),  ##改一下
                   ranges=IRanges(start=myDMP_gse_GABRA5$MAPINFO, end=myDMP_gse_GABRA5$MAPINFO),
                   strand=Rle(rep("*",nrow(myDMP_gse_GABRA5))),
                   betas=GABRA5_beta_gse)


cpgData_tcga <- GRanges(seqnames=Rle(rep("chr15",4)),  ##改一下
                        ranges=IRanges(start=myDMP_tcga_GABRA5$MAPINFO, end=myDMP_tcga_GABRA5$MAPINFO),
                        strand=Rle(rep("*",nrow(myDMP_tcga_GABRA5))),
                        betas=GABRA5_beta_tcga)

# extract data on CpGs in DMR
#cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])


# methylation data track
methTrack_gse <- DataTrack(range=cpgData, groups=meta_data_thymoma_kajiura$Myasthenia_gravis,genome = gen,
                           chromosome="chr15", ylim=c(0,0.6), col=pal,
                           type=c("a","p"), name="DNA Meth.\n(beta value)",
                           background.panel="white", legend=TRUE, cex.title=0.8,
                           cex.axis=0.8, cex.legend=0.8)

methTrack_tcga <- DataTrack(range=cpgData_tcga, groups=met_meta$MG,genome = gen,
                            chromosome="chr15", ylim=c(0,0.8), col=pal,
                            type=c("a","p"), name="DNA Meth.\n(beta value)",
                            background.panel="white", legend=TRUE, cex.title=0.8,
                            cex.axis=0.8, cex.legend=0.8)

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome="chr15",fill="darkgreen")


#Set up the track list and indicate the relative sizes of the different tracks. Finally, draw the plot using the plotTracks function (Figure 11).

tracks <- list(iTrack, gTrack, rTrack, islandTrack, methTrack_tcga, methTrack_gse)
tracks <- list(iTrack, gTrack, rTrack, islandTrack, methTrack_tcga)

sizes <- c(2,2,4,2,12,12) # set up the relative sizes of the tracks
plotTracks(tracks, from=27184000, to=27189000, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))








