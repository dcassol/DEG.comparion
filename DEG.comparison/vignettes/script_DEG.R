##-------------------------------------------------------------------
##   RNA-Seq - Comparison of DEG Methods  
##   Last update: June 2015
##   Daniela Cassol
##   Project_RNA-seq1: data set from (Jiao & Meyerowitz 2010)
##-------------------------------------------------------------------
## R code from vignette source 'script_DEG.Rnw'


#############################################
###DEG.comparison
#############################################
## library("DEG.comparison") # Loads the package
## library(help="DEG.comparison") # Lists all functions and classes 
## vignette("DEG.comparison") # Opens this PDF manual from R


##-------------------------------------------------------------------
##   Workflow: systemPipeR  
##-------------------------------------------------------------------

###############################################
###Workflow: systemPipeR - Structure of targets file: Single end (SE)
###############################################
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment.char = "#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")


# ###############################################
# ###Workflow: systemPipeR - Structure of param file and SYSargs container
# ###############################################
# parampath <- system.file("extdata", "tophat.param", package="DEG.comparison")
# read.delim(parampath, comment.char = "#")
# args <- systemArgs(sysma=parampath, mytargets=targetspath)
# args
# names(args)
# modules(args)
# cores(args)
# outpaths(args)[1]
# sysargs(args)[1]
# systemArgs(sysma=parampath, mytargets=targetspath, type="json")
# 
# 
# ###############################################
# ###Workflow: systemPipeR - Define environment settings and samples
# ###############################################
# trim.param <- system.file("extdata", "trim.param", package="DEG.comparison")
# trim.param <- read.delim(trim.param, comment.char = "#")
# args <- systemArgs(sysma="trim.param", mytargets=targetspath) #Construct SYSargs object from param and targets files.
# 
# 
# ###############################################
# ###Workflow: systemPipeR - Read Preprocessing
# ###############################################
# preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
# writeTargetsout(x=args, file="targets_trim.txt")
# 
# filterFct <- function(fq) {
#   filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
#   filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes low complexity reads
#   filter <- compose(filter1, filter2)
#   fq[filter(fq)]
# }
# preprocessReads(args=args, Fct="filterFct(fq)", batchsize=100000)
# 
# 
# ###############################################
# ###Workflow: systemPipeR - FASTQ quality report
# ###############################################
# fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
# pdf("./results/fastqReport3.pdf", height=18, width=4*length(fqlist))
# seeFastqPlot(fqlist)
# dev.off()
# 
# 
# ###############################################
# ###Workflow: systemPipeR - Alignment with Tophat2
# ###############################################
# args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
# moduleload(modules(args)) # Skip if module system is not available
# system("bowtie2-build ./data/TAIR10_chr_all.fas ./data/TAIR10_chr_all.fas")
# bampaths <- runCommandline(args=args)
# file.copy((paste0("./.BatchJobs.R"), ".")
# file.copy(paste0("./torque.tmpl"), ".")
# resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="10gb")
# reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", resourceList=resources)
# showStatus(reg)
# file.exists(outpaths(args))
# sapply(1:length(args), function(x) loadResult(reg, x)) # Works after job completion
#           
#           
# ###############################################
# ###Workflow: systemPipeR - Read and alignment count stats
# ###############################################
# read_statsDF <- alignStats(args)
# write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
read.table(system.file("extdata", "alignStats.xls", package="DEG.comparison"), header=TRUE)[1:4,]
#           
#           
# ###############################################
# ###Workflow: systemPipeR - Create symbolic links for viewing BAM files in IGV
# ###############################################
# symLink2bam(sysargs=args, htmldir=c("~/.html/", "cassol_GEN242/"),urlbase="http://biocluster.ucr.edu/~dcassol/", urlfile="IGVurl.txt")
#         
#           
# ###############################################
# ###Workflow: systemPipeR - Read counting for mRNA profiling experiments
# ###############################################
# ###Read counting for mRNA proCreate txdb (needs to be done only once)
# library(GenomicFeatures); library(BiocParallel)
# txdb <- makeTranscriptDbFromGFF(file="data/TAIR10_GFF3_genes.gff", format="gff3", dataSource="TAIR", species="Arabidopsis thaliana")
# saveDb(txdb, file="./data/tair10.sqlite")
#           
# ###Read counting with summarizeOverlaps in parallel mode with multiple cores
# txdb <- loadDb("./data/tair10.sqlite")
# eByg <- exonsBy(txdb, by="gene")
# bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
# multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
# # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE' 
# counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE))
# countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
# rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
# countDFeByg[1:4,1:12]
# write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
# rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
# write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
# rpkmDFeByg[1:4,1:7]
# 
# 
# ###############################################  
# ###Workflow: systemPipeR - Correlation analysis of samples / cluster and create phylogenetic tree
# ###############################################  
# library(ape)
# rpkmDFeBygpath <-paste0("./results/rpkmDFeByg.xls")
# rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
# rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
# d <- cor(rpkmDFeByg, method="spearman")
# hc <- hclust(as.dist(1-d))
# pdf("results/sample_tree.pdf")
# plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
# dev.off()

##-------------------------------------------------------------------
##   Workflow: DEG.comparison
##-------------------------------------------------------------------

##Data
countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="DEG.comparison")
countDFeByg <- read.delim(countDFeBygpath, row.names=1)
rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="DEG.comparison")
rpkmDFeByg <- read.delim(rpkmDFeBygpath, row.names=1)

##Comparasions
Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"), AP3.4_AP3.67=c("AP3.4A","AP3.4B", "AP3.67A" ,"AP3.67B"), 
              AG.4_AG.67=c("AG.4A", "AG.4B", "AG.67A","AG.67B"), AP1.4_AP3.4=c("AP1.4A","AP1.4B", "AP3.4A","AP3.4B"), 
              AP1.4_AG.4=c("AP1.4A", "AP1.4B", "AG.4A", "AG.4B"), AP3.4_AG.4=c("AP3.4A", "AP3.4B", "AG.4A", "AG.4B"), 
              AP1.67_AP3.67=c("AP1.67A","AP1.67B", "AP3.67A", "AP3.67B"), AP1.67_AG.67=c("AP1.67A","AP1.67B", "AG.67A", "AG.67B"), 
              AP3.67_AG.67=c("AP3.67A", "AP3.67B", "AG.67A", "AG.67B"))


################################################################  
###DEG1:  Simple Fold Change Method - RPKM  
################################################################ 
##Settings
Comp1 <- list(Factor=(Reduce(union, targets$Factor)), Sample=c(colnames(rpkmDFeByg)), 
              group=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))

Comp2 <- list(AP1.4_AP1.67=c("AP1.4", "AP1.67"), AP3.4_AP3.67=c("AP3.4", "AP3.67"), 
              AG.4_AG.67=c("AG.4", "AG.67"), AP1.4_AP3.4=c("AP1.4", "AP3.4"), 
              AP1.4_AG.4=c("AP1.4", "AG.4"), AP3.4_AG.4=c("AP3.4", "AG.4"), 
              AP1.67_AP3.67=c("AP1.67", "AP3.67"), AP1.67_AG.67=c("AP1.67", "AG.67"), 
              AP3.67_AG.67=c("AP3.67", "AG.67"))

##Compute mean values for replicates and logFC for comparisons
RPKM_FC <- run_RPKM (rpkmDFeByg, Comp1, Comp2)
write.table(RPKM_FC, "./results/RPKM_FC.xls", quote=FALSE, sep="\t", col.names = NA)
pdf("./results/DEG_list_RPKM.pdf") 
DEG_list_RPKM <- filterDEG_logFC(degDF=RPKM_FC, filter=c(Fold=2), method="RPKM")
dev.off()
DEG_list_RPKM$Summary[1:4,]


################################################################
###DEG2: edgeR  
################################################################
edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")
write.table(edgeDF, "results/edgeDF.xls", col.names=NA, quote=FALSE, sep="\t")
pdf("./results/DEG_list_edgeR.pdf") 
DEG_list_edgeR <- filterDEGnew(degDF=edgeDF, filter=c(Fold=2, FDR=1), method="edgeR")
dev.off()
DEG_list_edgeR$Summary[1:4,]


################################################################
###DEG3: DESeq2
################################################################
deseq2DF <- run_DESeq2(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE)
write.table(deseq2DF, "results/deseq2DF.xls", col.names=NA, quote=FALSE, sep="\t")
pdf("./results/DEG_list_DESeq2.pdf")
DEG_list_DESeq2 <- filterDEGnew(degDF=deseq2DF, filter=c(Fold=2, FDR=1), method="DESeq2")
dev.off()
DEG_list_DESeq2$Summary[1:4,]


################################################################
###DEG4: baySeq
################################################################
dim(countDFeByg)
bayseqDF <- run_BaySeq(countDFeByg, Comp3, number=27416)
write.table(bayseqDF, "results/bayseqDF.xls", col.names=NA, quote=FALSE, sep="\t")
pdf("./results/DEG_list_bayseqDF.pdf") 
DEG_list_bayseqDF <- filterDEG_FDR(degDF=bayseqDF, filter=c(FDR=1), method="BaySeq")
dev.off()
DEG_list_bayseqDF$Summary[1:4,]
       

################################################################
###DEG5: NBPSeq
################################################################
###NBPSeq.glm
#For each row of the input data matrix, nb.glm.test fits an NB log-linear regression model and performs large-sample tests for a one-dimensional regression coefficient.

NBPSeq.glmDF <- run_NBPSeq_glm (countDFeByg, Comp3)
write.table(NBPSeq.glmDF, "results/NBPSeq_glmDF.xls", col.names=NA, quote=FALSE, sep="\t")
pdf("./results/DEG_list_NBPSeq_glmDF.pdf")
DEG_list_NBPSeq.glmDF <- filterDEGnew(degDF=NBPSeq.glmDF, filter=c(Fold=2, FDR=1), method="NBPSeq.glm")
dev.off()
DEG_list_NBPSeq.glmDF$Summary[1:4,]

###NBPSeq.nbp.test
#nbp.test fits an NBP model to the RNA-Seq counts and performs Robinson and Smyth's exact NB test on each gene to assess differential gene expression between two groups.

NBPSeq.nbpDF <- run_NBPSeq_nbp (countDFeByg, Comp3)
write.table(NBPSeq.nbpDF, "results/NBPSeq.nbpDF.xls", col.names=NA, quote=FALSE, sep="\t")
pdf("./results/DEG_list_NBPSeq.nbpDF.pdf")
DEG_list_NBPSeq.nbpDF <- filterDEGnew(degDF=NBPSeq.nbpDF, filter=c(Fold=2, FDR=1), method="NBPSeq.nbp")
dev.off()
DEG_list_NBPSeq.nbpDF$Summary[1:4,]


################################################################
###DEG6: TSPM 
################################################################
TSPMDF <- run_TSPM(countDFeByg, Comp3)
write.table(TSPMDF, "results/TSPMDF.xls", col.names=NA, quote=FALSE, sep="\t")
pdf("./results/DEG_list_TSPM.pdf")
DEG_list_TSPM <- filterDEGnew(degDF=TSPMDF, filter=c(Fold=2, FDR=1), method="TSPM")
dev.off()
DEG_list_TSPM$Summary[1:4,]


################################################################
###DEG: Comparison
################################################################

################################################################ 
###Number of significantly differentially expressed 
################################################################ 
totalgenes <- system.file("extdata", "totalgenes.csv", package="DEG.comparison")
totalgenes <- read.delim (totalgenes, sep=",")
pdf("./results/TotalGenes.pdf")
plot <- ggplot(totalgenes, aes(Package, Number.of.significantly.differentially.expressed)) + geom_bar(aes(fill = DEG), stat="identity") + facet_wrap(~DEG, ncol=1)
print(plot) 
dev.off()


################################################################ 
###data.frame and list with all DEGs 
################################################################ 
      ###Up and Down
      aa1<- (DEG_list_RPKM$UporDown$"AP1.4_AP1.67"); aa2<- (DEG_list_RPKM$UporDown$"AP3.4_AP3.67"); aa3<- (DEG_list_RPKM$UporDown$"AG.4_AG.67"); aa4<- (DEG_list_RPKM$UporDown$"AP1.4_AP3.4"); aa5<- (DEG_list_RPKM$UporDown$"AP1.4_AG.4"); aa6<- (DEG_list_RPKM$UporDown$"AP3.4_AG.4"); aa7<- (DEG_list_RPKM$UporDown$"AP1.67_AP3.67"); aa8<- (DEG_list_RPKM$UporDown$"AP1.67_AG.67"); aa9<- (DEG_list_RPKM$UporDown$"AP3.67_AG.67"); 
      DEG_list_RPKM_UporDown <-Reduce(union, list(aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9))
      dfDEG_RPKM_UporDown <- data.frame(DEG_list_RPKM_UporDown); rownames(dfDEG_RPKM_UporDown) <- dfDEG_RPKM_UporDown[[1]]
      #
      bb1<- (DEG_list_edgeR$UporDown$"AP1.4.AP1.67"); bb2<- (DEG_list_edgeR$UporDown$"AP3.4.AP3.67"); bb3<- (DEG_list_edgeR$UporDown$"AG.4.AG.67"); bb4<- (DEG_list_edgeR$UporDown$"AP1.4.AP3.4"); bb5<- (DEG_list_edgeR$UporDown$"AP1.4.AG.4"); bb6<- (DEG_list_edgeR$UporDown$"AP3.4.AG.4"); bb7<- (DEG_list_edgeR$UporDown$"AP1.67.AP3.67"); bb8<- (DEG_list_edgeR$UporDown$"AP1.67.AG.67"); bb9<- (DEG_list_edgeR$UporDown$"AP3.67.AG.67"); 
      DEG_list_edgeR_UporDown <-Reduce(union, list(bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb8,bb9))
      dfDEG_edgeR_UporDown <- data.frame(DEG_list_edgeR_UporDown); rownames(dfDEG_edgeR_UporDown) <- dfDEG_edgeR_UporDown[[1]]
      #
      cc1<- (DEG_list_DESeq2$UporDown$"AP1.4.AP1.67"); cc2<- (DEG_list_DESeq2$UporDown$"AP3.4.AP3.67"); cc3<- (DEG_list_DESeq2$UporDown$"AG.4.AG.67"); cc4<- (DEG_list_DESeq2$UporDown$"AP1.4.AP3.4"); cc5<- (DEG_list_DESeq2$UporDown$"AP1.4.AG.4"); cc6<- (DEG_list_DESeq2$UporDown$"AP3.4.AG.4"); cc7<- (DEG_list_DESeq2$UporDown$"AP1.67.AP3.67"); cc8<- (DEG_list_DESeq2$UporDown$"AP1.67.AG.67"); cc9<- (DEG_list_DESeq2$UporDown$"AP3.67.AG.67"); 
      DEG_list_DESeq2_UporDown <-Reduce(union, list(cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9))
      dfDEG_DESeq2_UporDown <- data.frame(DEG_list_DESeq2_UporDown); rownames(dfDEG_DESeq2_UporDown) <- dfDEG_DESeq2_UporDown[[1]]
      #
      dd1<- (DEG_list_bayseqDF$UporDown$"AP1.4_AP1.67"); dd2<- (DEG_list_bayseqDF$UporDown$"AP3.4_AP3.67"); dd3<- (DEG_list_bayseqDF$UporDown$"AG.4_AG.67"); dd4<- (DEG_list_bayseqDF$UporDown$"AP1.4_AP3.4"); dd5<- (DEG_list_bayseqDF$UporDown$"AP1.4_AG.4"); dd6<- (DEG_list_bayseqDF$UporDown$"AP3.4_AG.4"); dd7<- (DEG_list_bayseqDF$UporDown$"AP1.67_AP3.67"); dd8<- (DEG_list_bayseqDF$UporDown$"AP1.67_AG.67"); dd9<- (DEG_list_bayseqDF$UporDown$"AP3.67_AG.67"); 
      DEG_list_bayseqDF_UporDown <-Reduce(union, list(dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9))
      dfDEG_bayseqDF_UporDown <- data.frame(DEG_list_bayseqDF_UporDown); rownames(dfDEG_bayseqDF_UporDown) <- dfDEG_bayseqDF_UporDown[[1]]
      #
      ee1<- (DEG_list_NBPSeq.glmDF$UporDown$"AP1.4_AP1.67"); ee2<- (DEG_list_NBPSeq.glmDF$UporDown$"AP3.4_AP3.67"); ee3<- (DEG_list_NBPSeq.glmDF$UporDown$"AG.4_AG.67"); ee4<- (DEG_list_NBPSeq.glmDF$UporDown$"AP1.4_AP3.4"); ee5<- (DEG_list_NBPSeq.glmDF$UporDown$"AP1.4_AG.4"); ee6<- (DEG_list_NBPSeq.glmDF$UporDown$"AP3.4_AG.4"); ee7<- (DEG_list_NBPSeq.glmDF$UporDown$"AP1.67_AP3.67"); ee8<- (DEG_list_NBPSeq.glmDF$UporDown$"AP1.67_AG.67"); ee9<- (DEG_list_NBPSeq.glmDF$UporDown$"AP3.67_AG.67"); 
      DEG_list_NBPSeq.glmDF_UporDown <-Reduce(union, list(ee1,ee2,ee3,ee4,ee5,ee6,ee7,ee8,ee9))
      dfDEG_NBPSeq.glmDF_UporDown <- data.frame(DEG_list_NBPSeq.glmDF_UporDown); rownames(dfDEG_NBPSeq.glmDF_UporDown) <- dfDEG_NBPSeq.glmDF_UporDown[[1]]
      #
      ff1<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP1.4_AP1.67"); ff2<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP3.4_AP3.67"); ff3<- (DEG_list_NBPSeq.nbpDF$UporDown$"AG.4_AG.67"); ff4<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP1.4_AP3.4"); ff5<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP1.4_AG.4"); ff6<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP3.4_AG.4"); ff7<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP1.67_AP3.67"); ff8<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP1.67_AG.67"); ff9<- (DEG_list_NBPSeq.nbpDF$UporDown$"AP3.67_AG.67"); 
      DEG_list_NBPSeq.nbpDF_UporDown <-Reduce(union, list(ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9))
      dfDEG_NBPSeq.nbpDF_UporDown <- data.frame(DEG_list_NBPSeq.nbpDF_UporDown); rownames(dfDEG_NBPSeq.nbpDF_UporDown) <- dfDEG_NBPSeq.nbpDF_UporDown[[1]]
      #
      gg1<- (DEG_list_TSPM$UporDown$"AP1.4_AP1.67"); gg2<- (DEG_list_TSPM$UporDown$"AP3.4_AP3.67"); gg3<- (DEG_list_TSPM$UporDown$"AG.4_AG.67"); gg4<- (DEG_list_TSPM$UporDown$"AP1.4_AP3.4"); gg5<- (DEG_list_TSPM$UporDown$"AP1.4_AG.4"); gg6<- (DEG_list_TSPM$UporDown$"AP3.4_AG.4"); gg7<- (DEG_list_TSPM$UporDown$"AP1.67_AP3.67"); gg8<- (DEG_list_TSPM$UporDown$"AP1.67_AG.67"); gg9<- (DEG_list_TSPM$UporDown$"AP3.67_AG.67"); 
      DEG_list_TSPM_UporDown <-Reduce(union, list(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9))
      dfDEG_TSPM_UporDown <- data.frame(DEG_list_TSPM_UporDown); rownames(dfDEG_TSPM_UporDown) <- dfDEG_TSPM_UporDown[[1]]
      #
      ###Up 
      hh1<- (DEG_list_RPKM$Up$"AP1.4_AP1.67"); hh2<- (DEG_list_RPKM$Up$"AP3.4_AP3.67"); hh3<- (DEG_list_RPKM$Up$"AG.4_AG.67"); hh4<- (DEG_list_RPKM$Up$"AP1.4_AP3.4"); hh5<- (DEG_list_RPKM$Up$"AP1.4_AG.4"); hh6<- (DEG_list_RPKM$Up$"AP3.4_AG.4"); hh7<- (DEG_list_RPKM$Up$"AP1.67_AP3.67"); hh8<- (DEG_list_RPKM$Up$"AP1.67_AG.67"); hh9<- (DEG_list_RPKM$Up$"AP3.67_AG.67"); 
      DEG_list_RPKM_Up <-Reduce(union, list(hh1,hh2,hh3,hh4,hh5,hh6,hh7,hh8,hh9))
      dfDEG_RPKM_Up <- data.frame(DEG_list_RPKM_Up); rownames(dfDEG_RPKM_Up) <- dfDEG_RPKM_Up[[1]]
      #
      ii1<- (DEG_list_edgeR$Up$"AP1.4.AP1.67"); ii2<- (DEG_list_edgeR$Up$"AP3.4.AP3.67"); ii3<- (DEG_list_edgeR$Up$"AG.4.AG.67"); ii4<- (DEG_list_edgeR$Up$"AP1.4.AP3.4"); ii5<- (DEG_list_edgeR$Up$"AP1.4.AG.4"); ii6<- (DEG_list_edgeR$Up$"AP3.4.AG.4"); ii7<- (DEG_list_edgeR$Up$"AP1.67.AP3.67"); ii8<- (DEG_list_edgeR$Up$"AP1.67.AG.67"); ii9<- (DEG_list_edgeR$Up$"AP3.67.AG.67"); 
      DEG_list_edgeR_Up <-Reduce(union, list(ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9))
      dfDEG_edgeR_Up <- data.frame(DEG_list_edgeR_Up); rownames(dfDEG_edgeR_Up) <- dfDEG_edgeR_Up[[1]]
      #
      jj1<- (DEG_list_DESeq2$Up$"AP1.4.AP1.67"); jj2<- (DEG_list_DESeq2$Up$"AP3.4.AP3.67"); jj3<- (DEG_list_DESeq2$Up$"AG.4.AG.67"); jj4<- (DEG_list_DESeq2$Up$"AP1.4.AP3.4"); jj5<- (DEG_list_DESeq2$Up$"AP1.4.AG.4"); jj6<- (DEG_list_DESeq2$Up$"AP3.4.AG.4"); jj7<- (DEG_list_DESeq2$Up$"AP1.67.AP3.67"); jj8<- (DEG_list_DESeq2$Up$"AP1.67.AG.67"); jj9<- (DEG_list_DESeq2$Up$"AP3.67.AG.67"); 
      DEG_list_DESeq2_Up <-Reduce(union, list(jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9))
      dfDEG_DESeq2_Up <- data.frame(DEG_list_DESeq2_Up); rownames(dfDEG_DESeq2_Up) <- dfDEG_DESeq2_Up[[1]]
      #
      ll1<- (DEG_list_NBPSeq.glmDF$Up$"AP1.4_AP1.67"); ll2<- (DEG_list_NBPSeq.glmDF$Up$"AP3.4_AP3.67"); ll3<- (DEG_list_NBPSeq.glmDF$Up$"AG.4_AG.67"); ll4<- (DEG_list_NBPSeq.glmDF$Up$"AP1.4_AP3.4"); ll5<- (DEG_list_NBPSeq.glmDF$Up$"AP1.4_AG.4"); ll6<- (DEG_list_NBPSeq.glmDF$Up$"AP3.4_AG.4"); ll7<- (DEG_list_NBPSeq.glmDF$Up$"AP1.67_AP3.67"); ll8<- (DEG_list_NBPSeq.glmDF$Up$"AP1.67_AG.67"); ll9<- (DEG_list_NBPSeq.glmDF$Up$"AP3.67_AG.67"); 
      DEG_list_NBPSeq.glmDF_Up <-Reduce(union, list(ll1,ll2,ll3,ll4,ll5,ll6,ll7,ll8,ll9))
      dfDEG_NBPSeq.glmDF_Up <- data.frame(DEG_list_NBPSeq.glmDF_Up); rownames(dfDEG_NBPSeq.glmDF_Up) <- dfDEG_NBPSeq.glmDF_Up[[1]]
      #
      mm1<- (DEG_list_NBPSeq.nbpDF$Up$"AP1.4_AP1.67"); mm2<- (DEG_list_NBPSeq.nbpDF$Up$"AP3.4_AP3.67"); mm3<- (DEG_list_NBPSeq.nbpDF$Up$"AG.4_AG.67"); mm4<- (DEG_list_NBPSeq.nbpDF$Up$"AP1.4_AP3.4"); mm5<- (DEG_list_NBPSeq.nbpDF$Up$"AP1.4_AG.4"); mm6<- (DEG_list_NBPSeq.nbpDF$Up$"AP3.4_AG.4"); mm7<- (DEG_list_NBPSeq.nbpDF$Up$"AP1.67_AP3.67"); mm8<- (DEG_list_NBPSeq.nbpDF$Up$"AP1.67_AG.67"); mm9<- (DEG_list_NBPSeq.nbpDF$Up$"AP3.67_AG.67"); 
      DEG_list_NBPSeq.nbpDF_Up <-Reduce(union, list(mm1,mm2,mm3,mm4,mm5,mm6,mm7,mm8,mm9))
      dfDEG_NBPSeq.nbpDF_Up <- data.frame(DEG_list_NBPSeq.nbpDF_Up); rownames(dfDEG_NBPSeq.nbpDF_Up) <- dfDEG_NBPSeq.nbpDF_Up[[1]]
      #
      nn1<- (DEG_list_TSPM$Up$"AP1.4_AP1.67"); nn2<- (DEG_list_TSPM$Up$"AP3.4_AP3.67"); nn3<- (DEG_list_TSPM$Up$"AG.4_AG.67"); nn4<- (DEG_list_TSPM$Up$"AP1.4_AP3.4"); nn5<- (DEG_list_TSPM$Up$"AP1.4_AG.4"); nn6<- (DEG_list_TSPM$Up$"AP3.4_AG.4"); nn7<- (DEG_list_TSPM$Up$"AP1.67_AP3.67"); nn8<- (DEG_list_TSPM$Up$"AP1.67_AG.67"); nn9<- (DEG_list_TSPM$Up$"AP3.67_AG.67"); 
      DEG_list_TSPM_Up <-Reduce(union, list(nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,nn9))
      dfDEG_TSPM_Up <- data.frame(DEG_list_TSPM_Up); rownames(dfDEG_TSPM_Up) <- dfDEG_TSPM_Up[[1]]
      #
      ###Down
      oo1<- (DEG_list_RPKM$Down$"AP1.4_AP1.67"); oo2<- (DEG_list_RPKM$Down$"AP3.4_AP3.67"); oo3<- (DEG_list_RPKM$Down$"AG.4_AG.67"); oo4<- (DEG_list_RPKM$Down$"AP1.4_AP3.4"); oo5<- (DEG_list_RPKM$Down$"AP1.4_AG.4"); oo6<- (DEG_list_RPKM$Down$"AP3.4_AG.4"); oo7<- (DEG_list_RPKM$Down$"AP1.67_AP3.67"); oo8<- (DEG_list_RPKM$Down$"AP1.67_AG.67"); oo9<- (DEG_list_RPKM$Down$"AP3.67_AG.67"); 
      DEG_list_RPKM_Down <-Reduce(union, list(oo1,oo2,oo3,oo4,oo5,oo6,oo7,oo8,oo9))
      dfDEG_RPKM_Down <- data.frame(DEG_list_RPKM_Down); rownames(dfDEG_RPKM_Down) <- dfDEG_RPKM_Down[[1]]
      #
      pp1<- (DEG_list_edgeR$Down$"AP1.4.AP1.67"); pp2<- (DEG_list_edgeR$Down$"AP3.4.AP3.67"); pp3<- (DEG_list_edgeR$Down$"AG.4.AG.67"); pp4<- (DEG_list_edgeR$Down$"AP1.4.AP3.4"); pp5<- (DEG_list_edgeR$Down$"AP1.4.AG.4"); pp6<- (DEG_list_edgeR$Down$"AP3.4.AG.4"); pp7<- (DEG_list_edgeR$Down$"AP1.67.AP3.67"); pp8<- (DEG_list_edgeR$Down$"AP1.67.AG.67"); pp9<- (DEG_list_edgeR$Down$"AP3.67.AG.67"); 
      DEG_list_edgeR_Down <-Reduce(union, list(pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9))
      dfDEG_edgeR_Down <- data.frame(DEG_list_edgeR_Down); rownames(dfDEG_edgeR_Down) <- dfDEG_edgeR_Down[[1]]
      #
      qq1<- (DEG_list_DESeq2$Down$"AP1.4.AP1.67"); qq2<- (DEG_list_DESeq2$Down$"AP3.4.AP3.67"); qq3<- (DEG_list_DESeq2$Down$"AG.4.AG.67"); qq4<- (DEG_list_DESeq2$Down$"AP1.4.AP3.4"); qq5<- (DEG_list_DESeq2$Down$"AP1.4.AG.4"); qq6<- (DEG_list_DESeq2$Down$"AP3.4.AG.4"); qq7<- (DEG_list_DESeq2$Down$"AP1.67.AP3.67"); qq8<- (DEG_list_DESeq2$Down$"AP1.67.AG.67"); qq9<- (DEG_list_DESeq2$Down$"AP3.67.AG.67"); 
      DEG_list_DESeq2_Down <-Reduce(union, list(qq1,qq2,qq3,qq4,qq5,qq6,qq7,qq8,qq9))
      dfDEG_DESeq2_Down <- data.frame(DEG_list_DESeq2_Down); rownames(dfDEG_DESeq2_Down) <- dfDEG_DESeq2_Down[[1]]
      #
      ss1<- (DEG_list_NBPSeq.glmDF$Down$"AP1.4_AP1.67"); ss2<- (DEG_list_NBPSeq.glmDF$Down$"AP3.4_AP3.67"); ss3<- (DEG_list_NBPSeq.glmDF$Down$"AG.4_AG.67"); ss4<- (DEG_list_NBPSeq.glmDF$Down$"AP1.4_AP3.4"); ss5<- (DEG_list_NBPSeq.glmDF$Down$"AP1.4_AG.4"); ss6<- (DEG_list_NBPSeq.glmDF$Down$"AP3.4_AG.4"); ss7<- (DEG_list_NBPSeq.glmDF$Down$"AP1.67_AP3.67"); ss8<- (DEG_list_NBPSeq.glmDF$Down$"AP1.67_AG.67"); ss9<- (DEG_list_NBPSeq.glmDF$Down$"AP3.67_AG.67"); 
      DEG_list_NBPSeq.glmDF_Down <-Reduce(union, list(ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9))
      dfDEG_NBPSeq.glmDF_Down <- data.frame(DEG_list_NBPSeq.glmDF_Down); rownames(dfDEG_NBPSeq.glmDF_Down) <- dfDEG_NBPSeq.glmDF_Down[[1]]
      #
      tt1<- (DEG_list_NBPSeq.nbpDF$Down$"AP1.4_AP1.67"); tt2<- (DEG_list_NBPSeq.nbpDF$Down$"AP3.4_AP3.67"); tt3<- (DEG_list_NBPSeq.nbpDF$Down$"AG.4_AG.67"); tt4<- (DEG_list_NBPSeq.nbpDF$Down$"AP1.4_AP3.4"); tt5<- (DEG_list_NBPSeq.nbpDF$Down$"AP1.4_AG.4"); tt6<- (DEG_list_NBPSeq.nbpDF$Down$"AP3.4_AG.4"); tt7<- (DEG_list_NBPSeq.nbpDF$Down$"AP1.67_AP3.67"); tt8<- (DEG_list_NBPSeq.nbpDF$Down$"AP1.67_AG.67"); tt9<- (DEG_list_NBPSeq.nbpDF$Down$"AP3.67_AG.67"); 
      DEG_list_NBPSeq.nbpDF_Down <-Reduce(union, list(tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9))
      dfDEG_NBPSeq.nbpDF_Down <- data.frame(DEG_list_NBPSeq.nbpDF_Down); rownames(dfDEG_NBPSeq.nbpDF_Down) <- dfDEG_NBPSeq.nbpDF_Down[[1]]
      #
      uu1<- (DEG_list_TSPM$Down$"AP1.4_AP1.67"); uu2<- (DEG_list_TSPM$Down$"AP3.4_AP3.67"); uu3<- (DEG_list_TSPM$Down$"AG.4_AG.67"); uu4<- (DEG_list_TSPM$Down$"AP1.4_AP3.4"); uu5<- (DEG_list_TSPM$Down$"AP1.4_AG.4"); uu6<- (DEG_list_TSPM$Down$"AP3.4_AG.4"); uu7<- (DEG_list_TSPM$Down$"AP1.67_AP3.67"); uu8<- (DEG_list_TSPM$Down$"AP1.67_AG.67"); uu9<- (DEG_list_TSPM$Down$"AP3.67_AG.67"); 
      DEG_list_TSPM_Down <-Reduce(union, list(uu1,uu2,uu3,uu4,uu5,uu6,uu7,uu8,uu9))
      dfDEG_TSPM_Down <- data.frame(DEG_list_TSPM_Down); rownames(dfDEG_TSPM_Down) <- dfDEG_TSPM_Down[[1]]
      #

################################################################ 
###Venn diagram
################################################################ 
###Union all comparison 

###Up
setlist2 <- list(edgeR=rownames(DEG_list_edgeR_Up), DESeq2=rownames(DEG_list_DESeq2_Up), TSPM=rownames(DEG_list_TSPM_Up), NBPSeq.glm=rownames(DEG_list_NBPSeq.glmDF_Up), NBPSeq.nbp=rownames(DEG_list_NBPSeq.nbpDF_Up) )
OLlist2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
counts2 <- sapply(OLlist2$Venn_List, length)
pdf("./results/Venn_diagram2.pdf")
vennPlot(counts=counts2, mymain="DEG Comparison Up")
dev.off()

###Down
setlist4 <- list(edgeR=rownames(DEG_list_edgeR_Down), DESeq2=rownames(DEG_list_DESeq2_Down), TSPM=rownames(DEG_list_TSPM_Down), NBPSeq.glm=rownames(DEG_list_NBPSeq.glmDF_Down), NBPSeq.nbp=rownames(DEG_list_NBPSeq.nbpDF_Down) )
OLlist4 <- overLapper(setlist=setlist4, sep="_", type="vennsets")
counts4 <- sapply(OLlist4$Venn_List, length)
pdf("./results/Venn_diagram4.pdf")
vennPlot(counts=counts4, mymain="DEG Comparison Down")
dev.off()
 
###Specific comparison
hh4<- (DEG_list_RPKM$Up$"AP1.4_AP3.4"); hh5<- (DEG_list_RPKM$Up$"AP1.4_AG.4"); hh6<- (DEG_list_RPKM$Up$"AP3.4_AG.4")
RPKMhh4 <- data.frame(hh4); rownames(RPKMhh4) <- RPKMhh4[[1]]
RPKMhh5 <- data.frame(hh5); rownames(RPKMhh5) <- RPKMhh5[[1]]
RPKMhh6 <- data.frame(hh6); rownames(RPKMhh6) <- RPKMhh6[[1]]
 
setlist6 <- list("AP1.4-AP3.4"=rownames(RPKMhh4), "AP1.4-AG.4"=rownames(RPKMhh5), "AP3.4-AG.4"=rownames(RPKMhh6))
OLlist6 <- overLapper(setlist=setlist6, sep="_", type="vennsets")
counts6 <- sapply(OLlist6$Venn_List, length)
pdf("./results/Venn_diagram_RPKM.pdf")
vennPlot(counts=counts6, mymain="DEG Comparison RPKM")
dev.off()
 
################################################################ 
###Scatterplot
################################################################ 
###Data: scatterDEG_logFC
      RPKM.S <- RPKM_FC[unlist(DEG_list_RPKM_UporDown),]; edgeR.S <- edgeDF[unlist(DEG_list_edgeR_UporDown),]; deseq2.S <- deseq2DF[unlist(DEG_list_DESeq2_UporDown),]
      NBPSeq.glm.S <- NBPSeq.glmDF[unlist(DEG_list_NBPSeq.glmDF_UporDown),]; NBPSeq.nbp.S <- NBPSeq.nbpDF[unlist(DEG_list_NBPSeq.nbpDF_UporDown),]; TSPM.S <- TSPMDF [unlist(DEG_list_TSPM_UporDown),]
      
      RPKMScatter <- data.frame(rownames(RPKM.S), RPKM.S$"AP3.67_AG.67_logFC", row.names=1) 
      edgeRScatter <- data.frame(rownames(edgeR.S), edgeR.S$"AP3.67.AG.67_logFC", row.names=1)
      deseq2Scatter <- data.frame(rownames(deseq2.S), deseq2.S$"AP3.67.AG.67_logFC" , row.names=1)
      NBPSeq.glmScatter <- data.frame(rownames(NBPSeq.glm.S), NBPSeq.glm.S$"AP3.67_AG.67_logFC" , row.names=1)
      NBPSeq.nbpScatter <- data.frame(rownames(NBPSeq.nbp.S), NBPSeq.nbp.S$"AP3.67_AG.67_logFC" , row.names=1)
      TSPMScatter <- data.frame(rownames(TSPM.S), TSPM.S$"AP3.67_AG.67_logFC", row.names=1)  
      
      scatterDEG1 <- merge(edgeRScatter, deseq2Scatter, by='row.names', all=TRUE); rownames(scatterDEG1) <- scatterDEG1[[1]] ; scatterDEG1 <- scatterDEG1 [,-1]
      scatterDEG2 <- merge(scatterDEG1, RPKMScatter, by='row.names', all=TRUE); rownames(scatterDEG2) <- scatterDEG2[[1]] ; scatterDEG2 <- scatterDEG2 [,-1]
      scatterDEG3 <- merge(scatterDEG2, NBPSeq.glmScatter, by='row.names', all=TRUE); rownames(scatterDEG3) <- scatterDEG3[[1]] ; scatterDEG3 <- scatterDEG3 [,-1]
      scatterDEG4 <- merge(scatterDEG3, NBPSeq.nbpScatter, by='row.names', all=TRUE); rownames(scatterDEG4) <- scatterDEG4[[1]] ; scatterDEG4 <- scatterDEG4 [,-1]
      scatterDEG5 <- merge(scatterDEG4, TSPMScatter, by='row.names', all=TRUE); rownames(scatterDEG5) <- scatterDEG5[[1]] ; scatterDEG5 <- scatterDEG5 [,-1]
      
      scatterDEG_logFC <- scatterDEG5
      colnames(scatterDEG_logFC) <- c("edgeR", "DESeq2", "RPKM", "NBPSeq.glm", "NBPSeq.nbp", "TSPM")
      scatterDEG_logFC[is.na(scatterDEG_logFC)] <- 0 

###Scatterplot
scatterDEG_logFC[1:4,]
source("./script/panel.cor.R")
pdf("./results/Scatterplot_AP3.67_AG.67_logFC.pdf")
pairs(scatterDEG_logFC, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main="Scatterplot Matrix log_FC")
dev.off()

###Data: scatterDEG_FDR
      bayseq.S <- bayseqDF [unlist(DEG_list_bayseqDF),]; edgeR.S <- edgeDF[unlist(DEG_list_edgeR_UporDown),]; deseq2.S <- deseq2DF[unlist(DEG_list_DESeq2_UporDown),]
      NBPSeq.glm.S <- NBPSeq.glmDF[unlist(DEG_list_NBPSeq.glmDF_UporDown),]; NBPSeq.nbp.S <- NBPSeq.nbpDF[unlist(DEG_list_NBPSeq.nbpDF_UporDown),]; TSPM.S <- TSPMDF [unlist(DEG_list_TSPM_UporDown),]
      
      bayseqScatter <- data.frame(rownames(bayseq.S), bayseq.S$"AP3.67_AG.67_FDR", row.names=1) 
      edgeRScatter <- data.frame(rownames(edgeR.S), edgeR.S$"AP3.67.AG.67_FDR", row.names=1)
      deseq2Scatter <- data.frame(rownames(deseq2.S), deseq2.S$"AP3.67.AG.67_FDR" , row.names=1)
      NBPSeq.glmScatter <- data.frame(rownames(NBPSeq.glm.S), NBPSeq.glm.S$"AP3.67_AG.67_FDR" , row.names=1)
      NBPSeq.nbpScatter <- data.frame(rownames(NBPSeq.nbp.S), NBPSeq.nbp.S$"AP3.67_AG.67_FDR" , row.names=1)
      TSPMScatter <- data.frame(rownames(TSPM.S), TSPM.S$"AP3.67_AG.67_FDR", row.names=1)  
      
      scatterDEG1 <- merge(edgeRScatter, deseq2Scatter, by='row.names', all=TRUE); rownames(scatterDEG1) <- scatterDEG1[[1]] ; scatterDEG1 <- scatterDEG1 [,-1]
      scatterDEG2 <- merge(scatterDEG1, bayseqScatter, by='row.names', all=TRUE); rownames(scatterDEG2) <- scatterDEG2[[1]] ; scatterDEG2 <- scatterDEG2 [,-1]
      scatterDEG3 <- merge(scatterDEG2, NBPSeq.glmScatter, by='row.names', all=TRUE); rownames(scatterDEG3) <- scatterDEG3[[1]] ; scatterDEG3 <- scatterDEG3 [,-1]
      scatterDEG4 <- merge(scatterDEG3, NBPSeq.nbpScatter, by='row.names', all=TRUE); rownames(scatterDEG4) <- scatterDEG4[[1]] ; scatterDEG4 <- scatterDEG4 [,-1]
      scatterDEG5 <- merge(scatterDEG4, TSPMScatter, by='row.names', all=TRUE); rownames(scatterDEG5) <- scatterDEG5[[1]] ; scatterDEG5 <- scatterDEG5 [,-1]
      
      scatterDEG_FDR <- scatterDEG5
      colnames(scatterDEG_FDR) <- c("edgeR", "DESeq2", "baySeq", "NBPSeq.glm", "NBPSeq.nbp", "TSPM")
      scatterDEG_FDR[is.na(scatterDEG_FDR)] <- 1

###Scatterplot
scatterDEG_FDR[1:4,]
source("./script/panel.cor.R")
pdf("./results/Scatterplot_AP3.67_AG.67_FDR.pdf")
pairs(scatterDEG_FDR, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main="Scatterplot Matrix FDR")
dev.off()

################################################################ 
###ROC Curve
################################################################
###Data
      bayseqComp <- bayseqDF[unlist(DEG_list_bayseqDF$UporDown$"AP3.67_AG.67"),]; bayseqComp <- bayseqComp[71]; colnames(bayseqComp) <- list("Bayseq")
      edgeComp <- edgeDF[unlist(DEG_list_edgeR$UporDown$"AP3.67.AG.67"),]; edgeComp <- edgeComp[45]; colnames(edgeComp) <- list("edgeR")
      deseq2Comp <- deseq2DF[unlist(DEG_list_DESeq2$UporDown$"AP3.67.AG.67"),]; deseq2Comp <- deseq2Comp[54]; colnames(deseq2Comp) <- list("DESeq2")
      NBPSeq.glmComp <- NBPSeq.glmDF[unlist(DEG_list_NBPSeq.glmDF$UporDown$"AP3.67_AG.67"),]; NBPSeq.glmComp <- NBPSeq.glmComp[63]; colnames(NBPSeq.glmComp) <- list("NBPSeq.glm")
      NBPSeq.nbpComp <- NBPSeq.nbpDF[unlist(DEG_list_NBPSeq.nbpDF$UporDown$"AP3.67_AG.67"),]; NBPSeq.nbpComp <- NBPSeq.nbpComp[27]; colnames(NBPSeq.nbpComp) <- list("NBPSeq.nbp")
      TSPMComp <- TSPMDF[unlist(DEG_list_TSPM$UporDown$"AP3.67_AG.67"),]; TSPMComp <- TSPMComp[27]; colnames(TSPMComp) <- list("TSPM")
      
      tmp <- merge (edgeComp, deseq2Comp, by='row.names', all=TRUE); rownames(tmp) <- tmp[[1]] ; tmp <- tmp [,-1]
      tmp2 <- merge (tmp, bayseqComp, by='row.names', all=TRUE); rownames(tmp2) <- tmp2[[1]] ; tmp2 <- tmp2 [,-1]
      tmp3 <- merge (tmp2, NBPSeq.glmComp, by='row.names', all=TRUE); rownames(tmp3) <- tmp3[[1]] ; tmp3 <- tmp3 [,-1]
      tmp4 <- merge (tmp3, NBPSeq.nbpComp, by='row.names', all=TRUE); rownames(tmp4) <- tmp4[[1]] ; tmp4 <- tmp4 [,-1]
      tmp5 <- merge (tmp4, TSPMComp, by='row.names', all=TRUE); rownames(tmp5) <- tmp5[[1]] ; tmp5 <- tmp5 [,-1]
      
      data.ROC <- tmp5
      tmp[is.na(tmp)] <- 1; 
      data.ROC <- merge(data.ROC, tmp[1], by='row.names', all=TRUE ); rownames(data.ROC) <- data.ROC[[1]] ; data.ROC <- data.ROC [,-1]
      colnames(data.ROC) <- c("edgeR", "DESeq2", "Bayseq", "NBPSeq.glm", "NBPSeq.nbp", "TSPM", "Common")
      write.table(data.ROC, "results/data_ROC.xls", col.names=NA, quote=FALSE, sep="\t")
      
      data.class <- data.ROC
      data.class[!is.na(data.class)] <- 0
      data.class[is.na(data.class)] <- 1
      write.table(data.class, "results/data.class.xls", col.names=NA, quote=FALSE, sep="\t")
      
      data.ROC1 <- data.ROC
      data.ROC1[,1:7][is.na(data.ROC1[,1:7])] <- 1 

### ROCR
library ("ROCR")
 
pdf("./results/ROC.pdf")
pred1 <- prediction(data.ROC1$Common, data.class$edgeR)
perf1 <- performance(pred1, "tpr", "fpr")
plot(perf1, avg= "threshold", col="black", lty=4,lwd= 2, main= "ROC curves compare with edgeR")
par(new = TRUE)
pred2 <- prediction(data.ROC1$Common, data.class$Bayseq)
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2, avg= "threshold", lty=4, lwd= 2, col="dodgerblue4")
par(new = TRUE)
pred3 <- prediction(data.ROC1$Common, data.class$NBPSeq.glm)
perf3 <- performance(pred3, "tpr", "fpr")
plot(perf3, avg= "threshold",  lty=4,lwd= 2, col="darkgreen")
par(new = TRUE)
pred4 <- prediction(data.ROC1$Common, data.class$NBPSeq.nbp)
perf4 <- performance(pred4, "tpr", "fpr")
plot(perf4, avg= "threshold",  lty=4,lwd= 2, col="darkviolet")
par(new = TRUE)
pred5 <- prediction(data.ROC1$Common, data.class$TSPM)
perf5 <- performance(pred5, "tpr", "fpr")
plot(perf5, avg= "threshold", lty=4, lwd= 2, col="firebrick")
par(new = TRUE)
pred6 <- prediction(data.ROC1$Common, data.class$DESeq2)
perf6 <- performance(pred6, "tpr", "fpr")
plot(perf6, avg= "threshold",  lty=4,lwd= 2, col="deeppink1")
legend("bottomright", legend=c("edgeR / 1", "BaySeq / 0.749", "NBPSeq.glm / 0.533", "NBPSeq.nbp / 0.872 ", "TSPM / 0.623", "DESeq2 / 0.974"),
       col=c("black", "dodgerblue4", "darkgreen", "darkviolet", "firebrick", "deeppink1"), lty=4, lwd=3, bty="n")
 
dev.off()

##calculate AUC
auc.tmp1 <- performance(pred1,"auc"); auc1 <- as.numeric(auc.tmp1@y.values)
auc.tmp2 <- performance(pred2,"auc"); auc2 <- as.numeric(auc.tmp2@y.values) 
auc.tmp3 <- performance(pred3,"auc"); auc3 <- as.numeric(auc.tmp3@y.values)
auc.tmp4 <- performance(pred4,"auc"); auc4 <- as.numeric(auc.tmp4@y.values) 
auc.tmp5 <- performance(pred5,"auc"); auc5 <- as.numeric(auc.tmp5@y.values)
auc.tmp6 <- performance(pred6,"auc"); auc6 <- as.numeric(auc.tmp6@y.values) 
 
 
################################################################ 
### Correlation dendrogram of methods
################################################################
library(ape)

dFDR <- cor(data.ROC1, method="spearman")
hc1 <- hclust(as.dist(1-dFDR))
pdf("./results/methodsCorrelation_FDR.pdf")
plot.phylo(as.phylo(hc1), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE,  main="Correlation FDR")
dev.off()

dlogFC <- cor(scatterDEG_logFC, method="spearman")
hc2 <- hclust(as.dist(1-dlogFC))
pdf("./results/methodsCorrelation_logFC.pdf")
plot.phylo(as.phylo(hc2), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, main="Correlation logFC")
dev.off()


