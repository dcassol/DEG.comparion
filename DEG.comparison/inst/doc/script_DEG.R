### R code from vignette source 'script_DEG.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: script_DEG.Rnw:42-44
###################################################
options(width=95)
unlink("test.db")


###################################################
### code chunk number 3: script_DEG.Rnw:74-75
###################################################
library(DEG.comparison)


###################################################
### code chunk number 4: script_DEG.Rnw:80-85
###################################################
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment.char = "#")
targets
cmp <- readComp(file=targetspath, format="matrix", delim="-")


###################################################
### code chunk number 5: script_DEG.Rnw:90-92
###################################################
parampath <- system.file("extdata", "tophat.param", package="DEG.comparison")
read.delim(parampath, comment.char = "#")


###################################################
### code chunk number 6: script_DEG.Rnw:95-97
###################################################
args <- systemArgs(sysma=parampath, mytargets=targetspath)
args


###################################################
### code chunk number 7: script_DEG.Rnw:100-105
###################################################
names(args)
modules(args)
cores(args)
outpaths(args)[1]
sysargs(args)[1]


###################################################
### code chunk number 8: script_DEG.Rnw:108-109
###################################################
systemArgs(sysma=parampath, mytargets=targetspath, type="json")


###################################################
### code chunk number 9: script_DEG.Rnw:115-118 (eval = FALSE)
###################################################
## trim.param <- system.file("extdata", "trim.param", package="DEG.comparison")
## trim.param <- read.delim(trim.param, comment.char = "#")
## args <- systemArgs(sysma="trim.param", mytargets=targetspath)


###################################################
### code chunk number 10: script_DEG.Rnw:124-127 (eval = FALSE)
###################################################
## preprocessReads(args=args, Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", 
##                 batchsize=100000, overwrite=TRUE, compress=TRUE)
## writeTargetsout(x=args, file="targets_trim.txt")


###################################################
### code chunk number 11: script_DEG.Rnw:132-139 (eval = FALSE)
###################################################
## filterFct <- function(fq) {
##   filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
## 	filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes low complexity reads 
## 	filter <- compose(filter1, filter2)
## 	fq[filter(fq)]
## }
## preprocessReads(args=args, Fct="filterFct(fq)", batchsize=100000)


###################################################
### code chunk number 12: script_DEG.Rnw:144-148 (eval = FALSE)
###################################################
## fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
## pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
## seeFastqPlot(fqlist)
## dev.off()


###################################################
### code chunk number 13: script_DEG.Rnw:160-162 (eval = FALSE)
###################################################
## args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
## sysargs(args)[1] # Command-line parameters for first FASTQ file


###################################################
### code chunk number 14: script_DEG.Rnw:165-174 (eval = FALSE)
###################################################
## moduleload(modules(args))
## system("bowtie2-build ./data/TAIR10_chr_all.fas ./data/TAIR10_chr_all.fas")
## bampaths <- runCommandline(args=args)
## file.copy(paste0("./.BatchJobs.R"), ".")
## file.copy(paste0("./torque.tmpl"), ".")
## resources <- list(walltime="20:00:00", nodes=paste0("1:ppn=", cores(args)), memory="16gb")
## reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", 
##                   resourceList=resources)
## showStatus(reg)


###################################################
### code chunk number 15: script_DEG.Rnw:177-179 (eval = FALSE)
###################################################
## file.exists(outpaths(args))
## sapply(1:length(args), function(x) loadResult(reg, x)) # Works after job completion


###################################################
### code chunk number 16: script_DEG.Rnw:184-186 (eval = FALSE)
###################################################
## read_statsDF <- alignStats(args=args) 
## write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 17: script_DEG.Rnw:188-189
###################################################
read.table(system.file("extdata", "alignStats.xls", package="DEG.comparison"), header=TRUE)[1:4,]


###################################################
### code chunk number 18: script_DEG.Rnw:194-197 (eval = FALSE)
###################################################
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "cassol_GEN242/"), 
##             urlbase="http://biocluster.ucr.edu/~dcassol/", urlfile="IGVurl.txt")
## 


###################################################
### code chunk number 19: script_DEG.Rnw:203-223 (eval = FALSE)
###################################################
## library("GenomicFeatures"); library(BiocParallel)
## txdb <- makeTranscriptDbFromGFF(file="data/TAIR10_GFF3_genes.gff", format="gff3", dataSource="TAIR", species="Arabidopsis thaliana")
## saveDb(txdb, file="./data/tair10.sqlite")
## 
## txdb <- loadDb("./data/tair10.sqlite")
## eByg <- exonsBy(txdb, by=c("gene"))
## bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
## multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
## # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE' 
## counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union",
##                                                          ignore.strand=TRUE, 
##                                                          inter.feature=TRUE, 
##                                                          singleEnd=TRUE))
## countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## rownames(countDFeByg) <- names(rowData(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
## countDFeByg[1:4,1:12]
## write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
## write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
## rpkmDFeByg[1:4,1:7]


###################################################
### code chunk number 20: script_DEG.Rnw:228-237 (eval = FALSE)
###################################################
## library(ape)
## rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="DEG.comparison")
## rpkmDFeByg <- read.delim(rpkmDFeBygpath, row.names=1)
## rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
## d <- cor(rpkmDFeByg, method="spearman")
## hc <- hclust(as.dist(1-d))
## pdf("results/sample_tree.pdf")
## plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
## dev.off()


###################################################
### code chunk number 21: script_DEG.Rnw:251-256 (eval = FALSE)
###################################################
## ##Data input
## countDFeBygpath <- system.file("extdata", "countDFeByg.xls", package="DEG.comparison")
## countDFeByg <- read.delim(countDFeBygpath, row.names=1)
## rpkmDFeBygpath <- system.file("extdata", "rpkmDFeByg.xls", package="DEG.comparison")
## rpkmDFeByg <- read.delim(rpkmDFeBygpath, row.names=1)


###################################################
### code chunk number 22: script_DEG.Rnw:261-271 (eval = FALSE)
###################################################
## ##Comparisons
## Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"), 
##               AP3.4_AP3.67=c("AP3.4A","AP3.4B", "AP3.67A" ,"AP3.67B"), 
##               AG.4_AG.67=c("AG.4A", "AG.4B", "AG.67A","AG.67B"), 
##               AP1.4_AP3.4=c("AP1.4A","AP1.4B", "AP3.4A","AP3.4B"), 
##               AP1.4_AG.4=c("AP1.4A", "AP1.4B", "AG.4A", "AG.4B"), 
##               AP3.4_AG.4=c("AP3.4A", "AP3.4B", "AG.4A", "AG.4B"), 
##               AP1.67_AP3.67=c("AP1.67A","AP1.67B", "AP3.67A", "AP3.67B"), 
##               AP1.67_AG.67=c("AP1.67A","AP1.67B", "AG.67A", "AG.67B"), 
##               AP3.67_AG.67=c("AP3.67A", "AP3.67B", "AG.67A", "AG.67B"))


###################################################
### code chunk number 23: script_DEG.Rnw:278-287 (eval = FALSE)
###################################################
## ##Settings
## Comp1 <- list(Factor=(Reduce(union, targets$Factor)), Sample=c(colnames(rpkmDFeByg)), 
##               group=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
## 
## Comp2 <- list(AP1.4_AP1.67=c("AP1.4", "AP1.67"), AP3.4_AP3.67=c("AP3.4", "AP3.67"), 
##               AG.4_AG.67=c("AG.4", "AG.67"), AP1.4_AP3.4=c("AP1.4", "AP3.4"), 
##               AP1.4_AG.4=c("AP1.4", "AG.4"), AP3.4_AG.4=c("AP3.4", "AG.4"), 
##               AP1.67_AP3.67=c("AP1.67", "AP3.67"), AP1.67_AG.67=c("AP1.67", "AG.67"), 
##               AP3.67_AG.67=c("AP3.67", "AG.67"))


###################################################
### code chunk number 24: script_DEG.Rnw:291-297 (eval = FALSE)
###################################################
## RPKM_FC <- run_RPKM (rpkmDFeByg, Comp1, Comp2)
## write.table(RPKM_FC, "./results/RPKM_FC.xls", quote=FALSE, sep="\t", col.names = NA)
## pdf("./results/DEG_list_RPKM.pdf") 
## DEG_list_RPKM <- filterDEG_logFC(degDF=RPKM_FC, filter=c(Fold=2), method="RPKM")
## dev.off()
## DEG_list_RPKM$Summary[1:4,]


###################################################
### code chunk number 25: script_DEG.Rnw:312-318 (eval = FALSE)
###################################################
## edgeDF <- run_edgeR(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")
## write.table(edgeDF, "results/edgeDF.xls", col.names=NA, quote=FALSE, sep="\t")
## pdf("./results/DEG_list_edgeR.pdf") 
## DEG_list_edgeR <- filterDEGnew(degDF=edgeDF, filter=c(Fold=2, FDR=1), method="edgeR")
## dev.off()
## DEG_list_edgeR$Summary[1:4,]


###################################################
### code chunk number 26: script_DEG.Rnw:330-336 (eval = FALSE)
###################################################
## deseq2DF <- run_DESeq2(countDF=countDFeByg, targets=targets, cmp=cmp[[1]], independent=FALSE)
## write.table(deseq2DF, "results/deseq2DF.xls", col.names=NA, quote=FALSE, sep="\t")
## pdf("./results/DEG_list_DESeq2.pdf")
## DEG_list_DESeq2 <- filterDEGnew(degDF=deseq2DF, filter=c(Fold=2, FDR=1), method="DESeq2")
## dev.off()
## DEG_list_DESeq2$Summary[1:4,]


###################################################
### code chunk number 27: script_DEG.Rnw:348-355 (eval = FALSE)
###################################################
## dim(countDFeByg)
## bayseqDF <- run_BaySeq(countDFeByg, Comp3, number=27416)
## write.table(bayseqDF, "results/bayseqDF.xls", col.names=NA, quote=FALSE, sep="\t")
## pdf("./results/DEG_list_bayseqDF.pdf") 
## DEG_list_bayseqDF <- filterDEG_FDR(degDF=bayseqDF, filter=c(FDR=1), method="BaySeq")
## dev.off()
## DEG_list_bayseqDF$Summary[1:4,]


###################################################
### code chunk number 28: script_DEG.Rnw:369-375 (eval = FALSE)
###################################################
## NBPSeq.glmDF <- run_NBPSeq_glm (countDFeByg, Comp3)
## write.table(NBPSeq.glmDF, "results/NBPSeq_glmDF.xls", col.names=NA, quote=FALSE, sep="\t")
## pdf("./results/DEG_list_NBPSeq_glmDF.pdf")
## DEG_list_NBPSeq.glmDF <- filterDEGnew(degDF=NBPSeq.glmDF, filter=c(Fold=2, FDR=1), method="NBPSeq.glm")
## dev.off()
## DEG_list_NBPSeq.glmDF$Summary[1:4,]


###################################################
### code chunk number 29: script_DEG.Rnw:386-392 (eval = FALSE)
###################################################
## NBPSeq.nbpDF <- run_NBPSeq_nbp (countDFeByg, Comp3)
## write.table(NBPSeq.nbpDF, "results/NBPSeq.nbpDF.xls", col.names=NA, quote=FALSE, sep="\t")
## pdf("./results/DEG_list_NBPSeq.nbpDF.pdf")
## DEG_list_NBPSeq.nbpDF <- filterDEGnew(degDF=NBPSeq.nbpDF, filter=c(Fold=2, FDR=1), method="NBPSeq.nbp")
## dev.off()
## DEG_list_NBPSeq.nbpDF$Summary[1:4,]


###################################################
### code chunk number 30: script_DEG.Rnw:404-410 (eval = FALSE)
###################################################
## TSPMDF <- run_TSPM(countDFeByg, Comp3)
## write.table(TSPMDF, "results/TSPMDF.xls", col.names=NA, quote=FALSE, sep="\t")
## pdf("./results/DEG_list_TSPM.pdf")
## DEG_list_TSPM <- filterDEGnew(degDF=TSPMDF, filter=c(Fold=2, FDR=1), method="TSPM")
## dev.off()
## DEG_list_TSPM$Summary[1:4,]


###################################################
### code chunk number 31: script_DEG.Rnw:422-428 (eval = FALSE)
###################################################
## totalgenes <- system.file("extdata", "totalgenes.csv", package="DEG.comparison")
## totalgenes <- read.delim (totalgenes, sep=",")
## pdf("./results/TotalGenes.pdf")
## plot <- ggplot(totalgenes, aes(Package, Number.of.significantly.differentially.expressed)) + geom_bar(aes(fill = DEG), stat="identity") + facet_wrap(~DEG, ncol=1)
## print(plot) 
## dev.off()


###################################################
### code chunk number 32: script_DEG.Rnw:439-481 (eval = FALSE)
###################################################
## ###UporDown
## List_RPKM_UporDown <- Reduce(union, (DEG_list_RPKM$UporDown))     
## RPKM_UporDown <- data.frame(List_RPKM_UporDown); rownames(RPKM_UporDown) <- RPKM_UporDown[[1]]
## List_edgeR_UporDown <- Reduce(union, (DEG_list_edgeR$UporDown))
## edgeR_UporDown <- data.frame(List_edgeR_UporDown); rownames(edgeR_UporDown) <- edgeR_UporDown[[1]]
## List_DESeq2_UporDown <- Reduce(union, (DEG_list_DESeq2$UporDown))
## DESeq2_UporDown <- data.frame(List_DESeq2_UporDown); rownames(DESeq2_UporDown) <- DESeq2_UporDown[[1]]
## List_bayseqDF_UporDown <- Reduce(union, (DEG_list_bayseqDF$UporDown))
## bayseqDF_UporDown <- data.frame(List_bayseqDF_UporDown); rownames(bayseqDF_UporDown) <- bayseqDF_UporDown[[1]]
## List_NBPSeq.glmDF_UporDown <- Reduce(union, (DEG_list_NBPSeq.glmDF$UporDown))
## NBPSeq.glmDF_UporDown <- data.frame(List_NBPSeq.glmDF_UporDown); rownames(NBPSeq.glmDF_UporDown) <- NBPSeq.glmDF_UporDown[[1]]
## List_NBPSeq.nbpDF_UporDown <- Reduce(union, (DEG_list_NBPSeq.nbpDF$UporDown))
## NBPSeq.nbpDF_UporDown <- data.frame(List_NBPSeq.nbpDF_UporDown); rownames(NBPSeq.nbpDF_UporDown) <- NBPSeq.nbpDF_UporDown[[1]]
## List_TSPM_UporDown <- Reduce(union, (DEG_list_TSPM$UporDown))
## TSPM_UporDown <- data.frame(List_TSPM_UporDown); rownames(TSPM_UporDown) <- TSPM_UporDown[[1]]
## ###Up 
## List_RPKM_Up <- Reduce(union, (DEG_list_RPKM$Up))     
## RPKM_Up <- data.frame(List_RPKM_Up); rownames(RPKM_Up) <- RPKM_Up[[1]]
## List_edgeR_Up <- Reduce(union, (DEG_list_edgeR$Up))
## edgeR_Up <- data.frame(List_edgeR_Up); rownames(edgeR_Up) <- edgeR_Up[[1]]
## List_DESeq2_Up <- Reduce(union, (DEG_list_DESeq2$Up))
## DESeq2_Up <- data.frame(List_DESeq2_Up); rownames(DESeq2_Up) <- DESeq2_Up[[1]]
## List_NBPSeq.glmDF_Up <- Reduce(union, (DEG_list_NBPSeq.glmDF$Up))
## NBPSeq.glmDF_Up <- data.frame(List_NBPSeq.glmDF_Up); rownames(NBPSeq.glmDF_Up) <- NBPSeq.glmDF_Up[[1]]
## List_NBPSeq.nbpDF_Up <- Reduce(union, (DEG_list_NBPSeq.nbpDF$Up))
## NBPSeq.nbpDF_Up <- data.frame(List_NBPSeq.nbpDF_Up); rownames(NBPSeq.nbpDF_Up) <- NBPSeq.nbpDF_Up[[1]]
## List_TSPM_Up <- Reduce(union, (DEG_list_TSPM$Up))
## TSPM_Up <- data.frame(List_TSPM_Up); rownames(TSPM_Up) <- TSPM_Up[[1]]
## ###Down
## List_RPKM_Down <- Reduce(union, (DEG_list_RPKM$Down))     
## RPKM_Down <- data.frame(List_RPKM_Down); rownames(RPKM_Down) <- RPKM_Down[[1]]
## List_edgeR_Down <- Reduce(union, (DEG_list_edgeR$Down))
## edgeR_Down <- data.frame(List_edgeR_Down); rownames(edgeR_Down) <- edgeR_Down[[1]]
## List_DESeq2_Down <- Reduce(union, (DEG_list_DESeq2$Down))
## DESeq2_Down <- data.frame(List_DESeq2_Down); rownames(DESeq2_Down) <- DESeq2_Down[[1]]
## List_NBPSeq.glmDF_Down <- Reduce(union, (DEG_list_NBPSeq.glmDF$Down))
## NBPSeq.glmDF_Down <- data.frame(List_NBPSeq.glmDF_Down); rownames(NBPSeq.glmDF_Down) <- NBPSeq.glmDF_Down[[1]]
## List_NBPSeq.nbpDF_Down <- Reduce(union, (DEG_list_NBPSeq.nbpDF$Down))
## NBPSeq.nbpDF_Down <- data.frame(List_NBPSeq.nbpDF_Down); rownames(NBPSeq.nbpDF_Down) <- NBPSeq.nbpDF_Down[[1]]
## List_TSPM_Down <- Reduce(union, (DEG_list_TSPM$Down))
## TSPM_Down <- data.frame(List_TSPM_Down); rownames(TSPM_Down) <- TSPM_Down[[1]]   
## 


###################################################
### code chunk number 33: script_DEG.Rnw:485-494 (eval = FALSE)
###################################################
## ##Up
## setlist2 <- list(edgeR=rownames(edgeR_Up), DESeq2=rownames(DESeq2_Up), 
##                  TSPM=rownames(TSPM_Up), NBPSeq.glm=rownames(NBPSeq.glmDF_Up), 
##                  NBPSeq.nbp=rownames(NBPSeq.nbpDF_Up))
## OLlist2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
## counts2 <- sapply(OLlist2$Venn_List, length)
## pdf("./results/Venn_diagram2.pdf")
## vennPlot(counts=counts2, mymain="DEG Comparison Up")
## dev.off()


###################################################
### code chunk number 34: script_DEG.Rnw:503-512 (eval = FALSE)
###################################################
## ##Down
## setlist4 <- list(edgeR=rownames(edgeR_Down), DESeq2=rownames(DESeq2_Down), 
##                  TSPM=rownames(TSPM_Down), NBPSeq.glm=rownames(NBPSeq.glmDF_Down), 
##                  NBPSeq.nbp=rownames(NBPSeq.nbpDF_Down))
## OLlist4 <- overLapper(setlist=setlist4, sep="_", type="vennsets")
## counts4 <- sapply(OLlist4$Venn_List, length)
## pdf("./results/Venn_diagram4.pdf")
## vennPlot(counts=counts4, mymain="DEG Comparison Down")
## dev.off()


###################################################
### code chunk number 35: script_DEG.Rnw:521-533 (eval = FALSE)
###################################################
## ###Specific comparison
## hh4<- (DEG_list_RPKM$Up$"AP1.4_AP3.4"); hh5<- (DEG_list_RPKM$Up$"AP1.4_AG.4"); hh6<- (DEG_list_RPKM$Up$"AP3.4_AG.4")
## RPKMhh4 <- data.frame(hh4); rownames(RPKMhh4) <- RPKMhh4[[1]]
## RPKMhh5 <- data.frame(hh5); rownames(RPKMhh5) <- RPKMhh5[[1]]
## RPKMhh6 <- data.frame(hh6); rownames(RPKMhh6) <- RPKMhh6[[1]]
## 
## setlist6 <- list("AP1.4-AP3.4"=rownames(RPKMhh4), "AP1.4-AG.4"=rownames(RPKMhh5), "AP3.4-AG.4"=rownames(RPKMhh6))
## OLlist6 <- overLapper(setlist=setlist6, sep="_", type="vennsets")
## counts6 <- sapply(OLlist6$Venn_List, length)
## pdf("./results/Venn_diagram_RPKM.pdf")
## vennPlot(counts=counts6, mymain="DEG Comparison RPKM")
## dev.off()


###################################################
### code chunk number 36: script_DEG.Rnw:545-553 (eval = FALSE)
###################################################
## ###Data: scatterDEG
## RPKM.S <- data.frame(RPKM_FC[unlist(List_RPKM_UporDown),])
## bayseq.S <- data.frame(bayseqDF [unlist(List_bayseqDF_UporDown),])
## edgeR.S <- data.frame(edgeDF[unlist(List_edgeR_UporDown),])
## NBPSeq.glm.S <- data.frame(NBPSeq.glmDF[unlist(List_NBPSeq.glmDF_UporDown),])
## deseq2.S <- data.frame(deseq2DF[unlist(List_DESeq2_UporDown),]) 
## TSPM.S <- data.frame(TSPMDF [unlist(List_TSPM_UporDown),])
## NBPSeq.nbp.S <- data.frame(NBPSeq.nbpDF[unlist(List_NBPSeq.nbpDF_UporDown),])


###################################################
### code chunk number 37: script_DEG.Rnw:556-578 (eval = FALSE)
###################################################
## ##logFC
## RPKM.logFC <- data.frame(rownames(RPKM.S), RPKM.S$"AP3.67_AG.67_logFC", row.names=1) 
## edgeR.logFC <- data.frame(rownames(edgeR.S), edgeR.S$"AP3.67.AG.67_logFC", row.names=1)
## deseq2.logFC <- data.frame(rownames(deseq2.S), deseq2.S$"AP3.67.AG.67_logFC" , row.names=1)
## NBPSeq.glm.logFC <- data.frame(rownames(NBPSeq.glm.S), NBPSeq.glm.S$"AP3.67_AG.67_logFC" , row.names=1)
## NBPSeq.nbp.logFC <- data.frame(rownames(NBPSeq.nbp.S), NBPSeq.nbp.S$"AP3.67_AG.67_logFC" , row.names=1)
## TSPM.logFC <- data.frame(rownames(TSPM.S), TSPM.S$"AP3.67_AG.67_logFC", row.names=1)  
##       
## scatterDEG1 <- merge(edgeR.logFC, deseq2.logFC, by='row.names', all=TRUE); rownames(scatterDEG1) <- scatterDEG1[[1]] ; scatterDEG1 <- scatterDEG1 [,-1]
## scatterDEG2 <- merge(scatterDEG1, RPKM.logFC , by='row.names', all=TRUE); rownames(scatterDEG2) <- scatterDEG2[[1]] ; scatterDEG2 <- scatterDEG2 [,-1]
## scatterDEG3 <- merge(scatterDEG2, NBPSeq.glm.logFC, by='row.names', all=TRUE); rownames(scatterDEG3) <- scatterDEG3[[1]] ; scatterDEG3 <- scatterDEG3 [,-1]
## scatterDEG4 <- merge(scatterDEG3, NBPSeq.nbp.logFC, by='row.names', all=TRUE); rownames(scatterDEG4) <- scatterDEG4[[1]] ; scatterDEG4 <- scatterDEG4 [,-1]
## scatterDEG5 <- merge(scatterDEG4, TSPM.logFC, by='row.names', all=TRUE); rownames(scatterDEG5) <- scatterDEG5[[1]] ; scatterDEG5 <- scatterDEG5 [,-1]
## scatterDEG_logFC <- scatterDEG5
## colnames(scatterDEG_logFC) <- c("edgeR", "DESeq2", "RPKM", "NBPSeq.glm", "NBPSeq.nbp", "TSPM")
## scatterDEG_logFC[is.na(scatterDEG_logFC)] <- 0 
## 
## ###Scatterplot
## scatterDEG_logFC[1:4,]
## pdf("./results/Scatterplot_AP367_AG67_logFC.pdf")
## pairs(scatterDEG_logFC, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main="Scatterplot Matrix log_FC")
## dev.off()


###################################################
### code chunk number 38: script_DEG.Rnw:587-610 (eval = FALSE)
###################################################
## ##FDR   
## bayseqFDR <- data.frame(rownames(bayseq.S), bayseq.S$"AP3.67_AG.67_FDR", row.names=1) 
## edgeRFDR <- data.frame(rownames(edgeR.S), edgeR.S$"AP3.67.AG.67_FDR", row.names=1)
## deseq2FDR <- data.frame(rownames(deseq2.S), deseq2.S$"AP3.67.AG.67_FDR" , row.names=1)
## NBPSeq.glmFDR <- data.frame(rownames(NBPSeq.glm.S), NBPSeq.glm.S$"AP3.67_AG.67_FDR" , row.names=1)
## NBPSeq.nbpFDR <- data.frame(rownames(NBPSeq.nbp.S), NBPSeq.nbp.S$"AP3.67_AG.67_FDR" , row.names=1)
## TSPMFDR <- data.frame(rownames(TSPM.S), TSPM.S$"AP3.67_AG.67_FDR", row.names=1)  
##       
## scatterDEG1 <- merge(edgeRFDR, deseq2FDR, by='row.names', all=TRUE); rownames(scatterDEG1) <- scatterDEG1[[1]] ; scatterDEG1 <- scatterDEG1 [,-1]
## scatterDEG2 <- merge(scatterDEG1, bayseqFDR, by='row.names', all=TRUE); rownames(scatterDEG2) <- scatterDEG2[[1]] ; scatterDEG2 <- scatterDEG2 [,-1]
## scatterDEG3 <- merge(scatterDEG2, NBPSeq.glmFDR, by='row.names', all=TRUE); rownames(scatterDEG3) <- scatterDEG3[[1]] ; scatterDEG3 <- scatterDEG3 [,-1]
## scatterDEG4 <- merge(scatterDEG3, NBPSeq.nbpFDR, by='row.names', all=TRUE); rownames(scatterDEG4) <- scatterDEG4[[1]] ; scatterDEG4 <- scatterDEG4 [,-1]
## scatterDEG5 <- merge(scatterDEG4, TSPMFDR, by='row.names', all=TRUE); rownames(scatterDEG5) <- scatterDEG5[[1]] ; scatterDEG5 <- scatterDEG5 [,-1]
##       
## scatterDEG_FDR <- scatterDEG5
## colnames(scatterDEG_FDR) <- c("edgeR", "DESeq2", "baySeq", "NBPSeq.glm", "NBPSeq.nbp", "TSPM")
## scatterDEG_FDR[is.na(scatterDEG_FDR)] <- 1
## 
## ###Scatterplot
## scatterDEG_FDR[1:4,]
## pdf("./results/Scatterplot_AP367_AG67_FDR.pdf")
## pairs(scatterDEG_FDR, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main="Scatterplot Matrix FDR")
## dev.off()


###################################################
### code chunk number 39: script_DEG.Rnw:626-640 (eval = FALSE)
###################################################
## ##Data
## data.ROC <- scatterDEG5
## scatterDEG1[is.na(scatterDEG1)] <- 1; 
## data.ROC <- merge(data.ROC, tmp[1], by='row.names', all=TRUE ); rownames(data.ROC) <- data.ROC[[1]] ; data.ROC <- data.ROC [,-1]
## colnames(data.ROC) <- c("edgeR", "DESeq2", "Bayseq", "NBPSeq.glm", "NBPSeq.nbp", "TSPM", "Common")
## write.table(data.ROC, "results/data_ROC.xls", col.names=NA, quote=FALSE, sep="\t")
##       
## data.class <- data.ROC
## data.class[!is.na(data.class)] <- 0
## data.class[is.na(data.class)] <- 1
## write.table(data.class, "results/data.class.xls", col.names=NA, quote=FALSE, sep="\t")
##       
## data.ROC1 <- data.ROC
## data.ROC1[,1:7][is.na(data.ROC1[,1:7])] <- 1 


###################################################
### code chunk number 40: script_DEG.Rnw:643-672 (eval = FALSE)
###################################################
## ### ROCR
## pdf("./results/ROC.pdf")
## pred1 <- prediction(data.ROC1$Common, data.class$edgeR)
## perf1 <- performance(pred1, "tpr", "fpr")
## plot(perf1, avg= "threshold", col="black", lty=4,lwd= 2, main= "ROC curves compare with edgeR")
## par(new = TRUE)
## pred2 <- prediction(data.ROC1$Common, data.class$Bayseq)
## perf2 <- performance(pred2, "tpr", "fpr")
## plot(perf2, avg= "threshold", lty=4, lwd= 2, col="dodgerblue4")
## par(new = TRUE)
## pred3 <- prediction(data.ROC1$Common, data.class$NBPSeq.glm)
## perf3 <- performance(pred3, "tpr", "fpr")
## plot(perf3, avg= "threshold",  lty=4,lwd= 2, col="darkgreen")
## par(new = TRUE)
## pred4 <- prediction(data.ROC1$Common, data.class$NBPSeq.nbp)
## perf4 <- performance(pred4, "tpr", "fpr")
## plot(perf4, avg= "threshold",  lty=4,lwd= 2, col="darkviolet")
## par(new = TRUE)
## pred5 <- prediction(data.ROC1$Common, data.class$TSPM)
## perf5 <- performance(pred5, "tpr", "fpr")
## plot(perf5, avg= "threshold", lty=4, lwd= 2, col="firebrick")
## par(new = TRUE)
## pred6 <- prediction(data.ROC1$Common, data.class$DESeq2)
## perf6 <- performance(pred6, "tpr", "fpr")
## plot(perf6, avg= "threshold",  lty=4,lwd= 2, col="deeppink1")
## legend("bottomright", legend=c("edgeR / 1", "BaySeq / 0.749", "NBPSeq.glm / 0.533", "NBPSeq.nbp / 0.872 ", "TSPM / 0.623", "DESeq2 / 0.974"),
##        col=c("black", "dodgerblue4", "darkgreen", "darkviolet", "firebrick", "deeppink1"), lty=4, lwd=3, bty="n")
##  
## dev.off()


###################################################
### code chunk number 41: script_DEG.Rnw:676-684 (eval = FALSE)
###################################################
## ##calculate AUC
## auc.tmp1 <- performance(pred1,"auc"); auc1 <- as.numeric(auc.tmp1@y.values)
## auc.tmp2 <- performance(pred2,"auc"); auc2 <- as.numeric(auc.tmp2@y.values) 
## auc.tmp3 <- performance(pred3,"auc"); auc3 <- as.numeric(auc.tmp3@y.values)
## auc.tmp4 <- performance(pred4,"auc"); auc4 <- as.numeric(auc.tmp4@y.values) 
## auc.tmp5 <- performance(pred5,"auc"); auc5 <- as.numeric(auc.tmp5@y.values)
## auc.tmp6 <- performance(pred6,"auc"); auc6 <- as.numeric(auc.tmp6@y.values) 
## 


###################################################
### code chunk number 42: script_DEG.Rnw:694-699 (eval = FALSE)
###################################################
## dFDR <- cor(data.ROC1, method="spearman")
## hc1 <- hclust(as.dist(1-dFDR))
## pdf("./results/methodsCorrelation_FDR.pdf")
## plot.phylo(as.phylo(hc1), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE,  main="Correlation FDR")
## dev.off()


###################################################
### code chunk number 43: script_DEG.Rnw:707-712 (eval = FALSE)
###################################################
## dlogFC <- cor(scatterDEG_logFC, method="spearman")
## hc2 <- hclust(as.dist(1-dlogFC))
## pdf("./results/methodsCorrelation_logFC.pdf")
## plot.phylo(as.phylo(hc2), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, main="Correlation logFC")
## dev.off()


###################################################
### code chunk number 44: sessionInfo
###################################################
toLatex(sessionInfo())


