##-------------------------------------------------------------------
##   Name: run_DEGs 
##   Date: June 2015
##   Daniela Cassol 
##-------------------------------------------------------------------

#######################################################################
## DEG1: Run RPKM: Simple Fold Change Method ##
#######################################################################
run_RPKM <- function (counts, mycomp1, mycomp2){
  sample <- mycomp1$Sample
  group <- mycomp1$group
  factor <- mycomp1$Factor
  myfun <- mean
  mycomp1 <- tapply(sample, group, list)
  names(mycomp1) <- sapply(mycomp1, paste, collapse="_")
  mymean <- sapply(mycomp1, function(x) apply(counts[, x, drop=FALSE], 1, myfun))
  colnames(mymean) <- factor
  mylog <- sapply(names(mycomp2), function(x) log2(mymean[,mycomp2[[x]][1]] / mymean[,mycomp2[[x]][2]]))
  is.na(mylog) <- sapply(mylog, is.infinite)
  colnames(mylog) <- paste(paste(colnames(mylog), sep="_", "logFC"))
  return(mylog)
}

#usage: 
# Comp1 <- list(Factor=(Reduce(union, targets$Factor)), Sample=c(colnames(rpkmDFeByg)), group=c(1,1,2,2))
# Comp2 <- list(AP1.4_AP1.67=c("AP1.4", "AP1.67"))
# RPKM_FC <- run_RPKM (rpkmDFeByg, Comp1, Comp2)

################################################################
###DEG4: Run baySeq
################################################################
run_BaySeq <- function (counts, mycomp3, number){
  if(require("parallel")) cl <- makeCluster(8) else cl <- NULL
  BaySeqDF <- data.frame(row.names = rownames(counts))
  counts <- as.matrix(counts) # SO IMPORTANT!!!!!!!!!
  for(i in seq(along=mycomp3))
  {
    if(length(mycomp3[[i]])==2){
      groups <- list(NDE = c(1,1), DE = c(1,2))
      replicates <-  c("simA", "simB")}
    if(length(mycomp3[[i]])==3){
      groups <- list(NDE = c(1,1,1), DE = c(1,2,2))
      replicates <-  c("simA", "simB", "simB")}
    if(length(mycomp3[[i]])==4){
      groups <- list(NDE = c(1,1,1,1), DE = c(1,1,2,2)) 
      replicates <-  c("simA", "simA", "simB", "simB")}
    
    CD <- new("countData", data = counts[ , mycomp3[[i]]], replicates = replicates, groups = groups)
    print(colnames(counts[ , mycomp3[[i]]]))
    libsizes(CD) <- getLibsizes(CD)
    CD@annotation <- data.frame(row.names = rownames(counts[ , mycomp3[[i]]]))
    CD <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl) 
    CD <- getLikelihoods(CD, pET = 'BIC', cl = cl)
    deg <- as.data.frame(cbind(topCounts(CD, group="DE", number=number))) #creates a data.frame the same in the original countDFeByg 
    colnames(deg)[colnames(deg) %in% c("FDR.DE")] <- c("FDR")
    colnames(deg) <- paste(paste(names(mycomp3[i]), colnames(deg), sep="_"))
    BaySeqDF <- cbind(BaySeqDF, deg[rownames(BaySeqDF), ])
  }
  return(BaySeqDF)
  showConnections(all = FALSE)
}

# Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"))
# bayseqDF <- run_BaySeq(countDFeByg, Comp3, number=27416)

################################################################
# myFUNCTION ###nb.glm.test 
################################################################
run_NBPSeq_glm <- function (counts, mycomp3){
  NBPSeqDF <- data.frame(row.names = rownames(counts))
  counts <- as.matrix(counts) 
  beta0 = c(NA, 0)
  for(i in seq(along=mycomp3)){
    
    if(length(mycomp3[[i]])==2){
      grp.ids = as.factor(c(1,2))
      x = model.matrix(~grp.ids)
      stop("only sample with biological replicates")}
    if(length(mycomp3[[i]])==3){
      grp.ids = as.factor(c(1,2,2))
      x = model.matrix(~grp.ids)
      stop("only sample with biological replicates")}
    if(length(mycomp3[[i]])==4){
      grp.ids = as.factor(c(1,1,2,2))
      x = model.matrix(~grp.ids)}
    
    fit <- NBPSeq::nb.glm.test(counts=counts[ , mycomp3[[i]]], x, beta0)
    print(colnames(counts[ , mycomp3[[i]]]))
    deg <- as.data.frame(cbind(fit$data$counts, fit$test.results$HOA))
    colnames(deg)[colnames(deg) %in% c("statistic", "q.values")] <- c("logFC", "FDR")
    colnames(deg) <- paste(paste(names(mycomp3[i]), colnames(deg), sep="_"))
    NBPSeqDF <- cbind(NBPSeqDF, deg[rownames(NBPSeqDF), ])  
  }
  return(NBPSeqDF)
  showConnections(all = FALSE)
}

#################################################################
## myFUNCTION # nbp.test 
#################################################################
run_NBPSeq_nbp <- function (counts, mycomp3){
  NBPSeqDF <- data.frame(row.names = rownames(counts))
  counts <- as.matrix(counts) 
  grp1 = 1
  grp2 = 2
  set.seed(999)
  for(i in seq(along=mycomp3)){
    
    if(length(mycomp3[[i]])==2){
      grp.ids = as.factor(c(1,2)) 
      stop("only sample with biological replicates")}
    if(length(mycomp3[[i]])==3){
      grp.ids = as.factor(c(1,2,2)) 
      stop("only sample with biological replicates")}
    if(length(mycomp3[[i]])==4){
      grp.ids = as.factor(c(1,1,2,2))}
    norm.factors = estimate.norm.factors(counts=counts[ ,mycomp3[[i]]], lib.sizes = colSums(counts[ ,mycomp3[[i]]]), method = "AH2010") 
    res <- NBPSeq::nbp.test(counts=counts[ ,mycomp3[[i]]], grp.ids, grp1, grp2, 
                            norm.factors = norm.factors, model.disp = "NBQ", lib.sizes = colSums(counts[ ,mycomp3[[i]]]), 
                            print.level = 3) 
    print(colnames(counts[ , mycomp3[[i]]]))
    deg <- as.data.frame(cbind(res$log.fc, res$p.values, res$q.values))
    colnames(deg) <- c("logFC", "p.values", "FDR")
    colnames(deg) <- paste(paste(names(mycomp3[i]), colnames(deg), sep="_"))
    NBPSeqDF <- cbind(NBPSeqDF, deg[rownames(NBPSeqDF), ])  
  }
  return(NBPSeqDF)
  showConnections(all = FALSE)
}

run_TSPM <- function (counts, mycomp3){
  TSPM_DF <- data.frame(row.names = rownames(counts))
  counts <- as.matrix(counts) 
  for(i in seq(along=mycomp3)){
    if(length(mycomp3[[i]])==2){
      x1= factor(c(1, 2))
      x0= (c(1, 1))
      stop("only sample with biological replicates")}
    if(length(mycomp3[[i]])==3){
      x1= factor(c(1, 2, 2))
      x0= (c(1, 1, 1))
      stop("only sample with biological replicates")}
    if(length(mycomp3[[i]])==4){
      x1= factor(c(1, 1, 2, 2))
      x0= (c(1, 1, 1, 1))}
    result <- TSPM(counts=counts[ , mycomp3[[i]]], x1, x0, lib.size = apply((counts[ , mycomp3[[i]]]),2,sum))
    print(colnames(counts[ , mycomp3[[i]]]))
    deg <- data.frame(result$log.fold.change, result$pvalues, result$padj)
    colnames(deg) <- c("logFC", "p.values", "FDR")
    colnames(deg) <- paste(paste(names(mycomp3[i]), colnames(deg), sep="_"))
    TSPM_DF <- cbind(TSPM_DF, deg)
  }
  return(TSPM_DF)
}
