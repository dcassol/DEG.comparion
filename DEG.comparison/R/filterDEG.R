##-------------------------------------------------------------------
##   Name: filterDEGnew
##   Function modified from 'systemPipeR' - (Thomas Girke. systemPipeR: NGS workflow and report generation environment, 28 June 2014. URL https://github.com/tgirke/systemPipeR.)    
##      Date: May 2015
##      Daniela Cassol
##-------------------------------------------------------------------

filterDEGnew <- function (degDF, filter, plot = TRUE, method) 
{
  pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
  log2FC <- degDF[, grep("_logFC$", colnames(degDF)), drop = FALSE]
  pf <- pval <= filter["FDR"]/100 & (log2FC >= log2(filter["Fold"]) | 
                                       log2FC <= -log2(filter["Fold"]))
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                                     x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  pf <- pval <= filter["FDR"]/100 & log2FC >= log2(filter["Fold"])
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                               x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  pf <- pval <= filter["FDR"]/100 & log2FC <= -log2(filter["Fold"])
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                                 x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  df <- data.frame(Comparisons = names(DEGlistUPorDOWN), Counts_Up_or_Down = sapply(DEGlistUPorDOWN, 
                                                                                    length), Counts_Up = sapply(DEGlistUP, length), Counts_Down = sapply(DEGlistDOWN, 
                                                                                                                                                         length))
  resultlist <- list(UporDown = DEGlistUPorDOWN, Up = DEGlistUP, 
                     Down = DEGlistDOWN, Summary = df)
  if (plot == TRUE) {
    mytitle <- paste("DEG Counts (", names(filter)[1], ": ", 
                     filter[1], " & ", names(filter)[2], ": ", filter[2], 
                     "%) ", method, sep = "") #just acrescent the method! 
    df_plot <- data.frame(Comparisons = rep(as.character(df$Comparisons), 
                                            2), Counts = c(df$Counts_Up, df$Counts_Down), Type = rep(c("Up", 
                                                                                                       "Down"), each = length(df[, 1])))
    p <- ggplot(df_plot, aes(Comparisons, Counts, fill = Type)) + 
      geom_bar(position = "stack", stat = "identity") + 
      coord_flip() + theme(axis.text.y = element_text(angle = 0, 
                                                      hjust = 1)) + ggtitle(mytitle)
    print(p)
  }
  return(resultlist)
}

filterDEG_logFC <- function (degDF, filter, plot = TRUE, method) 
{
  log2FC <- degDF[, grep("_logFC$", colnames(degDF)), drop = FALSE]
  pf <- (log2FC >= log2(filter["Fold"]) | log2FC <= -log2(filter["Fold"]))
  colnames(pf) <- gsub("_logFC", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                                     x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  pf <- log2FC >= log2(filter["Fold"])
  colnames(pf) <- gsub("_logFC", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUP <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                               x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  pf <- log2FC <= -log2(filter["Fold"])
  colnames(pf) <- gsub("_logFC", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                                 x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  df <- data.frame(Comparisons = names(DEGlistUPorDOWN), 
                   Counts_Up_or_Down = sapply(DEGlistUPorDOWN, length), 
                   Counts_Up = sapply(DEGlistUP, length), Counts_Down = sapply(DEGlistDOWN, length))
  resultlist <- list(UporDown = DEGlistUPorDOWN, Up = DEGlistUP, 
                     Down = DEGlistDOWN, Summary = df)
  library(ggplot2)
  if (plot == TRUE) {
    mytitle <- paste("DEG Counts (", names(filter)[1],  ": ", filter[1],
                     "%) ", method, sep = "")
    df_plot <- data.frame(Comparisons = rep(as.character(df$Comparisons), 2), 
                          Counts = c(df$Counts_Up, df$Counts_Down), Type = rep(c("Up", "Down"), 
                                                                               each = length(df[, 1])))
    p <- ggplot(df_plot, aes(Comparisons, Counts, fill = Type)) + 
      geom_bar(position = "stack", stat = "identity") + 
      coord_flip() + theme(axis.text.y = element_text(angle = 0, 
                                                      hjust = 1)) + ggtitle(mytitle)
    print(p)
  }
  return(resultlist)
}

filterDEG_FDR <- function (degDF, filter, plot = TRUE, method) 
{
  pval <- degDF[, grep("_FDR$", colnames(degDF)), drop = FALSE]
  pf <- pval <= filter["FDR"]/100 
  colnames(pf) <- gsub("_FDR", "", colnames(pf))
  pf[is.na(pf)] <- FALSE
  DEGlistUPorDOWN <- sapply(colnames(pf), function(x) rownames(pf[pf[, 
                                                                     x, drop = FALSE], , drop = FALSE]), simplify = FALSE)
  
  df <- data.frame(Comparisons = names(DEGlistUPorDOWN), Counts_Up_or_Down = sapply(DEGlistUPorDOWN, length))
  resultlist <- list(UporDown = DEGlistUPorDOWN, Summary = df)
  library(ggplot2)
  if (plot == TRUE) {
    mytitle <- paste("DEG Counts (", names(filter)[1],  ": ", filter[1],
                     "%) ", method, sep = "")
    df_plot <- data.frame(Comparisons = rep(as.character(df$Comparisons), 
                                            2), Counts = c(df$Counts_Up_or_Down), Type = rep(c("Up or Down"), each = length(df[, 1])))
    p <- ggplot(df_plot, aes(Comparisons, Counts, fill = Type)) + 
      geom_bar(position = "stack", stat = "identity") + coord_flip() +
      theme(axis.text.y = element_text(angle = 0, hjust = 1)) + ggtitle(mytitle)
    print(p)
  }
  return(resultlist)
}

