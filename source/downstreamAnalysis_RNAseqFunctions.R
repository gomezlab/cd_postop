#load required libraries
library(DESeq2) 
library(limma)
#library(sva)
library(ggrepel)
library(ggplot2) 
library(plotly)
library(plyr)
library(dplyr)
library(GGally)
library(tidyverse)
library(scales)
library(pheatmap)
library(reshape2)
#library(ComplexHeatmap)
library(circlize)
library(edgeR)
library(preprocessCore)
library(RColorBrewer) 
#library(M3C)
library(gplots)
library(VennDiagram)
#library(GeneOverlap)
#library(EnhancedVolcano)
#library(ggpubr)

#write some code to install all of the packages above if they are not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c('ggrepel', 'ggplot2', 'plotly', 'plyr', 'dplyr', 'GGally', 'tidyverse', 'scales', 'pheatmap', 'reshape2', 'circlize', 'RColorBrewer', 'gplots', 'VennDiagram'))
BiocManager::install(c('edgeR', 'preprocessCore', 'limma'))



# I load data into a data-frame-like structure called a tibble. This makes sure that each column data-type is properly specified to avoid downstream issues with plotting numerical data as characters. This function just converts the tibble back to a dataframe.
tibble2dataframe <- function(d) {
  r <- as.data.frame(d)
  rownames(r) <- r[, 1]
  r <- r[, -1]
  return(r)
}

# Differential analysis
differentialAnalysis <- function(counts, coldata, designFormula) {
	dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = designFormula)

	dds <- DESeq(dds) # Run differential analysis
	res <- results(dds) # Generate results table
	write.table(res[order(res$padj),], file="DESeq_results.txt", sep="\t", quote=FALSE) # Write results to text file.
}

# Note: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the formula and make sure the control level is the first level.
DESeqNormalize <- function(counts, coldata, designFormula) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = designFormula)
  
  dds <- estimateSizeFactors(dds)
  return(dds)
}

# if no "trans" argument is specified, variance stabilizing transformation is the default transformation.
# The use of "blind=FALSE" is to take the design formula into account when transforming the data. According to Mike, blinding the data to the design formula can result in an overestimation of varaiance. In my experience, this actually makes little difference to the downstream results, but it's best to include this as a precaution.
transformData <- function(dds, trans = "vst", blind = FALSE) {
  if (trans == "vst") {
    if (nrow(assay(dds)) < 1000) {
      vst <- varianceStabilizingTransformation(dds, blind)
      return(vst)
    }
    else {
      vst <- vst(dds, blind)
      return(vst)
    }
  }
  else if (trans == "rlog") {
    rld <- rlog(dds, blind)
    return(rld)
  }
}

# Using edgeR's cpm and proprocessCore's normalize quantiles to normalize and transform the data. This seems to be specified as an alternative method to process the raw ATAC data. Although it should produce similar results to DESeq norm/cst, this method is cited (specifically the normalize quantile portion).
quantileNorm <- function(countMatrix, priorCount, log = TRUE) {
  cpmLog <- cpm(countMatrix, log = log, prior.count = priorCount)
  
  cpmLog.QN <- as.data.frame(normalize.quantiles(cpmLog))
  rownames(cpmLog.QN) <- rownames(cpmLog); colnames(cpmLog.QN) <- colnames(cpmLog)
  return(cpmLog.QN)
}

# function to correct normalized/transformed counts after performing svaseq for the purposes of downstream analyses like pca of hierarchical clustering.
#'y' as the gene expresion matrix
#'mod' as the model matrix you sent to sva (the full model)
#'svs' as svobj$sv where svobj is the output from the sva or svaseq function
svaCorrect <- function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

# Limma function removeBatchEffect for batch correction. Batch and covariates can be specified.
#Batch2 is just in case there is another linear batch
#IECAve and TSS 

limmaCorrect <- function(transData, batch, design, covariate, batch2) {
  
  if(missing(design)) {
    design <- matrix(1,ncol(transData),1)
  }
  if(typeof(transData) == "S4") {
    transData <- assay(transData)
  }
  
  if(missing(batch2)) {
    if (missing(covariate)) {
      transData <- removeBatchEffect(x = transData, batch = batch, design = design)
      return(as.data.frame(transData))
    } else { # covariate specified
      transData <- removeBatchEffect(x = transData, batch = batch, covariates = covariate, design = design)
      return(as.data.frame(transData))
    }
  } else { #batch2 specificed
    if (missing(covariate)) {
      transData <- removeBatchEffect(x = transData, batch = batch, batch2 = batch2, design = design)
      return(as.data.frame(transData))
    } else { # covaraite specified
      transData <- removeBatchEffect(x = transData, batch = batch, batch2 = batch2, covariates = covariate, design = design)
      return(as.data.frame(transData))
    }
  }
}

# Function to fix rownames after filtering to avoid potential DESeq2 issues.
postFilterRownames <- function (dataframe, statusColumn=2) {
  
  # This whole section of code is to deal with getting the right colnames for DESeq2
  nibdCt <- 1
  cdCt <- 1
  coldatRownames <- vector('character')
  for (i in seq(from=1, nrow(dataframe), by=1)) {
    line <- dataframe[i,]
    coldatRownames <- c(coldatRownames, paste0(if (line[statusColumn]=="CD") "cd" else "nibd", if (line[statusColumn]=="CD") cdCt else nibdCt))
    if (line[statusColumn]=="CD") cdCt=cdCt+1 else nibdCt=nibdCt+1
  }
  
  #applying the new colnames to the 'passed' coldata 
  rownames(dataframe) <- coldatRownames
  return(dataframe)
  
}


# Filter function to isolate samples of interest based on some criteria in the coldata
filterData <- function(coldata, cts, filterString, matchingColumn="sample") {
  if (typeof(filterString) != "language") {
    stop("To make sure that the filter is parsed correctly, wrap the filterString in 'quote()'")
  }
  # COLDATA
  coldata.filter <- filter(coldata, !!filterString)
  # This is the fix the rownames incase this data is used for any DESeq related functions.
  rownames(coldata.filter) <- coldata.filter[,as.character(matchingColumn)]
  
  # CTS
  if (typeof(cts) == "S4") {cts.filter <- as.data.frame(assay(cts))} else {cts.filter <- cts}
  cts.filter <- cts.filter %>% dplyr::select(rownames(coldata.filter)) 
  
  return(list(coldata.filter,cts.filter))
}

coldata.filterData <- function(filteredData) {
  coldata <- as.data.frame(filteredData[1])
  return(coldata)
}

cts.filterData <- function (filteredData) {
  cts <- as.data.frame(filteredData[2], check.names=FALSE)
  return(cts)
}

# This function allows you to threshold data and pull out specific biotypes (e.g. lncRNA/protein-coding genes). Values must be set for thresholding.
# A normThres of 6 and sampleThres of 5 would mean that a gene is only retained in the output if it has a transformed count value of at least 6 across at least 5 samples. 
thresholdBiotype <- function(transData, normThres, sampleThres, biotype) {
  print(paste0("Number of genes in input data: ", nrow(transData)))
  
  if (!missing(biotype)) { 
  # biotype Filter
  transData$ensembl <- rownames(transData)
  transData <- join(transData, ens2gene, by=c("ensembl"))
  transData <- filter(transData, biotype == !!biotype)
  print(paste0("Number of ", biotype, " genes: ", nrow(transData)))
  rownames(transData) <- transData$ensembl
  transData$ensembl <- NULL; transData$biotype <- NULL; transData$gene_id <- NULL; transData$chr <- NULL; transData$start <- NULL; transData$stop <- NULL
  }
  
  if (!missing(normThres)) {
    # threshold
    idx <- rowSums(transData >= normThres) >= sampleThres
    transData <- transData[idx,]
  }
  print(paste0("Number of genes in output data: ", nrow(transData)))
  return(transData)
}

# Selecting genes based on decreaseing varaince across all samples. By default, this function will extract the top 1000 variable genes from the input.
varFilter <- function(transData, variableGenes=1000) {
  rv <- rowVars(transData)
  select <- order(rv, decreasing = TRUE)[seq_len(min(variableGenes, length(rv)))]
  transData <- transData[select, ]
}

# General function for PCA.
rnaSeqPCA <- function (object, coldata, plotLabs=FALSE, labs="sample", returnData=FALSE, title, color="status", manualColor="none", x="PC1", y="PC2", proportional=FALSE, savePlot=FALSE, saveTitle="PCAPlot", loadings = FALSE, loadingsTitle = "pcaLoadings.txt", outDir = "plots/pca", plotly=FALSE, legend = TRUE, plotRatio = 1, pointSize = 3, shape) {
  
  if (!all(rownames(coldata) %in% colnames(object))) {
    stop("rownames in coldata do not match the colnames in the counts")
  }
  
  pca <- prcomp(t(object))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  x_pc <- as.numeric(substr(x, 3, 3))
  y_pc <- as.numeric(substr(y, 3, 3))
  percentVar <- percentVar[c(x_pc,y_pc)]
  
  d <- data.frame(cbind(pca$x,coldata))
  
  p <- ggplot(data = d, aes_string(x = as.character(x), y = as.character(y), color = as.character(color), label = labs)) + 
    geom_point(size = pointSize) + 
    xlab(paste0(x,": ", percent(percentVar[1]), " Variance")) + 
    ylab(paste0(y,": ", percent(percentVar[2]), " Variance")) + 
    theme_minimal() +
    coord_fixed(ratio = plotRatio)
  if (!missing(shape)) {
    p <- ggplot(data = d, aes_string(x = as.character(x), y = as.character(y), color = as.character(color), label = labs, shape = as.character(shape))) + 
      geom_point(size = pointSize) + 
      xlab(paste0(x,": ", percent(percentVar[1]), " Variance")) + 
      ylab(paste0(y,": ", percent(percentVar[2]), " Variance")) + 
      theme_minimal() +
      coord_fixed(ratio = plotRatio)
  }
  
  if (loadings) {
    loadings <- as.data.frame(pca$rotation)
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
    loadingsOut <- loadings[,1:5]; loadingsOut$ensembl <- rownames(loadingsOut)
    loadingsOut <- join(loadingsOut, ens2gene, by="ensembl")
    write.table(loadingsOut, file = paste0(outDir,"/",loadingsTitle), sep="\t", quote=F, col.names = NA)
  }
  
  if (manualColor == "status") {p <- p + scale_color_manual(values = c("#0a9618", "#000000"))}
  else if (manualColor == "sex") {p <- p + scale_color_manual(values = c("#DC4AE7", "#DC7D00"))}
  else if (manualColor == "statusNA") {p <- p + scale_color_manual(values = c("#0a9618", "#000000", "#CBCBCB"))}
  else if (manualColor == "subtype") {p <- p + scale_color_manual(values = c("#0000FF", "#E60000", "#000000"))}
  else if (manualColor == "pcaExtremes") {p <- p + scale_color_manual(values = c("#6582cf", "#e36666", "#8f8f8f"))} 
  else if (manualColor == "CLsubtypes") {p <- p + scale_color_manual(values = c("#0BC0FE", "#EB8400"))} 
  else if (manualColor == "CLsubtypesNIBD") {p <- p + scale_color_manual(values = c("#000000", "#0BC0FE", "#EB8400"))}  
  
  if (!missing(title)) {p <- p + labs(title = title)}   
  if (plotLabs) {p <- p + geom_text_repel(size=2.5)}
  if (!legend) (p <- p + theme(legend.position = "none"))
  
  if (proportional) {
    print (p)
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
    aspectRatio <- as.numeric(percentVar[1])/as.numeric(percentVar[2])
    width <- 7*aspectRatio
    ggsave(paste0(outDir,"/",saveTitle,"_proportional.pdf"), device="pdf", height=7, width=width, useDingbats=FALSE)
  }
  
  if (savePlot) {
    print(p)
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE) 
    ggsave(paste0(outDir,"/",saveTitle,".pdf"), device="pdf", useDingbats=FALSE)
  }
  
  if (plotly) {return(ggplotly(p))}
  if (!returnData) {return(p)} else {print (p); return(d)}
  
}

#Function for multi PCA
rnaSeqMultiPCA <- function(object, coldata, color="status", savePlot = FALSE, saveTitle="PCAMultiPlot", outDir = "plots/pca") {
  
  pca <- prcomp(t(object))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  #x_pc <- as.numeric(substr(x, 3, 3))
  #y_pc <- as.numeric(substr(y, 3, 3))
  #percentVar <- percentVar[c(x_pc,y_pc)]
  
  d <- data.frame(cbind(pca$x,coldata))
  
  #### ggpairs ####
  ## Lines to sort out the axis labs for ggpairs ##
  colLabs <- vector()
  for (i in 1:5) {
    colLabs <- c(colLabs, paste0("PC", i, ": ", percent(percentVar[i])))
  }
  ## run ggpairs
  #plotTitle = paste0(datframe, "_", category)
  if (is.double(eval(parse(text=paste0("d", "$", color, "[1]"))))){
    print("here_cont")
    p <- ggpairs(d, columns = 1:5, 
                 columnLabels = colLabs, 
                 #title = plotTitle,
                 upper = list(continuous = wrap("points", alpha = 1), combo = "box"),
                 lower = list(discrete = wrap("points", alpha = 1, size=0.4), 
                              combo = wrap("dot", alpha = 1, size=0.4) ),
                 mapping = aes(color = get(color)))
    #p <- ggpairsColorAdjust
  }
  else {
    print("here_dis")
    p <- ggpairs(d, mapping = aes_string(color = color),
                 columns = 1:5, 
                 columnLabels = colLabs, 
                 #title = plotTitle,
                 upper = list(continuous = wrap("points", alpha = 1), combo = "box"),
                 lower = list(discrete = wrap("points", alpha = 1, size=0.4), 
                              combo = wrap("dot", alpha = 1, size=0.4) ),
                 legend = c(5,5))
    #p <- ggpairsColorAdjust(p, colors = c("#0a9618", "#000000"))
  }
  
  
  if (savePlot) {
    print(p)
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
    ggsave(paste0(outDir,"/",saveTitle,".pdf"), device="pdf")
  }
  else {
    return(p)
  }
}

# Function to adjust colors of plots in the multi pca Plot
# two colors are required for the continuous flag, and colors must match the number of factors for discrete values
# colors specfied in hex format as a character vector -> c("#FF0000", "#0000FF")
ggpairsColorAdjust <- function(plot, colors, continuous = FALSE) {
  #colours <- as.vector(colours)
  #print(colours)
  for(i in 1:plot$nrow) {
    for(j in 1:plot$ncol){
      if(continuous) {
        plot[i,j] <- plot[i,j] + 
        scale_color_gradient(high=colors[1], low=color[2])
      }
      else {
        plot[i,j] <- plot[i,j] + 
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors)  
      }
    }
  }
  plot$legend = 1
  #plot$title = plotTitle
  #plot$upper = "na"
  return(plot + theme(legend.position = "bottom"))
  #return(plot)
}

#Function for creating pca gifs
benchmarking <- function(thresholds, countMatrix, coldata, plotTitle = "var", delay = "1x20", gifTitle = "varThresholds", color = "status", barPlot=FALSE, manualColor="none", benchmark = "varGenes", variableGenes = 1000, sampleThres = 5, biotype="protein_coding", x = "PC1", y="PC2") {
  loop <- 10000
  thresholds <- as.vector(thresholds)
  countMatrix <- as.data.frame(countMatrix)
  dir.create(path = "tmp", recursive = TRUE)
  for (i in seq(thresholds[1], thresholds[2], thresholds[3])) {
    
    if (benchmark == "varGenes") {
      varCountMatrix <- varFilter(countMatrix, variableGenes = i)
    }
    else if (benchmark == "threshold") {
      thresCountMatrix <- thresholdBiotype(countMatrix, normThres = i, sampleThres = sampleThres, biotype=biotype)
      
      varCountMatrix <- varFilter(thresCountMatrix, variableGenes = variableGenes)
    }
    else{
      stop("benchmark should be specifed as either 'threshold' or 'varGenes'.")
    }
    
    object <- varCountMatrix
    coldata <- coldata
    
    pca <- prcomp(t(object))
    percentVar <- percent(pca$sdev^2/sum(pca$sdev^2))
    
    x_pc <- as.numeric(substr(x, 3, 3))
    y_pc <- as.numeric(substr(y, 3, 3))
    percentVar <- percentVar[c(x_pc,y_pc)]
    
    d <- data.frame(cbind(pca$x,coldata))
    if (!barPlot) {
      p <- ggplot(data = d, aes_string(x = as.character(x), y = as.character(y), color = as.character(color), label = "sample")) + 
        geom_point(size = 3) + xlab(paste0(x,": ", percentVar[1], " Variance")) + 
        ylab(paste0(y,": ", percentVar[2], " Variance")) + theme_minimal()  + 
        if (benchmark == "threshold") {
          labs(title = paste0(plotTitle,"_", nrow(object), "genes: ",i))
        }
        else {
          labs(title = paste0(plotTitle,": ",i))
        }

      if (manualColor == "status") {p <- p + scale_color_manual(values = c("#0a9618", "#000000"))}
      
      print(p)
    }
    else {
      percentVar <- (pca$sdev^2/sum(pca$sdev^2))*100
      percentVar <- as.data.frame(percentVar)
      percentVar$pcs <- as.numeric(rownames(percentVar))
      percentVar <- percentVar[1:10,]
      
      p <- ggplot(data=percentVar, aes(x=pcs, y=percentVar)) +
        geom_bar(stat="identity", fill="steelblue")+
        ylim(0,50) +
        theme_minimal() +
        theme(axis.text.x=element_blank())
      print(p)
    }
    
    ggsave(filename = paste0("tmp/",loop,"thresholds",i,".png"), device="png", dpi=150, scale=0.9)
    loop=loop+1
  }
  dir.create(path = "plots/thresholding", recursive = TRUE)
  system(paste0("convert -delay ",delay," tmp/*.png plots/thresholding/",gifTitle,"_",thresholds[1],"_",thresholds[2],".gif"))
  system("rm -rf tmp")
}

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  genesV2_unique <- unique(genesV2)
  colnames(genesV2_unique) <- c("gene","gene_id")
  
  #return(humanx)
  return(genesV2_unique)
}

#Haber or Comelli List overlap - Currently set up for a list with ensembl id in the rownames.
overlapAnalysis <- function(data, overlap="haber", annotatedList=T, plot=T, outDir = "plots/overlap", savePlot = T, saveTitle = paste0(overlap,"Overlap"), annoOut=T) {
  #setting up data for overlap
  if (overlap == "haber") {
    haberList <- read.csv("overlapAnalysis/haber_assignment_list.csv")
    
    mouseConvert <- convertMouseGeneList(haberList$gene)
    haberList <- join(haberList, mouseConvert, by="gene")
    haberList <- join(haberList, ens2gene, by = "gene_id")
    haberList <- na.omit(haberList)

    
    markerList <- haberList
    colnames(markerList) <- c("type", "gene", "gene_id", "ensembl", "biotype")
    
  } else if (overlap == "comelli") {
    comelliMarkers <- read.table("overlapAnalysis/190917_comelliMarkers_v2.txt", header = TRUE); colnames(comelliMarkers) <- c("gene_id", "markerLocation")
    
    comelliMarkers <- join(comelliMarkers, ens2gene, by="gene_id")
    comelliMarkers <- na.omit(comelliMarkers) 
    
    markerList <- comelliMarkers
    colnames(markerList) <- c("gene_id", "type", "ensembl", "biotype")
    
  } else if (overlap == "weiser") {
    weiserMarkers <- read.table("overlapAnalysis/190926_Weiser_ileumColon.txt", header = TRUE)
    
    weiserMarkers <- join(weiserMarkers, ens2gene, by="gene_id")
    weiserMarkers <- na.omit(weiserMarkers) 
    
    markerList <- weiserMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else if (overlap == "smillie") {
    smillieMarkers <- read.table("overlapAnalysis/190926_Smillie_epithelial.txt", header = TRUE)
    
    smillieMarkers <- join(smillieMarkers, ens2gene, by="gene_id")
    smillieMarkers <- na.omit(smillieMarkers) 
    
    markerList <- smillieMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else if (overlap == "subtypes") {
    subtypeMarkers <- read.table("overlapAnalysis/191004_colonSubtypes.txt", header = TRUE)
    
    subtypeMarkers <- join(subtypeMarkers, ens2gene, by="gene_id")
    subtypeMarkers <- na.omit(subtypeMarkers) 
    
    markerList <- subtypeMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else if (overlap == "wangAll") {
    wangMarkers <- read.table("overlapAnalysis/191204_Wang_allSegments.txt", header = TRUE) 
    
    wangMarkers <- join(wangMarkers, ens2gene, by="gene_id")
    wangMarkers <- na.omit(wangMarkers) 
    
    markerList <- wangMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else if (overlap == "wangColon") {
    wangMarkers <- read.table("overlapAnalysis/191204_Wang_colon.txt", header = TRUE)
    
    wangMarkers <- join(wangMarkers, ens2gene, by="gene_id")
    wangMarkers <- na.omit(wangMarkers) 
    
    markerList <- wangMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else if (overlap == "wangRectum") {
    wangMarkers <- read.table("overlapAnalysis/191204_Wang_rectum.txt", header = TRUE)
    
    wangMarkers <- join(wangMarkers, ens2gene, by="gene_id")
    wangMarkers <- na.omit(wangMarkers) 
    
    markerList <- wangMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else if (overlap == "wangIleum") {
    wangMarkers <- read.table("overlapAnalysis/191204_Wang_ileum.txt", header = TRUE)
    
    wangMarkers <- join(wangMarkers, ens2gene, by="gene_id")
    wangMarkers <- na.omit(wangMarkers) 
    
    markerList <- wangMarkers
    colnames(markerList) <- c("type", "gene_id", "ensembl", "biotype")
  } else {
    stop("the 'overlap' flag must be set to 'comelli', 'haber','weiser','smillie','subtypes','wangAll','wangColon','wangRectum', or 'wangIleum'")
  }
  
  annoList <- data.frame()
  plotCount <- vector()
  for (i in unique(markerList$type)){
    typeList <- as.data.frame(filter(markerList, type == i))
    typeNum <- nrow(typeList)
    #print(typeList)
    #print(typeNum)
    
    inputMatch <- typeList[as.character(typeList$ensembl) %in% rownames(data),]
    inputMatch <- na.omit(inputMatch)
    #print(inputMatch)
    
    annoList <- rbind(annoList, inputMatch)
    
    matchNumType <- nrow(inputMatch[!duplicated(inputMatch$gene),])
    matchNumData <- nrow(inputMatch[!duplicated(inputMatch$gene_id),])
    matchProp <- matchNumType/typeNum
    #print(matchProp)
    dataProp <- matchNumData/nrow(data)
    
    plotCount <- c(plotCount, dataProp)
    #print(plotCount)
    
    print(paste0(i, " percent Match: ", percent(matchProp), " - ", matchNumType, "/", typeNum, "; ", percent(dataProp), " - ", matchNumData, "/", nrow(data)))
    
  }
  colnames(annoList) <- colnames(markerList)
  
  if (plot) {
    
    plotData <- data.frame(Type=unique(markerList$type), Proportion=plotCount)
    #print(plotData)
    
    p <- ggplot(data=plotData, aes(x=Type, y=Proportion)) + ylim(0,1) +
      geom_bar(stat = "identity", fill="black") + theme_minimal() #+ ggsave(paste0(overlap,"Overlap.pdf"))
    if (overlap == "smillie") { p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))}
    print(p)
    
    if (savePlot) {
      dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
      ggsave(paste0(outDir,"/",saveTitle,".pdf"), device="pdf")
    }
    
    p_zoom <- ggplot(data=plotData, aes(x=Type, y=Proportion)) + 
      geom_bar(stat = "identity", fill="black") + theme_minimal() #+ ggsave(paste0(overlap,"Overlap_zoom.pdf"))
    if (overlap == "smillie") { p_zoom <- p_zoom + theme(axis.text.x = element_text(angle = 90, hjust = 1))}
    print(p_zoom)
    
    if (savePlot) {
      dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
      ggsave(paste0(outDir,"/",saveTitle,"_zoom.pdf"), device="pdf")
    }
    
  }
  
  if (annotatedList) {
    annoList_tmp <- data.frame(ensembl = rownames(data))
    annoList_tmp <- join(annoList_tmp, ens2gene, by="ensembl")
    annoList_tmp <- join(annoList_tmp, annoList, by="ensembl")
    if (annoOut) {write.csv(annoList_tmp, file=paste0(outDir, 
            "/", saveTitle, "_annoList.csv"), quote=FALSE)}
    return(annoList_tmp)

  }

}

#### tSNE wrapper
ben_tSNE <- function (mydata, K = FALSE, labels = FALSE, perplex = 15, printres = FALSE, 
          seed = FALSE, axistextsize = 30, legendtextsize = 30, dotsize = 6, 
          textlabelsize = 4) 
{
  if (seed != FALSE) {
    set.seed(seed)
  }
  if (K == FALSE && labels == FALSE) {
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, theta = 0, pca_center = FALSE, normalize = FALSE, eta = 10,
                         perplexity = perplex, verbose = TRUE, max_iter = 5000)
    scores <- data.frame(tsne$Y)
    p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
      geom_point(aes(colour = factor(rep(1,ncol(mydata)))), size = dotsize) + theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text.y = element_text(size = axistextsize, colour = "black"), 
            axis.text.x = element_text(size = axistextsize, colour = "black"), 
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) + 
      scale_colour_manual(values = c(`1` = "sky blue"))
    if (printres == TRUE) {
      message("printing tSNE to current directory...")
      png("TSNEpriorclustering.png", height = 20, width = 22, 
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (K != FALSE && labels == FALSE) {
    res <- mydata
    mydata <- res$realdataresults[[K]]$ordered_data
    annon <- res$realdataresults[[K]]$ordered_annotation
    annon$id <- row.names(annon)
    annon <- annon[match(colnames(mydata), annon$id), ]
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, theta = 0, pca_center = FALSE, normalize = FALSE, eta = 10,
                         perplexity = perplex, verbose = TRUE, max_iter = 5000)
    scores <- data.frame(tsne$Y)
    p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
      geom_point(aes(colour = factor(annon$consensuscluster)), size = dotsize) + theme_bw() + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
            axis.text.y = element_text(size = axistextsize, colour = "black"), 
            axis.text.x = element_text(size = axistextsize, colour = "black"), 
            axis.title.x = element_text(size = axistextsize), 
            axis.title.y = element_text(size = axistextsize), 
            legend.title = element_text(size = legendtextsize), 
            legend.text = element_text(size = legendtextsize)) + guides(colour = guide_legend(title = "Cluster"))
    if (printres == TRUE) {
      message("printing tSNE to current directory...")
      png("TSNEpostclustering.png", height = 20, width = 26, 
          units = "cm", res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  else if (K == FALSE && labels != FALSE) {
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, theta = 0, pca_center = FALSE, normalize = FALSE, eta = 10,
                         perplexity = perplex, verbose = TRUE, max_iter = 5000)
    scores <- data.frame(tsne$Y)
    p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
      geom_point(aes(colour = labels), size = dotsize) + theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text.y = element_text(size = axistextsize, colour = "black"), 
            axis.text.x = element_text(size = axistextsize, colour = "black"), 
            axis.title.x = element_text(size = axistextsize), 
            axis.title.y = element_text(size = axistextsize), 
            legend.title = element_text(size = legendtextsize), 
            legend.text = element_text(size = legendtextsize)) #+ guides(colour = guide_legend(title = "Cluster"))
    if (printres == TRUE) {
      message("printing tSNE to current directory...")
      png("TSNElabeled.png", height = 20, width = 26, units = "cm", 
          res = 900, type = "cairo")
      print(p)
      dev.off()
    }
  }
  return(p)
}

# function to scale covariates for DESeq2 (if needed upon DESeq warning message)
centerAndScale <- function(x) (x - mean(x)) / sd(x)

# Function to generate a basic pairwise venn diagram, plot/save the diagram and return the gene lists
vennPlot <- function(filteredRes1, filteredRes2, savePlot=TRUE, saveTitle="vennPlot", saveGeneLists=TRUE, listNames=c("A","B"), title="vennPlot", colors = c("red", "lightblue"), outDir = "plots/venn", verbose = TRUE) {
  
  geneLists <- list(filteredRes1$gene_id, filteredRes2$gene_id); 
  names(geneLists) <- listNames

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.plot <- venn.diagram(geneLists , NULL, fill=colors, alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=listNames, main=title)
  print(grid.draw(venn.plot))
  
  if (savePlot) {
    dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(outDir,"/",saveTitle,".pdf"))
    grid.draw(venn.plot)
    dev.off()
  }
  
  vennObject <- venn(geneLists, show.plot=FALSE)
  intersections <- attr(vennObject,"intersections")
  
  max.len <- max(length(unlist(intersections[1])), length(unlist(intersections[2])), length(unlist(intersections[3])))
  x <- c(unlist(as.vector(intersections[1]), use.names = FALSE), rep(NA, max.len - length(unlist(intersections[1]))))
  y <- c(unlist(as.vector(intersections[2]), use.names = FALSE), rep(NA, max.len - length(unlist(intersections[2]))))
  z <- c(unlist(as.vector(intersections[3]), use.names = FALSE), rep(NA, max.len - length(unlist(intersections[3]))))

  df <- data.frame(x = x, t = y, z = z)
  colnames(df) <- c(listNames, paste0(listNames[1], listNames[2]))
  if(saveGeneLists) {write.table(df, file=paste0(outDir,"/",saveTitle,".geneLists.txt"), quote = FALSE, sep="\t", col.names = NA)}
  
  go.obj <- newGeneOverlap(listA = filteredRes1$ensembl, listB = filteredRes2$ensembl, spec="hg19.gene")
  go.obj <- testGeneOverlap(go.obj)
  #print(paste0("Overlap Fisher's exact test. Overlap genes: ",length(go.obj@intersection),". p-value: ",go.obj@pval))
  
  if (verbose==TRUE) {print(go.obj)}
  
  return(df)

}