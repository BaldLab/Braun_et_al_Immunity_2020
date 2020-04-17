
# Function for single gene analysis of the TCGA
# will do DGE and OS survival plot as well as some other plots

# must set working directory to file location
DGE.and.OS.plot.function <- function(data.set.name, 
                                     data.loaded = FALSE, 
                                     rda.file.path, 
                                     goi = "ENSG00000150637", 
                                     gene.name = "DNAM1",
                                     genes.to.plot = NULL,
                                     perform.DEG.analysis = FALSE){
  
  # data.set.name = character of datasetID in TCGA database i.e "TCGA-GBM"
  # data.loaded = logical indicating if data has previously been downloaded - if true with prevent downloading data and instead just load the data from .rda file
  # rda.file.path = character pointing to the pathfile of the .rda file to load, must be set if data.load = T
  # genes.to.plot = provide named vector of gene IDs to plot expression and biplot against goi.
  
  
  require(rstudioapi)
  require(ggplot2)
  require(gplots)
  require(TCGAbiolinks)
  require(SummarizedExperiment)
  
  
  
  #Set wd to source file location
  primary.dir <- dirname(getActiveDocumentContext()$path)  
  
  setwd(primary.dir)
  #print(getwd())
  
  # Creating directories for output
  
  print("Creating directories")
  if(!dir.exists("Output")){dir.create("Output")}
  if(!dir.exists(paste0("Output/", data.set.name))){dir.create(paste0("Output/", data.set.name), recursive = T)}
  
  print("Changing working directory")
  
  working.dir <- paste0("Output/", data.set.name)
  setwd(working.dir)
  
  
  if(!dir.exists("Data")){dir.create("Data", recursive = T)}
  if(!dir.exists("Figures")){dir.create("Figures", recursive = T)}
  
  
  
  if(data.loaded == TRUE){
    
    print("loading dataset")
    
    setwd(primary.dir)                   
    
    load(file = rda.file.path)
    TCGA.raw.counts <- data
    rm(data)
    
    setwd(working.dir)
    
  }else{
    
    print("downloading dataset")
    
    
    query.raw.expression.hg38 <- GDCquery(project = data.set.name, 
                                          data.category = "Transcriptome Profiling",
                                          data.type = "Gene Expression Quantification", 
                                          experimental.strategy = "RNA-Seq",
                                          workflow.type = "HTSeq - Counts")
    
    
    
    GDCdownload(query.raw.expression.hg38, files.per.chunk = 6)
    
    TCGA.raw.counts <- GDCprepare(query = query.raw.expression.hg38,
                                  save = TRUE, 
                                  save.filename = paste0("Data/", data.set.name, "_Transcriptome_Profiling_raw_counts.rda"))
    
    
  }
  
  
  
  print("Normalising by GC content")
  
  dataNorm <- TCGAanalyze_Normalization(tabDF = TCGA.raw.counts, 
                                        geneInfo =  TCGAbiolinks::geneInfoHT,
                                        method = "gcContent")
  
  
  
  
  # filter genes to reduce low expressed genes and false positive results
  # quantile filter of genes
  print("Filtering dataset")
  
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile", 
                                    qnt.cut =  0.25)
  
  
  
  print("Writing normalised & Filtered data to file")
  
  # Writing Normalised data 
  write.table(dataFilt, paste0("Data/Normalised_filtered_", data.set.name, "_RNAseq_data.txt"), 
              sep = "\t", row.names = T, quote = F)
  
  
  
  # GOI subsetting
  goi.vector <- dataFilt[goi, ]
  
  
  Summary.value <- summary(goi.vector)
  lower.quartile <- ceiling(Summary.value[2])
  upper.quartile <- floor(Summary.value[5])
  
  High.exp.logic <- dataFilt[goi,] >= upper.quartile
  print(paste0("n with high expression of ", goi, " = ", sum(High.exp.logic)))
  
  Low.exp.logic <- dataFilt[goi,] <= lower.quartile
  print(paste0("n with low expression of ", goi, " = ", sum(Low.exp.logic)))
  
  
  samples.high <- dataFilt[, High.exp.logic]
  samples.low <- dataFilt[, Low.exp.logic]
  
  # Creating dataframe for heatmaps and export as this was used in Dif gene expression
  high.low.Df <- cbind(samples.high, samples.low)
  
  write.table(high.low.Df, paste0("Data/High_vs_low_", gene.name, "_samples_", data.set.name, "_RNAseq_data.txt"), 
              sep = "\t", row.names = T, quote = F)
  
  
  if(perform.DEG.analysis == TRUE){
    #Note: does batch correction need to be performed
    
    # Diff.expr.analysis (DEA) #should include year batch factor in analyssi 
    dataDEGs <- TCGAanalyze_DEA(mat1 = samples.low,
                                mat2 = samples.high,
                                Cond1type = paste0("Low", gene.name),
                                Cond2type = paste0("High", gene.name),
                                fdr.cut = 0.01 ,
                                pipeline = "edgeR",
                                logFC.cut = 1,
                                method = "glmLRT")
    
    
    print(paste0("number of DEGs = ", nrow(dataDEGs)))
    
    # write DEG table
    print("Writing table of differentially expressed genes")
    write.table(dataDEGs, paste0("Data/", data.set.name, "_DEG_HighvsLow_", gene.name, ".txt"), sep = "\t", quote = F)
    
    
  }
  
  
  
  # Plotting overall survival KM curve
  print("Downloading clinical data")
  
  
  clin.data <- GDCquery_clinic(data.set.name, "clinical")
  
  print("Identifying samples by Gene expression signature")
  # identifying samples by signature
  name.vec <- colnames(samples.high)
  x <- substr(name.vec, 1, 12)
  
  samples.high.sig.logic <- clin.data$submitter_id %in% x
  
  name.vec <- colnames(samples.low)
  x <- substr(name.vec, 1, 12)
  
  samples.low.sig.logic <- clin.data$submitter_id %in% x
  
  
  
  
  
  
  print("Denoting signature")
  # denoting signature
  clin.data$Signature[samples.high.sig.logic] <- paste0("High_", gene.name)
  
  clin.data$Signature[samples.low.sig.logic] <- paste0("Low_", gene.name)
  
  
  print("Generating OS plot")
  
  
  # Full plot with CI and risk table
  TCGAanalyze_survival(clin.data,
                       clusterCol = "Signature",
                       filename = paste0("Figures/", data.set.name, "survival.pdf"),
                       risk.table = TRUE,
                       color = c("Red", "Blue"),
                       main = paste0(data.set.name, gene.name, "High vs Low"), 
                       conf.int = TRUE,
                       height = 10, width=10)
  
  # Plot with Risk table but no CI
  TCGAanalyze_survival(clin.data,
                       clusterCol = "Signature",
                       filename = paste0("Figures/", data.set.name, "survival_noCI.pdf"),
                       risk.table = TRUE,
                       color = c("Red", "Blue"),
                       main = paste0(data.set.name, gene.name, "High vs Low"), 
                       conf.int = FALSE,
                       height = 10, width=10)
  
  # Plot with No risk table but with CI plotted
  TCGAanalyze_survival(clin.data,
                       clusterCol = "Signature",
                       filename = paste0("Figures/", data.set.name, "survival_noRisk.pdf"),
                       risk.table = FALSE,
                       color = c("Red", "Blue"),
                       main = paste0(data.set.name, gene.name, "High vs Low"), 
                       conf.int = TRUE,
                       height = 10, width=10)
  
  # Plot with No risk table or CI
  TCGAanalyze_survival(clin.data,
                       clusterCol = "Signature",
                       filename = paste0("Figures/", data.set.name, "survival_noCI_noRisk.pdf"),
                       risk.table = FALSE,
                       color = c("Red", "Blue"),
                       main = paste0(data.set.name, gene.name, "High vs Low"), 
                       conf.int = FALSE,
                       height = 10, width=10)
  
  
  
  
  
  
  
  if(perform.DEG.analysis == TRUE){
    ####################
    #### Heatmap of DEGs
    ####################
    print("Plotting heatmap of DEGs")
    
    DEG.logic <- rownames(high.low.Df) %in% rownames(dataDEGs)
    
    
    DEG.heatmap.df <- high.low.Df[DEG.logic, ]
    print(paste0("test: number of genes in heatmap df = ", nrow(DEG.heatmap.df)))
    print(paste0("test: number of genes in DEG df = ", nrow(dataDEGs)))
    
    
    pdf(paste0("Figures/Heatmap_DEGs_", gene.name, "_high_vs_low_", data.set.name, "_.pdf"))
    my.heatmap(DEG.heatmap.df, title = paste0(gene.name, "high vs low expression"), column.cluster = F)
    dev.off()
    
  }
  
  
  ########################################
  #### Sorted Gene Expression plots
  ########################################
  print("Plotting sorted gene expression")
  
  if(is.null(genes.to.plot)){
    genes.to.plot <- goi
    names(goi.vector) <- gene.name
  }
  
  
  for(i in 1:length(genes.to.plot)){
    tryCatch({
      plot.var <- dataFilt[genes.to.plot[i], ]
      
      pdf(paste0("Figures/", names(genes.to.plot)[i], "_expression_", data.set.name, "_.pdf"))
      my.plot(plot.var, title = names(genes.to.plot)[i])
      dev.off()
    }, error = function(e){cat("\n", "ERROR IN PLOTTTING ORDERED EXPRESSION:", conditionMessage(e), "\n")})
  }
  
  #########################
  #### Correlation biplots
  #########################
  print("Plotting correlation graphs")
  
  if(length(genes.to.plot) > 1){
    x.axis.var <- dataFilt[goi, ]
    
    for(i in 1:length(genes.to.plot)){
      
      if(goi != genes.to.plot[i]){
        tryCatch({
          y.axis.var <- dataFilt[genes.to.plot[i], ]
          
          pdf(paste0("Figures/", gene.name, "_vs_", names(genes.to.plot)[i], "_correlation_", data.set.name, "_.pdf"))
          Biplot.function(x.axis.var, y.axis.var, xtitle = gene.name, ytitle = names(genes.to.plot)[i])
          dev.off()
        }, error = function(e){cat("\n", "ERROR IN PLOTTTING biplot EXPRESSION:", conditionMessage(e), "\n")})
      }
    }
  }
  
  
  # move back to primary directory
  setwd(primary.dir)
  print(getwd())
  
  print("FINISHED analysis")
  
}


# additional functions for plots contained within TCGA analysis function
my.plot <- function(data.set, title){
  
  Var.sort <- sort(data.set, decreasing = T)
  
  Var.sort <- log2(Var.sort)
  
  median.value <- median(data.set)
  median.value <- log2(median.value)
  
  plot(Var.sort, type = "p", 
       col = "black",
       cex = 0.5,
       pch = 20,
       main = title,
       xlab = "Index",
       ylab = "Expression (log2)")
  
  abline(h = median.value, lwd = 1)
  
  mtext(paste0("n = ", length(data.set)), adj = 0)
  mtext(paste0("median = ", signif(median.value, 4)), adj = 1)
  
}

my.heatmap <- function(data, title.var = NULL, row.cluster = TRUE, column.cluster = TRUE, row.font.size = 0.4, col.font.size = 0.5){
  
  dissimfun <- function(x) {
    dist(x, method = "manhattan") # "euclidean", "maximum", "manhattan", "canberra", "binary"
  }
  
  clusterfun <- function(x) {
    hclust(x, method = "ward.D")
  }
  
  
  
  
  print(paste0("Sum of data pre.processing is.na = ", sum(is.na(data))))
  
  z <- scale(log2(data + 1), center = T, scale = T)
  zz <- scale(t(z), center = T, scale = T)
  
  
  print(paste0("Sum of data post.processing is.na = ",sum(is.na(zz)))) 
  
  
  
  heatmap.2(as.matrix(t(zz)),
            col = colorRampPalette(c("turquoise", "black", "red"))(100), #turquoise
            scale = c("column"),
            na.rm = TRUE, 
            trace = "none",
            Rowv = row.cluster,
            Colv = column.cluster, 
            distfun = dissimfun,
            hclustfun = clusterfun,
            dendrogram = "row",
            margins = c(10, 7),
            cexCol = col.font.size,
            cexRow = row.font.size,
            main = title.var,
            key = TRUE,
            keysize = 1.3, 
            na.color = "white",
            # breaks = break_scale, 
            # sepcolor = "black", 
            # rowsep = nrow(t(zz))/2, 
            srtCol = 90)
  
}

Biplot.function <- function(x, y, xtitle, ytitle){
  
  xlab.var <- xtitle
  ylab.var <- ytitle
  title.var <- paste0(xlab.var, " Vs. ", ylab.var)
  
  x <- log2(x)
  y <- log2(y)
  
  
  plot(x, y, type = "p", 
       col = "black",
       cex = 0.5,
       pch = 20,
       main = title.var,
       xlab = paste0(xlab.var, " (log2)"),
       ylab = paste0(ylab.var, " (log2)"))
  
  corr.value <- cor(x, y)
  corr.value <- signif(corr.value, 4)
  
  mtext(paste0("r = ", corr.value), adj = 0)
  
  
  y <- cor.test(x, y)
  corr.pval <- signif(y$p.value, 4)
  
  mtext(paste0("P = ", corr.pval), adj = 1)
  
}

# Example use of the function
# Function for single dataset plotting
# DGE.and.OS.plot.function("TCGA-DLBC", data.loaded = T, 
#                           rda.file.path = "Archive_Output/Data/DLBC/Raw_Counts/TCGA_DLBC_Transcriptome_Profiling_raw_counts.rda",
#                         genes.to.plot = goi.plot.vector)

