# parameters and functions to use in scritps for J15221


# hetamap colors
ramp <- colorRamp(c("#00BFC4", '#FFFFFF', "#F8766D"))
redBlue = rgb( ramp(seq(0, 1, length = 21)), max = 255)
# response colors
# labels for time points
timeLab  = c("Baseline", "C1D1", "Week 8")

# set colors
recistCol = c("Non-responder" = "red","Responder" = "blue", "NonResponder" = "red")
cbrCol = c("CR,PR,SD" = "blue","PD" = "red")
respCol = c("CR"="blue", "PR" = "#00BFC4", "SD" = 'yellow', "PD" = "red")
subtypeCol = c("HR+" = "red", "TNBC" = "blue")
biopsyCol = c(Breast = 'pink', Chest_Wall = 'grey', Liver = 'gold4', "Lung" = 'blue', Lymph_Node = 'cyan', Omentum = "brown")

timeCol = setNames(c("black", 'grey45', 'grey75'), c("Baseline","C1D1","Week_8"))

tnbcCol = c("BL1" = "yellow", "BL2" = "olivedrab3", "IM" = "burlywood4", "MSL" = "orange", "LAR" = "purple", "UNS" = "white", "NA" = "grey")

# colors for gene set heatmaps to differ from gene expression heatmaps
heatmapCol = colorRampPalette(rev(brewer.pal(n = 7, name =
                                               "RdYlBu")))(100)

##Create All Facet for plotting
CreateAllFacet <- function(df, col){
  df$facet <- df[[col]]
  temp <- df
  temp$facet <- "All"
  merged <-rbind(temp, df)
  
  # ensure the facet value is a factor
  merged[[col]] <- as.factor(merged[[col]])
  
  return(merged)
}
# paired and unpaired Wilcoxon test
# k is a column in dat
getWilcox = function(k, dat)
{
  mat <- data.frame(matrix(nrow = length(unique(dat$patientID)),
                           ncol = 3, dimnames = list(unique(dat$patientID),
                                                     levels(factor(dat$Time)))), check.names = F)
  for(i in levels(factor(dat$Time)))
  {
    d1 = dat %>% filter(Time == i)
    mat[d1$patientID, i] = d1[,k]
  }
  # add response
  m = dat %>% select(patientID,  RECIST) %>% filter(!(duplicated(patientID)))
  rownames(m) = m$patientID
  mat$RECIST = m[rownames(mat),"RECIST"]
  
  # running paired test
  # Baseline vs C1D1
  p1 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,2]), paired = T)$p.value
  #  C1D1 vs  C2D1
  p2 = wilcox.test(as.numeric(mat[,2]), as.numeric(mat[,3]), paired = T)$p.value
  #  Baseline vs  C2D1
  p3 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,3]), paired = T)$p.value
  
  # unpaired (by response)
  p4 = wilcox.test(as.numeric(mat[,1])~ mat[,"RECIST"])$p.value
  #  
  p5 = wilcox.test(as.numeric(mat[,2])~ mat[,"RECIST"])$p.value
  # 
  p6 = wilcox.test(as.numeric(mat[,3])~ mat[,"RECIST"])$p.value
  
  # mean of delta between time points
  d1 = mean(as.numeric(mat[,2]) -  as.numeric(mat[,1]), na.rm = T)
  d2 = mean(as.numeric(mat[,3]) -  as.numeric(mat[,2]), na.rm = T)
  d3 = mean(as.numeric(mat[,3]) -  as.numeric(mat[,1]), na.rm = T)
  
  
  return(c(pValue_BvsT1 = p1, 'mean(T1-B)' = d1, pValue_T1vsT2 = p2, 'mean(T2-T1)' = d2, pValue_BvsT2 = p3,'mean(T1-B)' = d3, 
           pValue_B_RvsNR = p4, pValue_T1_RvsNR = p5, pValue_T2_RvsNR = p6))
}

# paired and unpaired Wilcoxon test
# test difference between responders and non-responders
# k is a column in dat
getWilcoxDelta = function(k, dat)
{
  mat <- data.frame(matrix(nrow = length(unique(dat$patientID)),
                           ncol = 3, dimnames = list(unique(dat$patientID),
                                                     levels(factor(dat$Time)))), check.names = F)
  for(i in levels(factor(dat$Time)))
  {
    d1 = dat %>% filter(Time == i)
    mat[d1$patientID, i] = d1[,k]
  }
  # add response
  m = dat %>% select(patientID,  RECIST) %>% filter(!(duplicated(patientID)))
  rownames(m) = m$patientID
  mat$RECIST = m[rownames(mat),"RECIST"]
  # create treatment effect (calculate delta)
  # substract baseline from post run-in and week 8 
  mat$delta_T1_B =  mat[,2]-mat[,1]
  mat$delta_T2_B =  mat[,3]-mat[,1]
  # and week 8 from post run-in 
  mat$delta_T2_T1 =  mat[,3]-mat[,2]
  
  # running paired test
  # Baseline vs C1D1
  p1 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,2]), paired = T)$p.value
  #  C1D1 vs  C2D1
  p2 = wilcox.test(as.numeric(mat[,2]), as.numeric(mat[,3]), paired = T)$p.value
  #  Baseline vs  C2D1
  p3 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,3]), paired = T)$p.value
  
  # unpaired (by response)
  p4 = wilcox.test(as.numeric(mat[,1])~ mat[,"RECIST"])$p.value
  #  
  p5 = wilcox.test(as.numeric(mat[,2])~ mat[,"RECIST"])$p.value
  # 
  p6 = wilcox.test(as.numeric(mat[,3])~ mat[,"RECIST"])$p.value
  
  # p-value for delta between time points
  # compare differences in treatment effect in responders vs non-responders
  p7 = wilcox.test(delta_T1_B ~ RECIST, data = mat)$p.value
  p8 = wilcox.test(delta_T2_B ~ RECIST, data = mat)$p.value
  p9 = wilcox.test(delta_T2_T1 ~ RECIST, data = mat)$p.value
  
  # mean of delta between time points
  d1 = mean(mat$delta_T1_B, na.rm = T)
  d2 = mean(mat$delta_T2_B, na.rm = T)
  d3 = mean(mat$delta_T2_T1, na.rm = T)
  
  
  return(c(pValue_BvsT1 = p1, 'mean(T1-B)' = d1, # p-value and mean of difference between time points
           pValue_T1vsT2 = p2, 'mean(T2-T1)' = d2, 
           pValue_BvsT2 = p3,'mean(T1-B)' = d3, 
           # R vs NR for each time point
           pValue_B_RvsNR = p4, pValue_T1_RvsNR = p5, pValue_T2_RvsNR = p6,
           # difference between time points compared between R and NR
           pValue_d_T1_B_RvsNR = p7, pValue_d_T2_B_RvsNR = p8, 
           pValue_d_T2_T1_RvsNR = p9))
}

# splits a string with gene names and fixes names
# removeNeg is a flag to remove genes with negative sign
fixGeneNames = function(x, removeNeg = F)
{
  library(stringr)
  g = unlist(strsplit(x, ','))
  if (removeNeg) 
  {
    # find all genes with hyphen at the end
    match = grep('[-]$',g)
    # removes genes with hyphens at the end
    if (length(match) > 0) g = g[-match]
  }

  # remove hyphens in the end
  # find all genes with hyphen at the end
  match = grep('[-]$',g)
  # removes hyphens at the end
  g[match] = str_remove_all(g[match], '[-]$')
  # separates genes with dots
  g = unlist(strsplit(g, '.', fixed = T))
  # remove unnecessary words
  g = gsub("lo","",g, fixed = T)
  g = gsub("low/","",g, fixed = T)
  g = gsub("hi","",g, fixed = T)
  g = gsub("(migratory cDC2 only)","",g, fixed = T)
  # remove white spaces
  g = str_trim(g)
  # remove + at the end
  g = str_remove_all(g, '[+]$')
  # remove * at the end
  g = str_remove_all(g, '[*]')
  # convert to upper case
  str_to_upper(g)
  #  print(g)
}

# function to make barplots
plotBars = function(dat, pos = "stack", plot_avg = F, plot_sum = F)
{
  
  df = data.frame(dat, sampInfo[rownames(dat), c("RECIST","Time", "Type")],
                  patientID = sampInfo[rownames(dat),"Sample.Group.ID"], check.names = F)
  if (plot_sum){
    df_melt = melt(df, id.vars = c("RECIST","Time", "Type"), measure.vars = colnames(dat))
    # Stacked bar plot
    print(ggplot(df_melt, aes(fill=variable, y=as.numeric(value), x=Time)) + 
            geom_bar(position=pos, stat="identity") +
            ylab("Summarized Cell Proportions") +
            facet_wrap("RECIST"))
  }
  # plot average
  if (plot_avg){
    avg = aggregate(df[,colnames(dat)],by = list(df$Time, df$RECIST), FUN = mean, na.rm = T)
    colnames(avg)[1:2] = c("Time","RECIST")
    df_melt = melt(avg, id.vars = c("Time","RECIST"), measure.vars = colnames(dat))
    # Stacked bar plot
    print(ggplot(df_melt, aes(fill=variable, y=as.numeric(value), x=Time)) + 
            geom_bar(position=pos, stat="identity") +
            ylab("Average Cell Proportions") +
            facet_wrap("RECIST"))
    
  }
}

# functions for DE analysis

# do DEseq analysis with the specified design and return the specified res
getDEres = function(counts, annot, design, contrast, lfcShrink = T)
{
  # check if counts and annotation have the same samples
  if(identical(colnames(counts), rownames(annot))){
    s = colnames(counts)
  }else{
    s = intersect(colnames(counts), rownames(annot))
    cat("There is a different number of samples in counts and annottion.\n")
    print(s)
  }
  dds <- DESeqDataSetFromMatrix(countData = round(counts), 
                                colData = annot, 
                                design=as.formula(design))
  dds <- DESeq(dds)
  
  res <- results(dds, name = contrast)
  # apply lfcShrink
  if (lfcShrink)
    res <- lfcShrink(dds, coef=contrast, res=res)
  # remove NAs
  res <- res[!is.na(res$padj),]
  return(res)
  
}

# make barplot of fgsea results
fgseaBarplot = function(res, 
                        adjPvalThr = 0.05, 
                        main_title = "Pathways Enrichment Score",
                        sig_only = FALSE)
{
  ## Complete  Barplot
  res$adjPvalue <- ifelse(res$padj <= adjPvalThr, "significant", "non-significant")
  cols <- c("non-significant" = "grey", "significant" = "red")
  if(!sig_only)
  {
    p = ggplot(res, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      theme(axis.text.y = element_text(size=5)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score", title=main_title) + 
      theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
    print(p)
  }
  ## Sig Only  Barplot
  Sig_ORGResHALLMARK <- res[padj < adjPvalThr]
  if (nrow(Sig_ORGResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
  )} else {
    # prints the number of significant pathways
    cat("The number of significant pathways: ",nrow(Sig_ORGResHALLMARK))
    ggplot(Sig_ORGResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      theme(axis.text.y = element_text(size=5)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=paste0("Significant (pAdj<",adjPvalThr,")", main_title)) + 
      theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
  }
  
}
# print the number of pairs between time points
# the input table should have patientID and timepoint columns
printPairs = function(dat, patID = "patientID", timepoint = "timepoint")
{
  library(Matrix)
  times = levels(as.factor(dat[,timepoint]))
  # create a sparse matrix and populate with 1 
  mat = sparseMatrix(i = as.numeric(as.factor(dat[,patID])), 
                     j = as.numeric(as.factor(dat[,timepoint])), 
                     x = 1, 
                     dimnames = list(levels(as.factor(dat[,patID])),times))
  mat = as.data.frame(as.matrix(mat))
  
  cat("The number of paired samples between ",times[1]," and ", times [2],":",
      sum(rowSums(mat[,times[c(1,2)]])==2),"\n")

  cat("The number of paired samples between ",times[1]," and ", times [3],":",
      sum(rowSums(mat[,times[c(1,3)]])==2),"\n")
  
  cat("The number of paired samples between ",times[2]," and ", times [3],":",
      sum(rowSums(mat[,times[c(2,3)]])==2),"\n")
  
  return(mat)
  
}
# plot volcano plot and save top DE genes
plotVolcano = function(res, fdrThr = 0.05, logFCthr = 0.5, 
                       fileName = "topDE_genes.csv", main_title = "",
                       n = 10,  # add labels to top n genes to reduce overplotting
                       ...)
{
  TopDE <- res[res$padj < fdrThr &
                 abs(res$log2FoldChange) > logFCthr,]
  # order by absolute logFC
  TopDE = TopDE[order(abs(TopDE$log2FoldChange), decreasing = T),]
  # add labels to top N genes to reduce overplotting
  lab = row.names(TopDE)
#  print(lab)
  if (nrow(TopDE) > n) lab = row.names(TopDE)[1:n]
  p = EnhancedVolcano(res, lab = rownames(res), 
                      selectLab= lab,
                  x='log2FoldChange', y='padj', FCcutoff = logFCthr,
                  pCutoff = fdrThr, title = main_title, 
                  drawConnectors = T, maxoverlapsConnectors=50, 
                  labSize = 3,...)
  print(p)
  # print some statistics
  cat("The number of significantly DE genes:\n",
      nrow(TopDE),"\n")
  cat("The number of down- and up-regulated genes:\n",
      table(TopDE[,"log2FoldChange"] >0),"\n")
  TopDE = as.data.frame(TopDE)
  # list genes up and down regulated
  cat("The top up-regulated genes:\n",
      TopDE %>% filter(log2FoldChange > 0) %>% rownames(),"\n")
  cat("The top down-regulated genes:\n",
      TopDE %>% filter(log2FoldChange < 0) %>% rownames(),"\n")
  
  if (is.null(fileName)) return(rownames(TopDE))
  write.csv(res, fileName)
  return(rownames(TopDE))
}

# make heatmap
plotHA = function(anno, dat, topDE = rownames(dat), title = "", 
                  scale = T, ...)
{
  # create a data.frame with annotation to be used at the top of a heatmap
  dfAnnot = data.frame(Timepoint=anno[,'Timepoint'],
                       Type = anno[,'Subtype'],
                       Biopsy_site = anno[,'biopsy_site'],
                       RECIST=anno[,'RECIST'],
                       Response = factor(anno[,"Best_Resp"], levels = c("CR","PR","SD","PD")))
  # use that data.frame to create HeatmapAnnotation object to be added to the heatmap
  HA <- HeatmapAnnotation(df = dfAnnot, 
                          col = list(Timepoint = timeCol,
                                     Type = subtypeCol,
                                     Biopsy_site = biopsyCol,
                                     RECIST = recistCol,
                                     Response = respCol))
  # create heatmap with clustering by Pearson correlation for rows and columns
  d = dat[topDE,rownames(anno)]
  # scale data for the visualization
  if (scale) d = t(apply(d,1,scale))
  
  colnames(d) = rownames(anno)
  ht = Heatmap(d,top_annotation = HA,
          clustering_distance_rows = 'pearson', 
          clustering_distance_columns = 'pearson', 
          column_title = title, show_column_names = T,
          ,...)
 print(ht)  
}
# make heatmap
plotHA_pam50 = function(anno, dat, topDE = rownames(dat), title = "", 
                  scale = T, ...)
{
  # create a data.frame with annotation to be used at the top of a heatmap
  dfAnnot = data.frame(Timepoint=anno[,'Timepoint'],
                       Type = anno[,'Subtype'],
                       Biopsy_site = anno[,'biopsy_site'],
                       RECIST=anno[,'RECIST'],
                       Response = factor(anno[,"Best_Resp"], 
                                         levels = c("CR","PR","SD","PD")),
                       PAM50 = anno[,"assignment"])
  
  # use that data.frame to create HeatmapAnnotation object to be added to the heatmap
  HA <- HeatmapAnnotation(df = dfAnnot, 
                          col = list(Timepoint = timeCol,
                                     Type = subtypeCol,
                                     Biopsy_site = biopsyCol,
                                     RECIST = recistCol,
                                     Response = respCol,
                                     PAM50 = pam50::pam50_palette))
  # create heatmap with clustering by Pearson correlation for rows and columns
  d = dat[topDE,rownames(anno)]
  # scale data for the visualization
  if (scale) d = t(apply(d,1,scale))
  
  colnames(d) = rownames(anno)
  ht = Heatmap(d,top_annotation = HA,
          clustering_distance_rows = 'pearson', 
          clustering_distance_columns = 'pearson', 
          column_title = title, show_column_names = T,
          ,...)
  print(ht)
}

# make heatmap
plotHA_pam50_tnbc = function(anno, dat, topDE = rownames(dat), title = "", 
                        scale = T, ...)
{
  # create a data.frame with annotation to be used at the top of a heatmap
  dfAnnot = data.frame(Timepoint=anno[,'Timepoint'],
                       Type = anno[,'Subtype'],
                       Biopsy_site = anno[,'biopsy_site'],
                       RECIST=anno[,'RECIST'],
                       Response = factor(anno[,"Best_Resp"], 
                                         levels = c("CR","PR","SD","PD")),
                       PAM50 = anno[,"assignment"],
                       TNBCtype = anno[,"TNBCtype"])
  
  # use that data.frame to create HeatmapAnnotation object to be added to the heatmap
  HA <- HeatmapAnnotation(df = dfAnnot, 
                          col = list(Timepoint = timeCol,
                                     Type = subtypeCol,
                                     Biopsy_site = biopsyCol,
                                     RECIST = recistCol,
                                     Response = respCol,
                                     PAM50 = pam50::pam50_palette,
                                     TNBCtype = tnbcCol))
  # create heatmap with clustering by Pearson correlation for rows and columns
  d = dat[topDE,rownames(anno)]
  # scale data for the visualization
  if (scale) d = t(apply(d,1,scale))
  
  colnames(d) = rownames(anno)
  ht = Heatmap(d,top_annotation = HA,
               clustering_distance_rows = 'pearson', 
               clustering_distance_columns = 'pearson', 
               column_title = title, show_column_names = T,
               ,...)
  print(ht)
}


# make heatmap
plotHA_cbr = function(anno, dat, topDE = rownames(dat), 
                      title = "", scale = T,...)
{
  library(ComplexHeatmap)
  # create a data.frame with annotation to be used at the top of a heatmap
  dfAnnot = data.frame(Timepoint=anno[,'Timepoint'],
                       Type = anno[,'Subtype'],
                       Biopsy_site = anno[,'biopsy_site'],
                       RECIST=anno[,'RECIST'],
                       CBR = anno[,'CBR'],
                       Response = factor(anno[,"Best_Resp"], levels = c("CR","PR","SD","PD")))
  # use that data.frame to create HeatmapAnnotation object to be added to the heatmap
  HA <- HeatmapAnnotation(df = dfAnnot, col = list(Timepoint = timeCol, 
                                                   Type = subtypeCol,
                                                   Biopsy_site = biopsyCol, 
                                                   RECIST = recistCol,
                                                   CBR = recistCol,
                                                   Response = respCol))
  # create heatmap with clustering by Pearson correlation for rows and columns
d = dat[topDE,rownames(anno)]
# scale data for the visualization
if (scale) d = t(apply(d,1,scale))

colnames(d) = rownames(anno)
    Heatmap(d,top_annotation = HA,
          clustering_distance_rows = 'pearson', 
          clustering_distance_columns = 'pearson', 
          column_title = title, show_column_names = T,
          column_labels = colnames(d),...)
}

# creates boxplots for responders vs non-responders split by time points
# for CIBERSORT
boxplotSplitByTime_ciber = function(dat, i)
{
  df = data.frame(cell_proportion = as.numeric(dat[,i]), 
                  sampInfo[rownames(dat), c("RECIST","Time", "Type", "pubID", "biopsy_site")],
                  patientID = sampInfo[rownames(dat),"Sample.Group.ID"],
                  check.names = F)
  
  p= ggplot(df, aes(x=RECIST, y=cell_proportion)) + 
    geom_boxplot() + 
    ylab(i) + 
    stat_compare_means() + 
    geom_point() + 
    facet_wrap("Time")+ 
    theme_bw()
  
  print(p)
  return(df)
}
# for gene expression
boxplotSplitByTime_expr = function(dat, i, sampInfo, timepoint= "Timepoint")
{
  df = data.frame(expression = as.numeric(dat[i,]), 
                  sampInfo[colnames(dat), ],
                  check.names = F)
  
  p= ggplot(df, aes(x=RECIST, y=expression)) + 
    geom_boxplot() + 
    ylab(i) + 
    stat_compare_means() + 
    geom_point() + 
    facet_wrap(timepoint)+ 
    theme_bw()
  
  print(p)
}

# correlates expression of cell markers with cell proportions
# type1 is a cell type name in Edgar's list
# type2 is a cibersort cell type
plotCor = function(type1, type2)
{
  #check if the type is in markers list. If not, print the name
  if (type1 %in% rownames(markers))
  {
    # genes to correlate
    g = intersect(rownames(rnaData),fixGeneNames(markers[type1,2]) )
    if (length(g) == 0) return(NULL);
    corMat = matrix(nrow = length(g), ncol = 3, dimnames = list(g,timeLab))
    for(i in timeLab){
      # samples to correlate
      s = intersect(colnames(rnaData), sampInfo %>% filter(Time == i) %>% rownames())
      # combine expression and cibersort results
      dat = cbind(t(rnaData[g,s]), ciberRes[s, type2])
      # run correlation
      corRes = cor(dat, method = 'spearman')
      corMat[g,i] = corRes[g,length(g)+1]
    }
    # plot correlation coefficients cells vs genes
    p = pheatmap(corMat, scale = "none", main = paste0('Spearman correlation coefficients. ', type1), cluster_rows = F, cluster_cols = F, breaks = seq(-1, 1, length.out = 21), color = redBlue, legend_labels = "Spearman correlation")
    
    print(p)
    
  }else{print(type1)}
}

# function to fix sample names after TNBCtype classification
fix_sample_names = function(s){
  # trim white spaces
  s = trimws(s)
  s = gsub("R.","R-",s, fixed = T)
  s = gsub("Post.Run.in","Post Run-in",s, fixed = T)
  s = gsub("Week.8","Week 8",s, fixed = T)
  return(s)
}
# Custom sorting function using stringr
custom_sort <- function(x) {
  # Extract numeric part
  nums <- as.numeric(str_extract(x, "\\d+"))
  # Order based on numeric part
  order(nums)
}

