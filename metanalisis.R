library(mixOmics)
library(AnnotationDbi)
library(ArrayExpress)
library(Biobase)
library(BiocGenerics)
library(clusterProfiler)
library(dplyr)
library(factoextra)
library(futile.logger)
library(GEOquery)
library(ggplot2)
library(IRanges)
library(limma)
library(MGSDA)
library(org.Hs.eg.db)
library(plotly)
library(S4Vectors)
library(stringr)
library(VennDiagram)
library(metafor)
library(knitr)



datos <- load("./metanalysis/data/SN_datasets.RData")



# BOXPLOT FUNCTION
boxplot_ggplot <- function(data,
                           groups,
                           title = "Boxplot. Normalized expression",
                           bottom = FALSE,
                           save = NULL,
                           width = 800,
                           height = 600) {
  require(reshape2)
  
  box_data <- melt(data,
                   varnames = c("gene", "sample"),
                   as.is = TRUE)
  
  names(groups) <- colnames(data)
  box_data$Condition <- groups[(box_data$sample)]
  
  ## Boxplot
  g = ggplot(box_data, aes(x = sample, y = value, fill = Condition)) +
    geom_boxplot(outlier.size = 0.001) +
    xlab("Samples") +
    ylab("Expression") +
    ggtitle(title) +
    labs(fill = "Groups") +
    theme(
      axis.text.x = element_text(angle = 90, size = 0.5),
      axis.text.y = element_text(size = 15),
      title = element_text(size = 20),
      text = element_text(size = 15)
    )
  
  if (bottom == TRUE) {
    # Plot legend
    g =  g + theme(
      legend.position = "bottom",
      legend.key.size = unit(0.5, "line"),
      legend.title = element_text(size = 10)
    )
  }
  
  if (!is.null(save)) {
    png(filename = save,
        width = width,
        height = height)
    print(g)
    dev.off()
  }
  
  return (g)
  
}

# DENDOGRAM FUNCTION
treeclust_from_ggplot <- function(data,
                                  groups,
                                  title = "Clustering. Correlation distance",
                                  subtitle = NULL,
                                  bottom = FALSE,
                                  dist = "cor",
                                  save = NULL,
                                  width = 800,
                                  height = 600) {
  require(ggplot2)
  require(ggdendro)
  
  ## Create cluster
  if (dist == "cor") {
    correlacion <- cor(data)
    distancia <- as.dist((1 - correlacion) / 2)
  } else if (dist == "euclid") {
    distancia <- dist(t(data), method = "euclidean") # dist needs variables as columns
  } else{
    stop("Please specify an accepted distance. Options are 'cor' or 'euclid'")
  }
  
  cluster <- hclust(distancia)
  cluster$clase <- groups
  
  ##
  dendr <- ggdendro::dendro_data(cluster, type = "rectangle")
  ##
  clases <- as.character(cluster$clase)
  clust.df <- data.frame(label = cluster$labels, Condition = factor(clases))
  ##
  dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
  ##
  cbPalette <- c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )
  g = ggplot() +
    geom_segment(data = ggdendro::segment(dendr), aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    )) +
    geom_text(
      data = label(dendr),
      aes(
        x,
        y,
        label = label,
        hjust = 0,
        color = Condition
      ),
      size = 2
    ) +
    coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      # title = element_text(size = 20),
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_colour_hue() +
    scale_fill_hue()
  
  if (!is.null(subtitle)) {
    # Añadir subtítulo
    g <- g + labs(title = title, subtitle = subtitle)
  } else{
    g <- g + ggtitle(title) + theme()
  }
  
  if (bottom == TRUE) {
    #leyenda bajo el gráfico
    g =  g + theme(
      legend.position = "bottom",
      legend.key.size = unit(0.5, "line"),
      legend.title = element_text(size = 6)
    )
  }
  
  if (!is.null(save)) {
    png(filename = save,
        width = width,
        height = height)
    print(g)
    dev.off()
  }
  
  return (g)
}

# EXPLORATORY ANALYSIS FUNCTION
plotEDA <- function(study, l.studies) {
  # Create plots directory
  wd <- "./Plots/"
  dir.create(wd, recursive = TRUE)
  
  # Extract expression set from list
  array <- l.studies[[study]]
  
  # BARPLOT
  dplot <- pData(array)
  agg <- count(dplot, sex, brain_region, disease, .drop = FALSE)
  agg_ord <- mutate(agg,
                    disease = reorder(disease, -n, sum),
                    sex = reorder(sex, -n, sum))
  ggplot(agg_ord, aes(x = disease, y = n, fill = sex)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = study) +
    geom_text(
      aes(label = n),
      vjust = 1.6,
      position = position_dodge(0.9),
      size = 4
    ) +
    theme_minimal() +
    facet_wrap(~ brain_region) +
    theme(strip.text.x = element_text(size = 12),
          plot.background = element_rect(fill = "white"))
  
  ggsave(paste0(wd, study, "_bar.png"),
         width = 12,
         height = 5)
  
  # Create groups from phenotypic variables
  dplot$groups <- with(dplot, paste(sex, brain_region, disease, sep = "_"))
  
  # BOXPLOT
  boxplot_ggplot(exprs(array), dplot$groups, title = paste0(study, " - Boxplot"))
  ggsave(paste0(wd, study, "_box.png"),
         width = 12,
         height = 8)
  
  # DENDOGRAM
  treeclust_from_ggplot(exprs(array), dplot$groups)
  ggsave(paste0(wd, study, "_dendo.png"),
         width = 12,
         height = 10)
  
  # PCA
  res.pca <- prcomp(t(exprs(array)))
  fviz_pca_ind(
    res.pca,
    habillage = factor(dplot$groups),
    geom = "point",
    pointsize = 4,
    invisible = "quali"
  ) +
    theme(plot.background = element_rect(fill = "white"))
  
  ggsave(paste0(wd, study, "_PCA.png"),
         width = 12,
         height = 8)
  
  
}

diffexprs <- function(X, batch = F) {
  sexdisease = vector("list", ncol(X))
  for (i in seq_along(sexdisease)) {
    sexdisease[[i]] = paste(pData(X)[, "sex"][i], pData(X)[, "disease"][i], sep =
                              "")
  }
  
  sexdisease = factor(unlist(sexdisease))
  
  # Consider batch effect or not(default)
  if (batch == T) {
    batch = factor(pData(X)[, "batch"])
    design <- model.matrix(~ 0 + sexdisease + batch)
  } else {
    design <- model.matrix(~ 0 + sexdisease)
  }
  
  colnames(design)[1:4] <- levels(sexdisease)
  fit <- lmFit(X, design)
  cont.matrix <- makeContrasts(
    dif1 = femalecase - femalecontrol,
    dif2 = malecase - malecontrol,
    dif12 = (femalecase - femalecontrol) - (malecase - malecontrol),
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  return(fit2)
}


make_results <- function(arraysP, dir = wd) {
  if (!dir.exists(wd)) {
    dir.create(wd, recursive = TRUE)
  }
  
  
  # Empty lists of results
  ALLfit <- list()
  FCfit <- list()
  STfit <- list()
  SNfit <- list()
  
  
  for (study in names(arraysP)) {
    array <- arraysP[[study]]
    # All tissues DEA
    fit <- diffexprs(array)
    ALLfit[[study]] <- fit
    # Write significant
    for (i in c(1:3)) {
      cont <- topTable(
        fit,
        coef = i,
        adjust = "BH",
        sort = "none",
        n = Inf
      )
      cont <- cont[cont$adj.P.Val < 0.05, ]
      write.csv(cont, paste0(wd, study, "_ALL_cont", i, ".csv"))
    }
    # FC differential expression
    if ("FC" %in% pData(array)$brain_region) {
      filtered <- array[, pData(array)$brain_region == "FC"]
      fit <- diffexprs(filtered)
      FCfit[[study]] <- fit
      for (i in c(1:3)) {
        cont <- topTable(
          fit,
          coef = i,
          adjust = "BH",
          sort = "none",
          n = Inf
        )
        cont <- cont[cont$adj.P.Val < 0.05, ]
        write.csv(cont, paste0(wd, study, "_FC_cont", i, ".csv"))
      }
    }
    # ST differential expression
    if ("ST" %in% pData(array)$brain_region) {
      filtered <- array[, pData(array)$brain_region == "ST"]
      fit <- diffexprs(filtered)
      STfit[[study]] <- fit
      for (i in c(1:3)) {
        cont <- topTable(
          fit,
          coef = i,
          adjust = "BH",
          sort = "none",
          n = Inf
        )
        cont <- cont[cont$adj.P.Val < 0.05, ]
        write.csv(cont, paste0(wd, study, "_ST_cont", i, ".csv"))
      }
    }
    # SN differential expression
    if ("SN" %in% pData(array)$brain_region) {
      filtered <- array[, pData(array)$brain_region == "SN"]
      fit <- diffexprs(filtered)
      SNfit[[study]] <- fit
      for (i in c(1:3)) {
        cont <- topTable(
          fit,
          coef = i,
          adjust = "BH",
          sort = "none",
          n = Inf
        )
        cont <- cont[cont$adj.P.Val < 0.05, ]
        write.csv(cont, paste0(wd, study, "_SN_cont", i, ".csv"))
      }
    }
  }
  
  if (!dir.exists("./Output/")) {
    dir.create("./Output")
  }
  # Save lists of "fit" objects
  save(ALLfit, file = "./Output/ALLfit.RData")
  save(FCfit, file = "./Output/FCfit.RData")
  save(STfit, file = "./Output/STfit.RData")
  save(SNfit, file = "./Output/SNfit.RData")
  
  return(list(
    ALL = ALLfit,
    FCfit = FCfit,
    STfit = STfit,
    SNfit = SNfit
  ))
}

SE_array <- function(fit, contrast) {
  # STEP 0. Pre-processing previous data
  # ===============================================================
  ## Calculate SE
  #OPTION1: https://support.bioconductor.org/p/70175/ by Gordon Smyth/January Weiner
  #The effect sizes are contained in fit$coefficients
  summary(fit$coefficients)
  head(fit$coefficients)
  contrasts_names <- c("femalesHC", "malesHC", "diff_FemalesHC_MalesHC")
  contrast_name <- contrasts_names[contrast]
  #The standard errors can be obtained from 2 sub-options:
  # SE.coef <- sqrt(fit$s2.post) #JANUARY  (Here I have a problem when having several contrasts
  #                                        #because I have the same information for all contrasts)
  # head(SE.coef)
  # summary(SE.coef)
  SE.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled[, contrast] #GORDON
  #OPTION2: https://support.bioconductor.org/p/70175/ by  Steve Lianoglou  (SE PARECE A GORDON)
  allgenes <- topTable(
    fit,
    number = "all",
    confint = TRUE,
    adjust.method = "fdr",
    coef = contrast
  )
  allgenes[, "SE"] <- (allgenes[, "CI.R"] - allgenes[, "CI.L"]) / 3.92
  
  #final results
  table(rownames(SE.coef) == rownames(fit$coefficients))
  mat <- cbind(fit$coefficients[, contrast], SE.coef)
  colnames(mat) <- c("coef", "se.coef")
  head(mat)
  
  int <- intersect(rownames(allgenes), rownames(mat))
  length(int)
  res <- cbind(allgenes, mat[rownames(allgenes), ])
  head(res)
  dim(res)
  return(res)
}

# START
load("./metanalysis/data/SN_datasets.RData")
arraysP <- data
for (s in names(arraysP)) {
  plotEDA(s, arraysP)
}
metanalisis_fun <- function(SNfit, cont) {
  # Use the loaded list as function input
  EDs_sel <- lapply(SNfit, function(X)
    SE_array(fit = X, contrast = cont))
  
  contrasts_names <- c("femalesHC", "malesHC", "diff_FemalesHC_MalesHC")
  metaanalysis_name <- contrasts_names[cont]
  ## directories
  wd <- "./Results/MA/" #Output dir
  wd <- paste0(wd, metaanalysis_name, "/")
  dir.create(wd, recursive = TRUE)
  # STEP 1. Preparing input for meta-analysis: LOR and SE matrix
  # ===============================================================
  
  # we search a list including all unique ID genes for all studies
  genes <- NULL
  for (fi in EDs_sel) {
    genes <- c(genes, rownames(fi))
  }
  
  length (genes)
  genes <- unique (genes)
  length (genes)
  genes[grepl("[0-9]+-[A-Z][a-z]{2}", genes)]
  genes <- genes[!grepl("[0-9]+-[A-Z][a-z]{2}", genes)]
  genes[grepl("[0-9]+-[A-Z][a-z]{2}", genes)]
  genes <- base::sort(genes)
  
  
  ### generating matrix with all logFC for all studies
  mat.logFC <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel))
  rownames (mat.logFC) <- genes
  colnames (mat.logFC) <- names(EDs_sel)
  head (mat.logFC)
  
  for (i in 1:length(EDs_sel)) {
    co <- names(EDs_sel[i])
    res <- EDs_sel[[i]]
    logFC <- res$logFC
    names (logFC) <- (rownames (res))
    mat.logFC[, co] <- logFC[rownames(mat.logFC)]
  }
  ##### prueba
  #mat.logFC <- data.frame(Symbol=genes,row.names =genes)
  #for(co in names(EDs_sel1)){
  #  res <- EDs_sel1[[co]]
  #  res <- res[,c("Symbol","logFC")]
  #  colnames(res) <- c("Symbol",co)
  #  mat.logFC <- merge(mat.logFC,res,by = "Symbol",all.x = TRUE)
  #}
  #rownames(mat.logFC) <- genes
  #mat.logFC <- mat.logFC[,-1]
  head (mat.logFC)
  tail(mat.logFC)
  table (is.na(mat.logFC))
  dim (mat.logFC)
  
  # select genes included at least in 2 or more studies
  mat.logFC.NA <- is.na(mat.logFC)
  head(mat.logFC.NA)
  sum.NA <-  apply(mat.logFC.NA, 1, sum)
  table(sum.NA)
  min.sum.NA <- sum.NA < ncol(mat.logFC) - 1
  table(min.sum.NA)
  
  # filter by min.sum.NA
  mat.logFC <- mat.logFC[min.sum.NA == T, ]
  dim(mat.logFC)
  
  
  ### generating matrix with all SE for all studies
  mat.SE <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel))
  rownames (mat.SE) <- genes
  colnames (mat.SE) <- names(EDs_sel)
  head (mat.SE)
  
  
  # (SE FROM GORDON: se.coef)
  for (i in 1:length(EDs_sel)) {
    co <- gsub("_ED", "", names(EDs_sel[i]))
    res <- EDs_sel[[i]]
    SE <- res$se.coef
    names (SE) <- (rownames (res))
    mat.SE[, co] <- SE[rownames(mat.SE)]
  }
  
  
  head (mat.SE)
  tail(mat.SE)
  table (is.na(mat.SE))
  dim (mat.SE)
  
  # filter by min.sum.NA
  mat.SE <- mat.SE[min.sum.NA == T, ]
  dim(mat.SE)
  
  
  
  
  # STEP 2. Meta-analysis for genes
  # ===============================================================
  
  # suppose between-study variance is non-zero.
  # there are different methods to estimate this variance:
  # DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
  # Now we have logFC and SE  (not VARIANCE), so:
  # yi -> logFC   sei -> SE
  # result.lor <- rma(yi = mat.logFC[1, ],
  #                   sei = mat.SE[1, ],   #pay attention, not vi (varianze)
  #                   method = "DL") # DerSimonian-Laird.
  
  
  # explore the function to do the meta-analysis
  #?rma
  
  MA <- lapply(1:length(rownames(mat.logFC)), function(x) {
    rma(yi = mat.logFC[x, ],
        sei = mat.SE[x, ],
        method = "DL")
  })
  
  # MA <- lapply(1:length(rownames(mat.logFC)),
  #              function(x){rma(yi = mat.logFC[x, ],
  #                              sei = mat.SE[x, ],
  #                              method = "FE")})
  
  names (MA) <- rownames(mat.logFC)
  class (MA)
  length(MA)
  head (MA)
  MA[[1]]
  
  #result.logFC$pval      #p-value about logFC = 0
  #result.logFC$ci.lb     #IC down
  #result.logFC$ci.ub     #IC up
  #re sult.logFC$b         #estimation of combined logFC
  
  #data.frame including all detailed results:
  result_meta <- as.data.frame(do.call("rbind", lapply(MA, function(x) {
    c(x$ci.lb,
      x$b,
      x$ci.ub,
      x$pval,
      x$QE,
      x$QEp,
      x$se,
      x$tau2,
      x$I2,
      x$H2)
  })))
  
  colnames(result_meta) <- c(
    "lower_bound",
    "logFC",
    "upper_bound",
    "pvalue",
    "QE",
    "QEp",
    "SE",
    "tau2",
    "I2",
    "H2"
  )
  
  p.adjust.fdr <- stats::p.adjust(result_meta[, 4], method = "fdr")
  p.adjust.BY  <- stats::p.adjust(result_meta[, 4], method = "BY")
  result_meta <- round(cbind(result_meta, p.adjust.fdr, p.adjust.BY), 3)
  head(result_meta)
  
  
  # significant genes
  corte = 0.05
  table(result_meta[, "pvalue"] < corte)
  table(result_meta[, "p.adjust.fdr"] < corte)
  table(result_meta[, "p.adjust.BY"] < corte)
  
  
  # add number of studies where the gene is evaluated
  sum.NA <- sum.NA[sum.NA < ncol(mat.logFC) - 1]
  n.studies <-  ncol(mat.logFC) - sum.NA
  table(n.studies)
  n.studies <- n.studies [rownames(mat.logFC)]
  length(n.studies)
  result_meta[, "n.studies"]  <- n.studies
  head(result_meta)
  summary(result_meta$p.adjust.fdr)
  
  # no tenemos genes significativos en estos estudios y he subido el FDR para seleccioar
  # al menos un grupo con los que ver el funcionamiento del script
  
  #corte = 0.2
  sig.genes.df = result_meta[result_meta$p.adjust.fdr < corte, ]
  dim(sig.genes.df)
  
  write.table(
    x = sig.genes.df[order(sig.genes.df$p.adjust.fdr), ],
    file = paste0(wd, "sig.genes.tsv"),
    sep = "\t",
    quote = FALSE
  )
  write.table(
    x = result_meta[order(result_meta$p.adjust.fdr), ],
    file = paste0(wd, "all.genes.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
  all.genes <- result_meta[order(result_meta$p.adjust.fdr), ]
  # STEP 3. INFLUENCE AND SENSITIVITY ANALYSIS
  # ===============================================================
  
  #add 4 new variables about influence & sensitivity analysis:
  
  for (i in rownames(sig.genes.df)) {
    print(i)
    #define studies for each function (not NA)
    estudios <- colnames(mat.logFC)[!mat.logFC.NA[i, ]]
    
    #influence info 1:
    #number of studies where the sign of the logOR is the same  of the global logOR:
    sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi) == rep(sign(MA[[i]]$b), length(estudios)))
    
    #influence info 2: how many studies could be influencers?
    inf <- influence(MA[[i]])
    res <- paste(estudios[inf$is.infl], collapse = ",")
    sig.genes.df[i, "infl.nstudies"] <- ifelse(res == "", "non", res)
    
    #sensivity analysis
    l1 <- as.data.frame(leave1out(MA[[i]]))
    rownames(l1) <- estudios
    
    #1. p.value about differences between all estimates from leave one out
    #   and global estimate)
    sig.genes.df[i, "sensi.global"] <- t.test(x = l1$estimate, mu = as.numeric(MA[[i]]$b))$p.value
    #2. number of  studies where pvalue > 0.05
    # (we hope p-values < 0.05, significant estimates)
    res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")
    sig.genes.df[i, "sensi.specific"] <- ifelse(res2 == "", "all.p.values < 0.05", res2)
  }
  
  
  ## QUESTIONS TO ASSESS META-ANALYSIS FOR EACH FUNCTION:
  
  #1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
  table(sig.genes.df$infl.same.sign.logFC)
  
  #2. INFLUENCE STUDIES. How many functions including influence studies?
  table(sig.genes.df$infl.nstudies == "non")
  
  #3. SENSITIVITY. In global, are there many functions with differences in the estimate?
  table(sig.genes.df$sensi.global < 0.05)
  
  #4. SENSITIVITY.  How many functions including changes in the significance about
  # its new estimate  after leave1out?
  table(sig.genes.df$sensi.specific == "all.p.values < 0.05")
  
  
  #save final results:
  cat ("ID\t", file = paste0(wd, "sig.genes.df.txt"))
  write.table(
    sig.genes.df,
    file = paste0(wd, "sig.genes.df.txt"),
    sep = "\t",
    quote = F,
    append = TRUE,
    row.names = T
  )
  
  
  
  
  # STEP 4. Visualization of significant genes
  # ===============================================================
  
  #select significant functions to visualize:
  sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.05, ]
  #sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.2,]
  sig.results
  dim(sig.results)
  
  dir.create(paste0(wd, "plots"), recursive = TRUE)
  #setwd(paste0(wd,"plots"))
  
  selMethod <- "DL"
  
  for (i in 1:nrow(sig.results)) {
    mygenes <- rownames(sig.results)[i]
    res <- rma(yi = mat.logFC[mygenes, ],
               sei = mat.SE[mygenes, ],
               method = "DL")
    
    #FOREST PLOT
    # png (filename = paste("FOREST_", mygenes,".png", sep =""), width = 960 ,
    #      height = 960, res = 200)
    png (
      filename = paste0(wd, "plots/", gsub("-", "_", mygenes), "_FOREST", ".png"),
      width = 960 ,
      height = 960,
      res = 200
    )
    forest(
      res,
      slab = toupper(colnames(mat.logFC)),
      #Nombre de los estudios
      xlab = "logFC",
      cex = 0.7,
      mlab = paste(selMethod, "Model for All Studies", sep = " "),
      border = "black",
      #Color del borde del rombo
      col = "red",
      #Color del rombo
      main = paste("\n", mygenes, sep = "")
    )
    text(9, -3, "logFC [IC 95%]", pos = 2, cex = 0.7)
    dev.off()
    
    #FUNNEL PLOT
    png (
      filename = paste0(wd, "plots/", gsub("-", "_", mygenes), "_FUNNEL", ".png"),
      width = 960 ,
      height = 960,
      res = 200
    )
    par(mfrow = c(2, 2))
    funnel(
      res,
      main = "Standard Error",
      back = "darkslategray1",
      xlab = paste("logFC (", mygenes, ")", sep = "")
    )
    funnel(
      res,
      yaxis = "vi",
      main = "Sampling Variance",
      back = "darkslategray1",
      xlab = paste("logFC (", mygenes, ")", sep = "")
    )
    funnel(
      res,
      yaxis = "seinv",
      main = "Inverse Standard Error",
      back = "darkslategray1",
      xlab = paste("logFC (", mygenes, ")", sep = "")
    )
    funnel(
      res,
      yaxis = "vinv",
      main = "Inverse Sampling Variance",
      back = "darkslategray1",
      xlab = paste("logFC (", mygenes, ")", sep = "")
    )
    par(mfrow = c(1, 1))
    dev.off()
    
    #INFLUENCE PLOTS
    # That shows various diagnostic measures
    png (
      filename = paste0(wd, "plots/", gsub("-", "_", mygenes), "_INFLUENCE", ".png"),
      width = 960 ,
      height = 960,
      res = 200
    ) ##CAMBIAR
    inf <- influence(res)
    #plot(inf, plotfb = T)#"plotfb" is not a graphical parameter
    plot(inf)
    dev.off()
    
  }
  
  # STEP 5. Generating report
  # ===============================================================
  sig.genes.df <- sig.genes.df[order(sig.genes.df$p.adjust.fdr), ]
  save(sig.genes.df,
       result_meta,
       file = paste0(wd, metaanalysis_name, ".RData"))
  
  # Function to create multiple tabs
  make.tabs <- function(sig.genes.df) {
    res <- NULL
    for (g in rownames(sig.genes.df)) {
      file_name <- gsub("-", "_", g)
      res <- c(
        res,
        '### ',
        g,
        '{-} \n',
        "**Statistics of ",
        g,
        " meta-analisys** \n",
        "```{r, fig.align='center'}",
        '\n',
        "kable(sig.genes.df['",
        g,
        "',1:11])",
        '\n',
        '```',
        '\n',
        "[Gene information](https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        g,
        ") \n\n",
        "**Forest plot** \n",
        "```{r, fig.align='center'}",
        '\n',
        'knitr::include_graphics("',
        'plots/',
        file_name,
        '_FOREST.png")\n',
        '```',
        '\n',
        "**Funnel plot** \n",
        "```{r, fig.align='center'}",
        '\n',
        'knitr::include_graphics("',
        'plots/',
        file_name,
        '_FUNNEL.png")\n',
        '```',
        '\n',
        "**Incluence plot** \n",
        "```{r, fig.align='center'}",
        '\n',
        'knitr::include_graphics("',
        'plots/',
        file_name,
        '_INFLUENCE.png")\n',
        '```',
        '\n\n'
      )
    }
    return(res)
  }
  
  
  # Create the Rmd to knit
  ns <- nrow(sig.genes.df)
  cat(
    '---
title: "Meta-analysis of genes [DRAFT]"
output:
  html_document:
    toc: false
    toc_float: false
    code_folding: hide
    number_sections: true
    theme: spacelab
---
## ',
    metaanalysis_name,
    ' {.tabset .tabset-pills -}

```{r, warning=F, message=F}
library(dplyr)
library(knitr)
load("',
    metaanalysis_name,
    '.RData")
```  \n
### Significant results {-}  \n',
    "```{r, fig.align='center'}",
    '\n',
    "kable(sig.genes.df[,1:11], caption='Statistics of ",
    ns,
    " significant genes')",
    '\n',
    '```',
    '\n\n',
    make.tabs(sig.genes.df),
    "\n\n",
    '### sessionInfo {-}  \n',
    "```{r, fig.align='center'}",
    '\n',
    "date()",
    "\n",
    "sessionInfo()",
    '\n',
    '```',
    '\n\n',
    sep = "",
    file = paste0(wd, metaanalysis_name, ".Rmd")
  )
  dir.create(wd, recursive = TRUE)
  
  #setwd(wd)
  rmarkdown::render(paste0(wd, metaanalysis_name, ".Rmd"))
  ret <- list(all = all.genes, sig = sig.genes.df)
  return(ret)
  
}

wd <- "./Results/DEA/"

results_uni_block <- make_results(arraysP, dir = wd)

SNfit <- results_uni_block$SNfit

conts <- 1:3
res_final <-  lapply(conts, function(x)
  metanalisis_fun(SNfit, cont = x))

######################################
### 05_meta-analysis.R
### Author: Francisco García-García
### Modified by: Adolfo López-Cerdán
### Contact: fgarcia@cipf.es
###          adlpecer@gmail.com
######################################




## load data
