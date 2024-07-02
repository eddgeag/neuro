library(dplyr)
library(mixOmics)
library(limma)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(caret)
library(ggVennDiagram)
library(ggpubr)
library(cowplot)
library(ggstatsplot)
library(faraway)
library(ggExtra)
library(psych)
##=====Build Functions for the analysis===========
fun_load <- function(path_data, path_results_metanalysis) {
  ## load data
  #
  # path_data <- "./metanalysis/data/SN_datasets.RData"
  # path_results_metanalysis <- "./Results/DEA"
  load(file = path_data)
  
  datos <- data
  metanalisis_results_files <- list.files(path_results_metanalysis,
                                          pattern = ".csv",
                                          full.names = T)
  metanalisis_results <- lapply(metanalisis_results_files, read.csv)
  
  ## lets identidy the features from the metanalysis
  
  features_meta <- unique(unlist(lapply(metanalisis_results, function(x)
    x$Symbol)))
  
  ret <- list(
    metanalisis_results = metanalisis_results,
    features_meta = features_meta,
    datps = datos
  )
  return(ret)
  
}

prepare_data <- function(datos, metanalisis_results) {
  pheno <- lapply(datos, function(X)
    as.data.frame(X@phenoData@data))
  
  pheno_df <- as.data.frame(Reduce("rbind", pheno))
  
  study_number <- vector(mode = "character", length = dim(pheno_df)[1])
  
  n <- length(study_number) #total n
  
  s  <- length(pheno) # number of studies
  study_name <- paste0("study", 1:s)
  factor_study <- vector(mode = "list", length = s)
  names(factor_study) <- study_name
  #i in n observations
  for (s_m in 1:s) {
    r <- dim(pheno[[s_m]])[1] ## number of study
    study_m <- rep(study_name[s_m], r)
    factor_study[[s_m]] <- study_m
    
  }
  
  pheno_df$study <- as.factor(unlist(factor_study))
  
  ## obtain the data
  
  all_datasets <- vector(mode = "list", length = s)
  
  for (s_m in 1:s) {
    study_m <- datos[[s_m]]
    data_sm <- study_m@assayData$exprs
    all_datasets[[s_m]] <- data_sm
  }
  
  
  ## features on each data set
  
  all_features <- lapply(all_datasets, rownames)
  
  
  ## features in common
  
  all_features <- Reduce(intersect, all_features)
  
  ## for each data set lets keep the common features
  
  all_datasets_common_f <- lapply(all_datasets, function(X)
    X[rownames(X) %in% all_features, ])
  
  ## concatenate dataframes
  
  X <- as.data.frame(t(Reduce("cbind", all_datasets_common_f)))
  
  ## response categorical variables
  
  ## how many variables there are and the levels
  
  pheno_df_ <- pheno_df[, -ncol(pheno_df)]
  l <- length(colnames(pheno_df_))
  
  ## any unique value ?
  
  unique_values <- lapply(pheno_df_, unique)
  value <- lapply(unique_values, length)
  w <- which(value > 1)
  pheno_df_ <- pheno_df_[, w]
  pheno_df_ <- lapply(pheno_df_, as.factor)
  aux <- as.factor(unlist(Reduce("interaction", pheno_df_)))
  pheno_df_$interaction <- gsub("[.]", "_", aux)
  pheno_df_ <- as.data.frame(pheno_df_)
  interaccion <- pheno_df_$interaction
  design <- as.data.frame(model.matrix(~ 0 + interaccion))
  colnames(design) <- gsub("interaccion", "", colnames(design))
  design <- as.data.frame(design)
  ## make contrasts
  cont.matrix <- data.frame(
    dif1 = design$female_case - design$female_control,
    dif2 = design$male_case - design$male_control,
    dif12 = (design$female_case - design$female_control) -
      (design$male_case - design$male_control)
  )
  Y_list = list(cont1 = cont.matrix$dif1,
                cont2 = cont.matrix$dif2,
                cont3 = cont.matrix$dif12)
  ret <- list(
    X = all_datasets_common_f,
    Y = Y_list ,
    study = pheno_df$study,
    pheno_df = pheno_df
  )
  
  return(ret)
  
  
  
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
    design <- model.matrix( ~ 0 + sexdisease + batch)
  } else {
    design <- model.matrix( ~ 0 + sexdisease)
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

fun_extract_Accuracy <- function(models, data_, i) {
  predictions <- predict(models[[i]],
                         newdata = data_[[i]]$test$X,
                         study.test = data_[[i]]$test$study)
  
  clases_pc <- as.data.frame(predictions$MajorityVote$mahalanobis.dist)
  clases_pc <- as.data.frame(lapply(clases_pc, as.factor))
  Ytest <- as.factor(data_splitted_contrast$female_contrast$test$Y)
  
  accuracy <- mean(Reduce("c", lapply(clases_pc, function(x)
    confusionMatrix(data = x, Ytest))))
  
  return(accuracy)
  
  
}
fun_split <- function(X, Y, study, p) {
  set.seed(123456)
  n <- nrow(X)
  idx <- sample(1:n, round(n * p))
  X_train <- X[idx, ]
  X_test <- X[-idx, ]
  Y_train <- Y[idx]
  Y_test <- Y[-idx]
  study_train <- study[idx]
  study_test <- study[-idx]
  
  ret <- list(
    train = list(X = X_train, Y = Y_train, study = study_train),
    test = list(X = X_test, Y = Y_test, study = study_test)
  )
  return(ret)
  
  
}
fun_auc <- function(model, data, filename) {
  y_test <- as.factor(data$test$Y)
  p <- auroc(
    model,
    newdata = data$test$X,
    outcome.test = y_test,
    study.test = data$test$study
  )
  
  ggsave(filename = filename, plot = p[[2]])
  return(p)
  
}


##=== Load data=====


path_data <- "./metanalysis/data/SN_datasets.RData"
path_results_metanalysis <- "./Results/DEA/"


pre_data <- fun_load(path_data = path_data, path_results_metanalysis = path_results_metanalysis)
data_ready <- prepare_data(pre_data$datps, pre_data$metanalisis_results)
input_parafac <- data_ready$X

x1 <- data_ready$X[[1]]
x2 <- data_ready$X[[2]]
x3 <- data_ready$X[[3]]
x4 <- data_ready$X[[4]]
data_list <- list(x1, x2, x3, x4)

X <- t(as.matrix(Reduce("cbind", data_list)))
Y <- data_ready$study <- data_ready$study
set.seed(123456)
study <- data_ready$study
data_splitted_contrast <- lapply(1:3, function(y)
  fun_split(X, data_ready$Y[[y]], study = study, p = 0.7))

##====Load models and Search for the models performance ####
final_models <- readRDS("./final_models_1.rds")
names(data_splitted_contrast) <- names(final_models) <- c("female_contrast", "males_contrast", "diff_contrast")
### get predictions
out_dir_performance <- "./Output_MINT/Downstream_analysis/performance"
if (!dir.exists(out_dir_performance)) {
  dir.create(out_dir_performance, recursive = T)
}

#### crear curvas roc
out_file_performance <- c(
  file.path(out_dir_performance, "femalesContrast.jpeg"),
  file.path(out_dir_performance, "malesContrast.jpeg"),
  file.path(out_dir_performance, "diff_contrast.jpeg")
)
lapply(
  1:3,
  FUN = function(i)
    fun_auc(final_models[[i]], data = data_splitted_contrast[[i]], filename = out_file_performance[i])
)

##====Compute the final models with all subjects according to the tuning ####

models <- lapply(1:3, function(i)
  mint.splsda(
    X,
    data_ready$Y[[i]],
    ncomp = final_models[[i]]$ncomp,
    keepX = final_models[[i]]$keepX,
    study = study
  ))


##====score plots  =======
name_factor <- list(
  females = c("Control", "Other", "Case"),
  males = c("Control", "Other", "Case"),
  diff = c("DiffFemales", "DiffMales")
)
titles <- c("Females_Control_vs_Case_vs _Others",
            "Males_Control_vs_Case_vs_Others",
            "Diff")
pplots <- vector("list", length = 3)
for (i in 1:3) {
  model_i <- models[[i]]
  X_variates <- model_i$variates$X
  variance_explained <- model_i$prop_expl_var
  Y <- as.factor(model_i$Y)
  levels(Y) <- name_factor[[i]]
  variates <-  data.frame(X_variates, subjects  = Y)
  pplots[[i]] <- ggplot(variates, aes(comp1, comp2, color = subjects)) +
    geom_point() + xlab(paste("PC1", round(100 * variance_explained$X$`all data`[1]), "%")) +
    ylab(paste("PC2", round(100 * variance_explained$X$`all data`[2]), "%")) + ggtitle(titles[i])
  
  
}
## Save the plots
score_dir <- "./Output_MINT/Downstream_analysis/score_plots"
if (!dir.exists(score_dir)) {
  dir.create(score_dir)
}
lapply(1:3, function(x)
  ggsave(pplots[[x]], filename = file.path(score_dir, paste0(titles[x], ".jpeg"))))


##===== Venn Diagrams #####
feno <- data_ready$pheno_df
names(models) <- c("femalesHC", "malesHC", "diff_FemalesHC_MalesHC")

### ====Venn Diagram 1st contrast ====
## Reconstrucción
R <- models$femalesHC$variates$X %*% t(models$femalesHC$loadings$X)
R <- R[, colSums(R  != 0) > 0]
myPhenoData <- AnnotatedDataFrame(feno)
# Create ExpressionSet
myExpressionSet <- ExpressionSet(assayData = t(R), phenoData = myPhenoData)
pData(myExpressionSet)
fit <- diffexprs(myExpressionSet)
cont_ <- topTable(
  fit,
  coef = 1,
  adjust = "BH",
  sort = "none",
  n = Inf
)
## read all genes from MA
all_genes.path <- file.path("Results", "MA", names(models)[1], "all.genes.tsv")
sig_ma3 <- read.delim(all_genes.path)
## obtain sig genes DE ups and downs MA
ups_ma <- rownames(sig_ma3)[(sig_ma3$logFC > 0) &
                              (sig_ma3$p.adjust.fdr < 0.05)]
downs_ma <- rownames(sig_ma3)[(sig_ma3$logFC < 0) &
                                (sig_ma3$p.adjust.fdr < 0.05)]
## joint info
sig_ma3$ID <- rownames(sig_ma3)
cont_$ID <- rownames(cont_)
contrib <- plotLoadings(models$femalesHC, contrib = "max", method = "median")$X
contrib$ID <- rownames(contrib)
### first joint mint
mint_contrib <- merge(cont_, contrib, by = "ID")
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == 1, "Females_Case", NA)
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == 0,
                                     "Others",
                                     mint_contrib$GroupContrib_)
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == -1,
                                     "Females_Control",
                                     mint_contrib$GroupContrib_)
## genes by contrib

### extract ups downs por contribucion by group

#### Contrast 1: Females

#### UPS
ups_case <- mint_contrib[(mint_contrib$logFC > 0) &
                           (mint_contrib$adj.P.Val < 0.05) &
                           (mint_contrib$GroupContrib_ == "Females_Case"), ]

ups_control <- mint_contrib[(mint_contrib$logFC > 0) &
                              (mint_contrib$adj.P.Val < 0.05) &
                              (mint_contrib$GroupContrib_ == "Females_Control"), ]

ups_others <- mint_contrib[(mint_contrib$logFC > 0) &
                             (mint_contrib$adj.P.Val < 0.05) &
                             (mint_contrib$GroupContrib_ == "Others"), ]
#### DOWNS

dwons_other <- mint_contrib[(mint_contrib$logFC < 0) &
                              (mint_contrib$adj.P.Val < 0.05) &
                              (mint_contrib$GroupContrib_ == "Others"), ]

dwons_case <- mint_contrib[(mint_contrib$logFC < 0) &
                             (mint_contrib$adj.P.Val < 0.05) &
                             (mint_contrib$GroupContrib_ == "Females_Case"), ]

dwons_control <- mint_contrib[(mint_contrib$logFC < 0) &
                                (mint_contrib$adj.P.Val < 0.05) &
                                (mint_contrib$GroupContrib_ == "Females_Control"), ]

## print histogram for each group contribution

genes_desglosados <- list(
  UPS_case = ups_case,
  UPS_control = ups_control,
  UPS_other = ups_others,
  DOWNS_other = dwons_other,
  DOWNS_case = dwons_case,
  DOWNS_control = dwons_control
)

## hich groups does not have DE
hist(genes_desglosados$UPS_case$X.1)
hist(genes_desglosados$UPS_control$X1)
hist(genes_desglosados$UPS_other$X0)
hist(genes_desglosados$DOWNS_other$X0)
hist(genes_desglosados$DOWNS_case$X.1)
hist(genes_desglosados$DOWNS_control$X1)


## comparativa ups downs MINT y MA en general sin distincion de grupo
thresh <- 0
input_venn_1 <- list(
  UPS_MINT = c(
    genes_desglosados$UPS_case$ID[abs(genes_desglosados$UPS_case$X1) > thresh],
    genes_desglosados$UPS_control$ID[abs(genes_desglosados$UPS_control$X.1) > thresh],
    genes_desglosados$UPS_other$ID[abs(genes_desglosados$UPS_other$X0) > thresh]
  ),
  DOWNS_MINT = c(
    genes_desglosados$DOWNS_control$ID[abs(genes_desglosados$DOWNS_control$X.1) >
                                         thresh],
    genes_desglosados$DOWNS_other$ID[abs(genes_desglosados$DOWNS_other$X0) > thresh],
    genes_desglosados$DOWNS_case$ID[abs(genes_desglosados$DOWNS_case$X1) > thresh]
  ),
  UPS_MA = ups_ma,
  DOWNS_MA = downs_ma
)

input_venn_1.females <- input_venn_1
p1 <- ggVennDiagram(
  input_venn_1,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("MINT vs MA: Threshold", thresh))

#

input_venn_2 <-  list(
  UPS_Case = genes_desglosados$UPS_case$ID[abs(genes_desglosados$UPS_case$X1) >
                                             thresh],
  UPS_Control = genes_desglosados$UPS_control$ID[abs(genes_desglosados$UPS_control$X.1) > thresh],
  DOWNS_Control =  genes_desglosados$DOWNS_control$ID[abs(genes_desglosados$DOWNS_control$X.1) >
                                                        thresh],
  DOWNS_case =  genes_desglosados$DOWNS_case$ID[abs(genes_desglosados$DOWNS_case$X1) > thresh]
)

input_venn_2.females <- input_venn_2

### UPS y DOWNS por grupo MINT
p2 <- ggVennDiagram(
  input_venn_2,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("MINT, by group of interest: Threshold", thresh))

# titulo con threshold y DE LF y contraste
### Interseccion para ver concordancia

gg_combined <- ggarrange(p1, p2, common.legend = T, legend = "left")

gg_combined_with_title <- annotate_figure(
  gg_combined,
  top = text_grob(
    "Female Control vs Case p.adj <0.05",
    size = 12,
    face = "bold"
  ),
  left = NULL,
  right = NULL,
  bottom = NULL
)

gg_combined_with_title

dir.venn <- file.path("./Output_MINT/Downstream_analysis/VennDiagrams/femalesHC")
if (!dir.exists(dir.venn)) {
  dir.create(dir.venn, recursive = T)
}
ggsave(plot = gg_combined_with_title,
       filename = file.path(dir.venn, "Venn_Diagram_FemalesHC.jpeg"))

input_venn_3_up <- lapply(input_venn_1, function(x)
  intersect(x, ups_ma))

input_venn_3_down <- lapply(input_venn_1, function(x)
  intersect(x, downs_ma))

### ====Venn Diagram 2nd contrast ====
## Reconstrucción
R <- models$malesHC$variates$X %*% t(models$malesHC$loadings$X)
R <- R[, colSums(R  != 0) > 0]
myPhenoData <- AnnotatedDataFrame(feno)
# Create ExpressionSet
myExpressionSet <- ExpressionSet(assayData = t(R), phenoData = myPhenoData)
pData(myExpressionSet)
fit <- diffexprs(myExpressionSet)
cont_ <- topTable(
  fit,
  coef = 2,
  adjust = "BH",
  sort = "none",
  n = Inf
)
## read all genes from MA
all_genes.path <- file.path("Results", "MA", names(models)[2], "all.genes.tsv")
sig_ma3 <- read.delim(all_genes.path)
## obtain sig genes DE ups and downs MA
ups_ma <- rownames(sig_ma3)[(sig_ma3$logFC > 0) &
                              (sig_ma3$p.adjust.fdr < 0.05)]
downs_ma <- rownames(sig_ma3)[(sig_ma3$logFC < 0) &
                                (sig_ma3$p.adjust.fdr < 0.05)]
## joint info
sig_ma3$ID <- rownames(sig_ma3)
cont_$ID <- rownames(cont_)
contrib <- plotLoadings(models$malesHC, contrib = "max", method = "median")$X
contrib$ID <- rownames(contrib)
### first joint mint
mint_contrib <- merge(cont_, contrib, by = "ID")
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == 1, "Males_Case", NA)
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == 0,
                                     "Others",
                                     mint_contrib$GroupContrib_)
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == -1,
                                     "Males_Control",
                                     mint_contrib$GroupContrib_)
## genes by contrib

### extract ups downs por contribucion by group

#### Contrast 2: Males

#### UPS
ups_case <- mint_contrib[(mint_contrib$logFC > 0) &
                           (mint_contrib$adj.P.Val < 0.05) &
                           (mint_contrib$GroupContrib_ == "Males_Case"), ]

ups_control <- mint_contrib[(mint_contrib$logFC > 0) &
                              (mint_contrib$adj.P.Val < 0.05) &
                              (mint_contrib$GroupContrib_ == "Males_Control"), ]

ups_others <- mint_contrib[(mint_contrib$logFC > 0) &
                             (mint_contrib$adj.P.Val < 0.05) &
                             (mint_contrib$GroupContrib_ == "Others"), ]
#### DOWNS

dwons_other <- mint_contrib[(mint_contrib$logFC < 0) &
                              (mint_contrib$adj.P.Val < 0.05) &
                              (mint_contrib$GroupContrib_ == "Others"), ]

dwons_case <- mint_contrib[(mint_contrib$logFC < 0) &
                             (mint_contrib$adj.P.Val < 0.05) &
                             (mint_contrib$GroupContrib_ == "Males_Case"), ]

dwons_control <- mint_contrib[(mint_contrib$logFC < 0) &
                                (mint_contrib$adj.P.Val < 0.05) &
                                (mint_contrib$GroupContrib_ == "Males_Control"), ]

## print histogram for each group contribution

genes_desglosados <- list(
  UPS_case = ups_case,
  UPS_control = ups_control,
  UPS_other = ups_others,
  DOWNS_other = dwons_other,
  DOWNS_case = dwons_case,
  DOWNS_control = dwons_control
)

## hich groups does not have DE
hist(genes_desglosados$UPS_case$X.1)
hist(genes_desglosados$UPS_control$X1)
hist(genes_desglosados$UPS_other$X0)
hist(genes_desglosados$DOWNS_other$X0)
hist(genes_desglosados$DOWNS_case$X.1)
hist(genes_desglosados$DOWNS_control$X1)


## comparativa ups downs MINT y MA en general sin distincion de grupo
thresh <- 0
input_venn_1 <- list(
  UPS_MINT = c(
    genes_desglosados$UPS_case$ID[abs(genes_desglosados$UPS_case$X1) > thresh],
    genes_desglosados$UPS_control$ID[abs(genes_desglosados$UPS_control$X.1) > thresh],
    genes_desglosados$UPS_other$ID[abs(genes_desglosados$UPS_other$X0) > thresh]
  ),
  DOWNS_MINT = c(
    genes_desglosados$DOWNS_control$ID[abs(genes_desglosados$DOWNS_control$X.1) >
                                         thresh],
    genes_desglosados$DOWNS_other$ID[abs(genes_desglosados$DOWNS_other$X0) > thresh],
    genes_desglosados$DOWNS_case$ID[abs(genes_desglosados$DOWNS_case$X1) > thresh]
  ),
  UPS_MA = ups_ma,
  DOWNS_MA = downs_ma
)

input_venn_1.male <- input_venn_1
p1 <- ggVennDiagram(
  input_venn_1,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("MINT vs MA: Threshold", thresh))

#

input_venn_2 <-  list(
  UPS_Case = genes_desglosados$UPS_case$ID[abs(genes_desglosados$UPS_case$X1) >
                                             thresh],
  UPS_Control = genes_desglosados$UPS_control$ID[abs(genes_desglosados$UPS_control$X.1) > thresh],
  DOWNS_Control =  genes_desglosados$DOWNS_control$ID[abs(genes_desglosados$DOWNS_control$X.1) >
                                                        thresh],
  DOWNS_case =  genes_desglosados$DOWNS_case$ID[abs(genes_desglosados$DOWNS_case$X1) > thresh]
)


### UPS y DOWNS por grupo MINT
p2 <- ggVennDiagram(
  input_venn_2,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("MINT, by group of interest: Threshold", thresh))

# titulo con threshold y DE LF y contraste
### Interseccion para ver concordancia
input_venn_2.male <- input_venn_2

gg_combined <- ggarrange(p1, p2, common.legend = T, legend = "left")

gg_combined_with_title <- annotate_figure(
  gg_combined,
  top = text_grob(
    "Males: Control vs Case. p.adj <0.05",
    size = 16,
    face = "bold"
  ),
  left = NULL,
  right = NULL,
  bottom = NULL
)

gg_combined_with_title

dir.venn <- file.path("./Output_MINT/Downstream_analysis/VennDiagrams/malesHC")
if (!dir.exists(dir.venn)) {
  dir.create(dir.venn, recursive = T)
}
ggsave(plot = gg_combined_with_title,
       filename = file.path(dir.venn, "Venn_Diagram_MalesHC.jpeg"))
#####====common genes UPS====
### Interseccion para ver scatterplot
input_venn_3 <- lapply(input_venn_1, function(x)
  intersect(x, ups_ma))
### this plot just for checking
ggVennDiagram(
  input_venn_3,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  
) + ggtitle("UPS by group of interest: MINT")

### scatter plot of common genes UP and DOWN with logFC

joint_all.UPS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, ups_ma), c("ID", "logFC", "X.1", "X0", "X1")])
ups_ma.df <- sig_ma3[(sig_ma3$logFC > 0) &
                       (sig_ma3$p.adjust.fdr < 0.05), ]
downs_ma.df <- sig_ma3[(sig_ma3$logFC < 0) &
                         (sig_ma3$p.adjust.fdr < 0.05), ]

joint_all.DOWNS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, downs_ma), c("ID", "logFC", "X.1", "X0", "X1")])

joint_with_MA.UPS <- lapply(joint_all.UPS, function(X)
  inner_join(X, ups_ma.df, by = "ID"))

joint_with_MA.DOWNS <- lapply(joint_all.DOWNS, function(X)
  inner_join(X, downs_ma.df, by = "ID"))

## extract common UPS genes

common_UPS <- data.frame(
  MA = c(
    joint_with_MA.UPS$UPS_case$logFC.y,
    joint_with_MA.UPS$UPS_control$logFC.y,
    joint_with_MA.UPS$UPS_other$logFC.y
  ),
  MINT = c(
    joint_with_MA.UPS$UPS_case$logFC.x,
    joint_with_MA.UPS$UPS_control$logFC.x,
    joint_with_MA.UPS$UPS_other$logFC.x
  ),
  Category_UPS = as.factor(c(
    rep("UPS_case", length(joint_with_MA.UPS$UPS_case$logFC.y)),
    rep(
      "UPS_control",
      length(joint_with_MA.UPS$UPS_control$logFC.y)
    ),
    rep("UPS_other", length(joint_with_MA.UPS$UPS_other$logFC.y))
  ))
)

rownames(common_UPS) <- c(
  joint_with_MA.UPS$UPS_case$ID,
  joint_with_MA.UPS$UPS_control$ID,
  joint_with_MA.UPS$UPS_other$ID
)


corr_MA_MINT <- with(
  common_UPS,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_UPS, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_UPS),
    fill = (Category_UPS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between UPS LogFC"
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals

dir.concordance <- file.path("Output_MINT", "Concordance", "MalesHC", "UPS")
if (!dir.exists(dir.concordance)) {
  dir.create(dir.concordance, recursive = T)
}
ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "LinearRegression_complete.jpeg")
)


p1 <- grouped_ggscatterstats(common_UPS, MINT, MA, grouping.var = Category_UPS)
p1
ggsave(
  plot = p1,
  filename = file.path(dir.concordance, "LinearRegression_complete_bygroup.jpeg")
)
## in factor level order; probably better
# Outliers found!!!

####====Diagnosis Analysis=====
model_outliers <- lm(MINT ~ MA, data = common_UPS)
plot(fitted(model_outliers),
     abs(residuals(model_outliers)),
     xlab = "Predict values",
     ylab = "|Residuals|")
summary(lm(sqrt(abs(
  residuals(model_outliers)
)) ~ fitted(model_outliers)))
## varianza es constante.. creo
hist(residuals(model_outliers), main = paste("Shapiro Test: " , round(
  shapiro.test(model_outliers$residuals)$p.value, 10
)))
## leverage points
hatv <- hatvalues(model_outliers)
head(sort(hatv, decreasing = T))
p <- length(model_outliers$coefficients) # k+1
n <- length(model_outliers$fitted.values)
leverage.mean <- p / n # (k+1)/n
leverage_genes <- which(hatv > 2 * leverage.mean)
plot(hatv, type = "h")
abline(h = 2 * leverage.mean, col = "red")
genes <- row.names(common_UPS)
halfnorm(hatv,
         labs = genes,
         nlab = 4,
         ylab = "Leverage")
### Valores atipicos
stud <- rstudent(model_outliers)
head(sort(abs(stud), decreasing = TRUE))
outlier_genes <- which(abs(stud) > 2)
plot(stud, type = "h")
abline(h = -2, col = "red")
abline(h = 0)
abline(h = 2, col = "red")
grlib <- n - p - 1
outlier_genes2 <- which(abs(stud) > abs(qt(0.05 / (2 * n), grlib)))
## cobservaciones influyentes
cook <- cooks.distance(model_outliers)
genes_cook <- names(cook)
halfnorm(cook,
         nlab = 3,
         labs = genes_cook,
         ylab = "Distancia de Cook")
# Cook's D plot
# identify D values > 4/(n-k-1)
plot(model_outliers, which = 4)
abline(h = 4 / ((n - p - 2)), col = "red")
influcence_genes_thresh <- c(4 / ((n - p - 2)))
### valores leverage atipicos e influyentes:
influence_genes <- cook[cook > influcence_genes_thresh]

atypical_DE <- unique(c(
  names(influence_genes)
  ,
  names(outlier_genes2)
  ,
  names(outlier_genes)
  ,
  names(leverage_genes)
))
cook.ordered <- cook[names(cook) %in% atypical_DE]

cook_data <- data.frame(
  labels = names(cook.ordered),
  cooks_distance = cook.ordered,
  index = 1:length(cook.ordered)
)

# Create Cook's distance plot with ggplot
# Create Cook's distance plot with ggplot
p <- ggplot(cook_data, aes(x = factor(index, labels = labels), y = cooks_distance)) +
  geom_segment(aes(xend = factor(index, labels = labels), yend = 0), color = "blue") +
  labs(title = "Cook's Distances", x = "Observation", y = "Cook's Distance") + geom_hline(yintercept = influcence_genes_thresh,
                                                                                          colour = "red",
                                                                                          linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(plot = p,
       filename = file.path(dir.concordance, "Outliers_UP.jpeg"))
####=== Linear analysis without influence points ====


common_UPS <- data.frame(
  MA = c(
    joint_with_MA.UPS$UPS_case$logFC.y,
    joint_with_MA.UPS$UPS_control$logFC.y,
    joint_with_MA.UPS$UPS_other$logFC.y
  ),
  MINT = c(
    joint_with_MA.UPS$UPS_case$logFC.x,
    joint_with_MA.UPS$UPS_control$logFC.x,
    joint_with_MA.UPS$UPS_other$logFC.x
  ),
  Category_UPS = as.factor(c(
    rep("UPS_case", length(joint_with_MA.UPS$UPS_case$logFC.y)),
    rep(
      "UPS_control",
      length(joint_with_MA.UPS$UPS_control$logFC.y)
    ),
    rep("UPS_other", length(joint_with_MA.UPS$UPS_other$logFC.y))
  ))
)

rownames(common_UPS) <- c(
  joint_with_MA.UPS$UPS_case$ID,
  joint_with_MA.UPS$UPS_control$ID,
  joint_with_MA.UPS$UPS_other$ID
)


common_UPS$GENE <- rownames(common_UPS)
common_UPS_v2 <- common_UPS[rownames(common_UPS) %in% atypical_DE, ]

corr_MA_MINT <- with(
  common_UPS_v2,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_UPS_v2, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_UPS),
    fill = (Category_UPS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between LogFC UPS w.o Outliers"
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals

ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "LinearRegression_wo_outliers_UP.jpeg")
)

#####====common genes DOWNS====

### Interseccion para ver scatterplot
input_venn_3 <- lapply(input_venn_1, function(x)
  intersect(x, downs_ma))
ggVennDiagram(
  input_venn_3,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle("DOWNS by group of interest:MINT")

### scatter plot of common genes UP and DOWN with logFC

joint_all.UPS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, ups_ma), c("ID", "logFC", "X.1", "X0", "X1")])
ups_ma.df <- sig_ma3[(sig_ma3$logFC > 0) &
                       (sig_ma3$p.adjust.fdr < 0.05), ]
downs_ma.df <- sig_ma3[(sig_ma3$logFC < 0) &
                         (sig_ma3$p.adjust.fdr < 0.05), ]

joint_all.DOWNS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, downs_ma), c("ID", "logFC", "X.1", "X0", "X1")])

joint_with_MA.UPS <- lapply(joint_all.UPS, function(X)
  inner_join(X, ups_ma.df, by = "ID"))

joint_with_MA.DOWNS <- lapply(joint_all.DOWNS, function(X)
  inner_join(X, downs_ma.df, by = "ID"))

## extract common DOWNS genes

common_DOWNS <- data.frame(
  MA = c(
    joint_with_MA.DOWNS$DOWNS_case$logFC.y,
    joint_with_MA.DOWNS$DOWNS_control$logFC.y,
    joint_with_MA.DOWNS$DOWNS_other$logFC.y
  ),
  MINT = c(
    joint_with_MA.DOWNS$DOWNS_case$logFC.x,
    joint_with_MA.DOWNS$DOWNS_control$logFC.x,
    joint_with_MA.DOWNS$DOWNS_other$logFC.x
  ),
  Category_DOWNS = as.factor(c(
    rep(
      "DOWNS_case",
      length(joint_with_MA.DOWNS$DOWNS_case$logFC.y)
    ),
    rep(
      "DOWNS_control",
      length(joint_with_MA.DOWNS$DOWNS_control$logFC.y)
    ),
    rep(
      "DOWNS_other",
      length(joint_with_MA.DOWNS$DOWNS_other$logFC.y)
    )
  ))
)

rownames(common_DOWNS) <- c(
  joint_with_MA.DOWNS$DOWNS_other$ID,
  joint_with_MA.DOWNS$DOWNS_case$ID,
  joint_with_MA.DOWNS$DOWNS_control$ID
)

corr_MA_MINT <- with(
  common_DOWNS,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_DOWNS, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_DOWNS),
    fill = (Category_DOWNS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between DOWNS LogFC "
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals

q <- grouped_ggscatterstats(common_DOWNS, MINT, MA, grouping.var = Category_DOWNS)

dir.concordance <- file.path("Output_MINT", "Concordance", "MalesHC", "DOWNS")
if (!dir.exists(dir.concordance)) {
  dir.create(dir.concordance, recursive = T)
}

ggsave(
  plot = q,
  filename = file.path(dir.concordance,
                       "LinearRegression_complete_wo_outliers.jpeg")
)
ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance,
                       "LinearRegression_complete_bygroup_wo_outliers.jpeg")
)

# Outliers found!!!

####====Diagnosis Analysis=====
model_outliers <- lm(MINT ~ MA, data = common_DOWNS)
plot(fitted(model_outliers),
     abs(residuals(model_outliers)),
     xlab = "Predict values",
     ylab = "|Residuals|")
summary(lm(sqrt(abs(
  residuals(model_outliers)
)) ~ fitted(model_outliers)))
## varianza es constante.. creo
hist(residuals(model_outliers), main = paste("Shapiro Test: " , round(
  shapiro.test(model_outliers$residuals)$p.value, 10
)))
## leverage points
hatv <- hatvalues(model_outliers)
head(sort(hatv, decreasing = T))
p <- length(model_outliers$coefficients) # k+1
n <- length(model_outliers$fitted.values)
leverage.mean <- p / n # (k+1)/n
leverage_genes <- which(hatv > 2 * leverage.mean)
plot(hatv, type = "h")
abline(h = 2 * leverage.mean, col = "red")
genes <- row.names(common_DOWNS)
halfnorm(hatv,
         labs = genes,
         nlab = 4,
         ylab = "Leverage")
### Valores atipicos
stud <- rstudent(model_outliers)
head(sort(abs(stud), decreasing = TRUE))
outlier_genes <- which(abs(stud) > 2)
plot(stud, type = "h")
abline(h = -2, col = "red")
abline(h = 0)
abline(h = 2, col = "red")
grlib <- n - p - 1
outlier_genes2 <- which(abs(stud) > abs(qt(0.05 / (2 * n), grlib)))
## cobservaciones influyentes
cook <- cooks.distance(model_outliers)
names_cook_genes <- names(cook)
halfnorm(cook,
         nlab = 3,
         labs = names_cook_genes,
         ylab = "Distancia de Cook")
# Cook's D plot
# identify D values > 4/(n-k-1)
plot(model_outliers, which = 4)
abline(h = 4 / ((n - p - 2)), col = "red")
influcence_genes_thresh <- c(4 / ((n - p - 2)))
### valores leverage atipicos e influyentes:
influence_genes <- cook[cook > influcence_genes_thresh]

atypical_DE <- unique(c(
  names(influence_genes)
  ,
  names(outlier_genes2)
  ,
  names(outlier_genes)
  ,
  names(leverage_genes)
))
cook.ordered <- cook[names(cook) %in% atypical_DE]

cook_data <- data.frame(observation = names(cook.ordered),
                        cooks_distance = cook.ordered)

# Create Cook's distance plot with ggplot
cook.ordered <- cook[names(cook) %in% atypical_DE]

cook_data <- data.frame(
  labels = names(cook.ordered),
  cooks_distance = cook.ordered,
  index = 1:length(cook.ordered)
)

# Create Cook's distance plot with ggplot
# Create Cook's distance plot with ggplot
p <- ggplot(cook_data, aes(x = factor(index, labels = labels), y = cooks_distance)) +
  geom_segment(aes(xend = factor(index, labels = labels), yend = 0), color = "blue") +
  labs(title = "Cook's Distances", x = "Observation", y = "Cook's Distance") + geom_hline(yintercept = influcence_genes_thresh,
                                                                                          colour = "red",
                                                                                          linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(plot = p,
       filename = file.path(dir.concordance, "Outliers_DOWNS.jpeg"))

####=== Linear analysis without influence points ====
common_DOWNS$GENE <- rownames(common_DOWNS)
common_DOWNS_v2 <- common_DOWNS[rownames(common_DOWNS) %in% atypical_DE, ]

corr_MA_MINT <- with(
  common_DOWNS_v2,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_DOWNS_v2, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_DOWNS),
    fill = (Category_DOWNS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between LogFC DOWNS"
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals

ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "LinearRegression_DOWNS_wo_outliers.jpeg")
)

### ====Venn Diagram 3RD contrast ====
## Reconstrucción
R <- models$diff_FemalesHC_MalesHC$variates$X %*% t(models$diff_FemalesHC_MalesHC$loadings$X)
R <- R[, colSums(R  != 0) > 0]
myPhenoData <- AnnotatedDataFrame(feno)
# Create ExpressionSet
myExpressionSet <- ExpressionSet(assayData = t(R), phenoData = myPhenoData)
pData(myExpressionSet)
fit <- diffexprs(myExpressionSet)
cont_ <- topTable(
  fit,
  coef = 3,
  adjust = "BH",
  sort = "none",
  n = Inf
)
## read all genes from MA
all_genes.path <- file.path("Results", "MA", names(models)[3], "all.genes.tsv")
sig_ma3 <- read.delim(all_genes.path)
## obtain sig genes DE ups and downs MA
ups_ma <- rownames(sig_ma3)[(sig_ma3$logFC > 0) &
                              (sig_ma3$p.adjust.fdr < 0.05)]
downs_ma <- rownames(sig_ma3)[(sig_ma3$logFC < 0) &
                                (sig_ma3$p.adjust.fdr < 0.05)]
## joint info
sig_ma3$ID <- rownames(sig_ma3)
cont_$ID <- rownames(cont_)
contrib <- plotLoadings(models$diff_FemalesHC_MalesHC,
                        contrib = "max",
                        method = "median")$X
contrib$ID <- rownames(contrib)
### first joint mint
mint_contrib <- merge(cont_, contrib, by = "ID")
mint_contrib$GroupContrib_ <- ifelse(mint_contrib$GroupContrib == -1, "FemalesHC", "MalesHC")
## genes by contrib

### extract ups downs por contribucion by group

#### Contrast 2: Males

#### UPS
ups_case <- mint_contrib[(mint_contrib$logFC > 0) &
                           (mint_contrib$adj.P.Val < 0.05) &
                           (mint_contrib$GroupContrib_ == "FemalesHC"), ]

ups_control <- mint_contrib[(mint_contrib$logFC > 0) &
                              (mint_contrib$adj.P.Val < 0.05) &
                              (mint_contrib$GroupContrib_ == "MalesHC"), ]

#### DOWNS


dwons_case <- mint_contrib[(mint_contrib$logFC < 0) &
                             (mint_contrib$adj.P.Val < 0.05) &
                             (mint_contrib$GroupContrib_ == "FemalesHC"), ]

dwons_control <- mint_contrib[(mint_contrib$logFC < 0) &
                                (mint_contrib$adj.P.Val < 0.05) &
                                (mint_contrib$GroupContrib_ == "MalesHC"), ]

## print histogram for each group contribution

genes_desglosados <- list(
  UPS_femalesHC = ups_case,
  UPS_malesHC = ups_control,
  DOWNS_femalesHC = dwons_case,
  DOWNS_malesHC = dwons_control
)

## hich groups does not have DE
hist(genes_desglosados$UPS_femalesHC$X.1)
hist(genes_desglosados$UPS_malesHC$X1)
hist(genes_desglosados$DOWNS_femalesHC$X.1)
hist(genes_desglosados$DOWNS_malesHC$X1)


## comparativa ups downs MINT y MA en general sin distincion de grupo
thresh <- 0
input_venn_1 <- list(
  UPS_MINT = c(
    genes_desglosados$UPS_femalesHC$ID[abs(genes_desglosados$UPS_femalesHC$X1) > thresh],
    genes_desglosados$UPS_malesHC$ID[abs(genes_desglosados$UPS_malesHC$X.1) > thresh]
  ),
  DOWNS_MINT = c(
    genes_desglosados$DOWNS_femalesHC$ID[abs(genes_desglosados$DOWNS_femalesHC$X.1) >
                                           thresh],
    genes_desglosados$DOWNS_malesHC$ID[abs(genes_desglosados$DOWNS_malesHC$X1) > thresh]
  ),
  UPS_MA = ups_ma,
  DOWNS_MA = downs_ma
)

input_venn_1.diff <- input_venn_1

p1 <- ggVennDiagram(
  input_venn_1,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("MINT vs MA: Threshold", thresh))

#

input_venn_2 <-  list(
  UPS_femalesHC = genes_desglosados$UPS_femalesHC$ID[abs(genes_desglosados$UPS_femalesHC$X1) >
                                                       thresh],
  UPS_malesHC = genes_desglosados$UPS_malesHC$ID[abs(genes_desglosados$UPS_malesHC$X.1) > thresh],
  DOWNS_femalesHC =  genes_desglosados$DOWNS_femalesHC$ID[abs(genes_desglosados$DOWNS_femalesHC$X.1) >
                                                            thresh],
  DOWNS_malesHC =  genes_desglosados$DOWNS_malesHC$ID[abs(genes_desglosados$DOWNS_malesHC$X1) > thresh]
)


### UPS y DOWNS por grupo MINT
p2 <- ggVennDiagram(
  input_venn_2,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("MINT, by group of interest: Threshold", thresh))

# titulo con threshold y DE LF y contraste
### Interseccion para ver concordancia
input_venn_2.diff <- input_venn_2

gg_combined <- ggarrange(p1, p2, common.legend = T, legend = "left")

gg_combined_with_title <- annotate_figure(
  gg_combined,
  top = text_grob(
    "Diff FemalesHC vs MalesHC p.adj <0.05",
    size = 12,
    face = "bold"
  ),
  left = NULL,
  right = NULL,
  bottom = NULL
)

gg_combined_with_title
dir.venn <- file.path("./Output_MINT/Downstream_analysis/VennDiagrams/diff_femalesHC_malesHC")
if (!dir.exists(dir.venn)) {
  dir.create(dir.venn, recursive = T)
}
ggsave(plot = gg_combined_with_title,
       filename = file.path(dir.venn, "venn_diagram.jpeg"))
#####====common genes UPS====
### Interseccion para ver scatterplot
input_venn_3 <- lapply(input_venn_1, function(x)
  intersect(x, ups_ma))

ggVennDiagram(
  input_venn_3,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2

) + ggtitle("UPS by group of interest:MINT")

### scatter plot of common genes UP and DOWN with logFC

joint_all.UPS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, ups_ma), c("ID", "logFC", "X.1", "X1")])
ups_ma.df <- sig_ma3[(sig_ma3$logFC > 0) &
                       (sig_ma3$p.adjust.fdr < 0.05), ]
downs_ma.df <- sig_ma3[(sig_ma3$logFC < 0) &
                         (sig_ma3$p.adjust.fdr < 0.05), ]

joint_all.DOWNS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, downs_ma), c("ID", "logFC", "X.1", "X1")])

joint_with_MA.UPS <- lapply(joint_all.UPS, function(X)
  inner_join(X, ups_ma.df, by = "ID"))

joint_with_MA.DOWNS <- lapply(joint_all.DOWNS, function(X)
  inner_join(X, downs_ma.df, by = "ID"))

## extract common UPS genes

common_UPS <- data.frame(
  MA = c(
    joint_with_MA.UPS$UPS_femalesHC$logFC.y,
    joint_with_MA.UPS$UPS_malesHC$logFC.y
  ),
  MINT = c(
    joint_with_MA.UPS$UPS_femalesHC$logFC.x,
    joint_with_MA.UPS$UPS_malesHC$logFC.x
  ),
  Category_UPS = as.factor(c(
    rep(
      "UPS_femalesHC",
      length(joint_with_MA.UPS$UPS_femalesHC$logFC.y)
    ), rep(
      "UPS_malesHC",
      length(joint_with_MA.UPS$UPS_malesHC$logFC.y)
    )
  ))
)

rownames(common_UPS) <- c(joint_with_MA.UPS$UPS_femalesHC$ID,
                          joint_with_MA.UPS$UPS_malesHC$ID)

corr_MA_MINT <- with(
  common_UPS,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_UPS, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_UPS),
    fill = (Category_UPS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between UPS LogFC"
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals

dir.concordance <- file.path("Output_MINT", "Concordance", "diff_FemalesHC_MalesHC", "UPS")
if (!dir.exists(dir.concordance)) {
  dir.create(dir.concordance, recursive = T)
}

ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "linear_regression_UPS.jpeg")
)
q <- grouped_ggscatterstats(common_UPS, MINT, MA, grouping.var = Category_UPS)
q
ggsave(
  plot = q,
  filename = file.path(dir.concordance, "Grouped_linear_regression_UPS.jpeg")
)

# Outliers found!!!

####====Diagnosis Analysis=====
model_outliers <- lm(MINT ~ MA, data = common_UPS)
plot(fitted(model_outliers),
     abs(residuals(model_outliers)),
     xlab = "Predict values",
     ylab = "|Residuals|")
summary(lm(sqrt(abs(
  residuals(model_outliers)
)) ~ fitted(model_outliers)))
## varianza es constante.. creo
hist(residuals(model_outliers), main = paste("Shapiro Test: " , round(
  shapiro.test(model_outliers$residuals)$p.value, 10
)))
## leverage points
hatv <- hatvalues(model_outliers)
head(sort(hatv, decreasing = T))
p <- length(model_outliers$coefficients) # k+1
n <- length(model_outliers$fitted.values)
leverage.mean <- p / n # (k+1)/n
leverage_genes <- which(hatv > 2 * leverage.mean)
plot(hatv, type = "h")
abline(h = 2 * leverage.mean, col = "red")
genes <- row.names(common_UPS)
halfnorm(hatv,
         labs = genes,
         nlab = 4,
         ylab = "Leverage")
### Valores atipicos
stud <- rstudent(model_outliers)
head(sort(abs(stud), decreasing = TRUE))
outlier_genes <- which(abs(stud) > 2)
plot(stud, type = "h")
abline(h = -2, col = "red")
abline(h = 0)
abline(h = 2, col = "red")
grlib <- n - p - 1
outlier_genes2 <- which(abs(stud) > abs(qt(0.05 / (2 * n), grlib)))
## cobservaciones influyentes
cook <- cooks.distance(model_outliers)
genes_cook <- names(cook)
halfnorm(cook,
         nlab = 3,
         labs = genes_cook,
         ylab = "Distancia de Cook")
# Cook's D plot
# identify D values > 4/(n-k-1)
plot(model_outliers, which = 4)
abline(h = 4 / ((n - p - 2)), col = "red")
influcence_genes_thresh <- c(4 / ((n - p - 2)))
### valores leverage atipicos e influyentes:
influence_genes <- cook[cook > influcence_genes_thresh]

atypical_DE <- unique(c(
  names(influence_genes)
  ,
  names(outlier_genes2)
  ,
  names(outlier_genes)
  ,
  names(leverage_genes)
))
# Create Cook's distance plot with ggplot
cook.ordered <- cook[names(cook) %in% atypical_DE]

cook_data <- data.frame(
  labels = names(cook.ordered),
  cooks_distance = cook.ordered,
  index = 1:length(cook.ordered)
)

# Create Cook's distance plot with ggplot
# Create Cook's distance plot with ggplot
p <- ggplot(cook_data, aes(x = factor(index, labels = labels), y = cooks_distance)) +
  geom_segment(aes(xend = factor(index, labels = labels), yend = 0), color = "blue") +
  labs(title = "Cook's Distances", x = "Observation", y = "Cook's Distance") + geom_hline(yintercept = influcence_genes_thresh,
                                                                                          colour = "red",
                                                                                          linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(plot = p,
       filename = file.path(dir.concordance, "Outliers_UPS.jpeg"))

####=== Linear analysis without influence points ====




common_UPS$GENE <- rownames(common_UPS)
common_UPS_v2 <- common_UPS[rownames(common_UPS) %in% atypical_DE, ]

corr_MA_MINT <- with(
  common_UPS_v2,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_UPS_v2, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_UPS),
    fill = (Category_UPS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between LogFC\n UPS w.o Outliers"
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals
ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "Lienar_Regression_WO_Outliers_UPS.jpeg")
)
#####====common genes DOWNS====

### Interseccion para ver scatterplot
input_venn_3 <- lapply(input_venn_1, function(x)
  intersect(x, downs_ma))
ggVennDiagram(
  input_venn_3,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle("DOWNS by group of interest:MINT")

### scatter plot of common genes UP and DOWN with logFC

joint_all.UPS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, ups_ma), c("ID", "logFC", "X.1", "X1")])
ups_ma.df <- sig_ma3[(sig_ma3$logFC > 0) &
                       (sig_ma3$p.adjust.fdr < 0.05), ]
downs_ma.df <- sig_ma3[(sig_ma3$logFC < 0) &
                         (sig_ma3$p.adjust.fdr < 0.05), ]

joint_all.DOWNS <- lapply(genes_desglosados, function(x)
  x[x$ID %in% intersect(x$ID, downs_ma), c("ID", "logFC", "X.1", "X1")])

joint_with_MA.UPS <- lapply(joint_all.UPS, function(X)
  inner_join(X, ups_ma.df, by = "ID"))

joint_with_MA.DOWNS <- lapply(joint_all.DOWNS, function(X)
  inner_join(X, downs_ma.df, by = "ID"))

## extract common DOWNS genes

common_DOWNS <- data.frame(
  MA = c(
    joint_with_MA.DOWNS$DOWNS_femalesHC$logFC.y,
    joint_with_MA.DOWNS$DOWNS_malesHC$logFC.y
  ),
  MINT = c(
    joint_with_MA.DOWNS$DOWNS_femalesHC$logFC.x,
    joint_with_MA.DOWNS$DOWNS_malesHC$logFC.x
  ),
  Category_DOWNS = as.factor(c(
    rep(
      "DOWNS_femalesHC",
      length(joint_with_MA.DOWNS$DOWNS_femalesHC$logFC.y)
    ),
    rep(
      "DOWNS_malesHC",
      length(joint_with_MA.DOWNS$DOWNS_malesHC$logFC.y)
    )
  ))
)

rownames(common_DOWNS) <- c(joint_with_MA.DOWNS$DOWNS_femalesHC$ID,
                            joint_with_MA.DOWNS$DOWNS_malesHC$ID)

corr_MA_MINT <- with(
  common_DOWNS,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_DOWNS, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_DOWNS),
    fill = (Category_DOWNS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between DOWNS LogFC "
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals

dir.concordance <- file.path("Output_MINT",
                             "Concordance",
                             "diff_FemalesHC_MalesHC",
                             "DOWNS")
if (!dir.exists(dir.concordance)) {
  dir.create(dir.concordance, recursive = T)
}
ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "Linear_Regression_DOWNS.jpeg")
)

q <- grouped_ggscatterstats(common_DOWNS, MINT, MA, grouping.var = Category_DOWNS)
q
ggsave(
  plot = q,
  filename = file.path(dir.concordance, "Linear_Regression_DOWNS_GROUPED.jpeg")
)

##aqui tambien
# Outliers found!!!

####====Diagnosis Analysis=====
model_outliers <- lm(MINT ~ MA, data = common_DOWNS)
plot(fitted(model_outliers),
     abs(residuals(model_outliers)),
     xlab = "Predict values",
     ylab = "|Residuals|")
summary(lm(sqrt(abs(
  residuals(model_outliers)
)) ~ fitted(model_outliers)))
## varianza es constante.. creo
hist(residuals(model_outliers), main = paste("Shapiro Test: " , round(
  shapiro.test(model_outliers$residuals)$p.value, 10
)))
## leverage points
hatv <- hatvalues(model_outliers)
head(sort(hatv, decreasing = T))
p <- length(model_outliers$coefficients) # k+1
n <- length(model_outliers$fitted.values)
leverage.mean <- p / n # (k+1)/n
leverage_genes <- which(hatv > 2 * leverage.mean)
plot(hatv, type = "h")
abline(h = 2 * leverage.mean, col = "red")
genes <- row.names(common_DOWNS)
halfnorm(hatv,
         labs = genes,
         nlab = 4,
         ylab = "Leverage")
### Valores atipicos
stud <- rstudent(model_outliers)
head(sort(abs(stud), decreasing = TRUE))
outlier_genes <- which(abs(stud) > 2)
plot(stud, type = "h")
abline(h = -2, col = "red")
abline(h = 0)
abline(h = 2, col = "red")
grlib <- n - p - 1
outlier_genes2 <- which(abs(stud) > abs(qt(0.05 / (2 * n), grlib)))
## cobservaciones influyentes
cook <- cooks.distance(model_outliers)
names_cook_genes <- names(cook)
halfnorm(cook,
         nlab = 3,
         labs = names_cook_genes,
         ylab = "Distancia de Cook")
# Cook's D plot
# identify D values > 4/(n-k-1)
plot(model_outliers, which = 4)
abline(h = 4 / ((n - p - 2)), col = "red")
influcence_genes_thresh <- c(4 / ((n - p - 2)))
### valores leverage atipicos e influyentes:
influence_genes <- cook[cook > influcence_genes_thresh]

atypical_DE <- unique(c(
  names(influence_genes)
  ,
  names(outlier_genes2)
  ,
  names(outlier_genes)
  ,
  names(leverage_genes)
))
# Create Cook's distance plot with ggplot
cook.ordered <- cook[names(cook) %in% atypical_DE]

cook_data <- data.frame(
  labels = names(cook.ordered),
  cooks_distance = cook.ordered,
  index = 1:length(cook.ordered)
)

# Create Cook's distance plot with ggplot
p <- ggplot(cook_data, aes(x = factor(index, labels = labels), y = cooks_distance)) +
  geom_segment(aes(xend = factor(index, labels = labels), yend = 0), color = "blue") +
  labs(title = "Cook's Distances", x = "Observation", y = "Cook's Distance") + geom_hline(yintercept = influcence_genes_thresh,
                                                                                          colour = "red",
                                                                                          linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(plot = p,
       filename = file.path(dir.concordance, "Outliers_DOWNS.jpeg"))

####=== Linear analysis without influence points ====
common_DOWNS$GENE <- rownames(common_DOWNS)
common_DOWNS_v2 <- common_DOWNS[rownames(common_DOWNS) %in% atypical_DE, ]

corr_MA_MINT <- with(
  common_DOWNS_v2,
  corr.test(
    MINT,
    MA,
    use = "pairwise",
    method = "pearson",
    adjust = "holm",
    alpha = .05,
    ci = TRUE,
    minlength = 5,
    normal = TRUE
  )
)
CI = round(corr_MA_MINT$ci.adj, 4)
PVAL = round(corr_MA_MINT$p.adj, 4)
N = corr_MA_MINT$n
T_ = round(corr_MA_MINT$t, 4)

p <- ggplot(common_DOWNS_v2, aes(MINT, MA)) +
  geom_point(aes(
    colour = (Category_DOWNS),
    fill = (Category_DOWNS)
  ), shape = 21, size = 1) +
  geom_smooth(
    method = lm,
    na.rm = TRUE,
    fullrange = TRUE,
    aes(group = 1),
    colour = "black"
  ) + ggtitle(
    paste(
      "CI =",
      CI[1],
      ",",
      CI[2],
      ",",
      CI[3],
      "\np.val =",
      PVAL,
      "N =",
      N,
      "T =",
      T_,
      "\nLinear Regression between LogFC DOWNS"
    )
  )

p_with_marginals <- ggMarginal(p, type = "histogram", bins = 30)
p_with_marginals
ggsave(
  plot = p_with_marginals,
  filename = file.path(dir.concordance, "LR_wo_Outliers_DOWNS.jpeg")
)

#==== Other Venn Diagrams ======
##===FEMALES=====
## UPS diagram diff women

input_list <- list(UPS_Diff = input_venn_2.diff$UPS_femalesHC,
                   UPS_case = input_venn_2.females$UPS_Case,
                   UPS_Control = input_venn_2.females$UPS_Control)



p1 <- ggVennDiagram(
  input_list,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("Up DE Genes", thresh))
p1

## downs diagram diff women

input_list <- list(DOWNS_Diff = input_venn_2.diff$DOWNS_femalesHC,
                   DOWNS_case = input_venn_2.females$DOWNS_case,
                   DOWNSControl = input_venn_2.females$DOWNS_Control)



p2 <- ggVennDiagram(
  input_list,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("Down DE  Genes", thresh))
p2

gg_combined <- ggarrange(p1, p2, common.legend = T, legend = "left")

gg_combined_with_title <- annotate_figure(
  gg_combined,
  top = text_grob(
    "Females: p.adj <0.05, threshold = 0. \n by MINT",
    size = 16,
    face = "bold"
  ),
  left = NULL,
  right = NULL,
  bottom = NULL
)
gg_combined_with_title

ggsave(plot = gg_combined_with_title,
       filename = "./Output_MINT/Downstream_analysis/VennDiagrams/femalesHC/UPS_DOWNS_females.jpeg")


##===MALES=====
input_list <- list(DOWNS_Diff = input_venn_2.diff$DOWNS_malesHC,
                   DOWNS_case = input_venn_2.male$DOWNS_case,
                   DOWNS_Control = input_venn_2.male$DOWNS_Control)
p1 <- ggVennDiagram(
  input_list,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("Down Genes", thresh))
p2
## UPS diagram diff males


input_list <- list(UPS_Diff = input_venn_2.diff$UPS_malesHC,
                   UPS_case = input_venn_2.male$UPS_Case,
                   UPS_Control = input_venn_2.male$UPS_Control)

p2 <- ggVennDiagram(
  input_list,
  title = "Venn Diagram",
  show.plot = TRUE,
  label = "both",
  set_size = 2,
  label_size = 3,
  label_percent_digit = 2
  # Para mostrar el gráfico en la ventana de gráficos
) + ggtitle(paste("Up Genes", thresh))
p2
gg_combined <- ggarrange(p1, p2, common.legend = T, legend = "left")

gg_combined_with_title <- annotate_figure(
  gg_combined,
  top = text_grob(
    "Males: p.adj <0.05, threshold = 0. \n by MINT",
    size = 16,
    face = "bold"
  ),
  left = NULL,
  right = NULL,
  bottom = NULL
)

ggsave(plot = gg_combined_with_title,
       filename = "./Output_MINT/Downstream_analysis/VennDiagrams/malesHC/UPS_male.jpeg")
## downs diagram diff women

