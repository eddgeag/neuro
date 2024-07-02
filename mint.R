library(mixOmics)
library(ggplot2)
library(limma)
library(dplyr)
library(Biobase)

#====Define Functions====
### function to load data and prepare
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
### prepare the data for the analysis
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

### run MINT
fun_mint <- function(X, Y, study, contrast) {
  ### First we need to set the optimal number of component
  ### set 20, afterwards, the algorithm will choose
  model_1 <- mint.splsda(X,
                         Y,
                         ncomp = 20,
                         study = study,
                         scale = T)
  model_1.perf <- perf(model_1)
  perf_dir <- "./Output_MINT/performance/tuning"
  
  if (!dir.exists(perf_dir)) {
    dir.create(perf_dir, recursive = T)
  }
  jpeg(filename = file.path(perf_dir, paste0(contrast, ".jpeg")))
  plot(model_1.perf)
  dev.off()
  ### set the maximium number of componentes chosen
  ncomp <- max(model_1.perf$choice.ncomp)
  ### if the component chosen is 1 set to two
  ### maybe will reduce the performance ?
  ### Thinking to change this
  if (ncomp < 2) {
    ncomp <- 2
  }
  
  p <- ncol(X)
  ## set parameter sparsity 
  ## grid search from 500 to the total genes by steps of 1000 
  secuencia <- seq(500, p, 1000)
  tune_model <- tune.mint.splsda(
    X,
    Y,
    ncomp = ncomp,
    study = study,
    scale = T,
    test.keepX = secuencia,
    progressBar = T,
    auc = T,
    tol = 1e-06,
    max.iter = 1000
  )
  ## here we were planning to search for a better grid
  # p.max <- max(tune_model$choice.keepX)
  # p.min <- min(tune_model$choice.keepX)
  # secuenica <- seq(10, p.max, round(p.min / 10))
  # tune_model_2 <- tune.mint.splsda(
  #   X,
  #   Y,
  #   ncomp = ncomp,
  #   study = study,
  #   scale = T,
  #   test.keepX = secuencia,
  #   progressBar = T,
  #   auc = T,
  #  max.iter = 1000
  # )
  ### final mint model
  final_mint <- mint.splsda(
    X,
    Y,
    ncomp = ncomp,
    study = study,
    keepX = tune_model$choice.keepX, ##  optimal number of variables
    scale = T
  )
  
  
  
  return(final_mint)
}

### function to split on train and test
fun_split <- function(X, Y, study, p) {
  n <- nrow(X)
  idx <- sample(n, round(n * p))
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

#====Downstream analysis======
## paths to metanalysis
path_data <- "./metanalysis/data/SN_datasets.RData"
path_results_metanalysis <- "./Results/DEA/"

## call the loading data functions
pre_data <- fun_load(path_data = path_data, path_results_metanalysis = path_results_metanalysis)
data_ready <- prepare_data(pre_data$datps, pre_data$metanalisis_results)
input_parafac <- data_ready$X

x1 <- data_ready$X[[1]]
x2 <- data_ready$X[[2]]
x3 <- data_ready$X[[3]]
x4 <- data_ready$X[[4]]
data_list <- list(x1, x2, x3, x4)

X <- t(as.matrix(Reduce("cbind", data_list))) # Concatenate the studies

Y <- data_ready$Y
set.seed(123456) ## set seed for reproducibility
study <- data_ready$study
## for each contrast split the data
data_splitted_contrast <- lapply(1:3, function(y)
  fun_split(X, data_ready$Y[[y]], study = study, p = 0.7))
final_models <-  vector("list", 3)


contrastnames <- c("Females", "Males", "Diff")
## compute the final models for each contrast
final_models <- lapply(1:3, function(x)
  fun_mint(
    X = data_splitted_contrast[[x]]$train$X,
    Y = data_splitted_contrast[[x]]$train$Y,
    study = data_splitted_contrast[[x]]$train$study,
    contrast = contrastnames[x]
  ))


## save the mddels
saveRDS(final_models, "./final_models_1.rds")
