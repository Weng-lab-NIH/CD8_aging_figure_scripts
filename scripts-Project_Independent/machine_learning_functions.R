
Generate_Multivariate_Regression <- function(
  matrix_train, 
  labels_train, 
  model_type, 
  alpha = 0.5, 
	n_cores, 
	path_summary = './data/interim/machine_learning/multivariate/', 
	name_model = '..') {
  
	require(biglasso)
	require(tidyverse)
	require(MASS)
	require(Matrix)

  version <- '6.3.2'
  dir.create(path_summary)
  
  # Generate a big.matrix object for each set:
  p.bm.train <- as.big.matrix(matrix_train)

  # PERFORM REGRESSION ANALYSIS:
  
  # ElasticNet Model:
  if (model_type == 'ENET') {
  ptm <- proc.time()
  enet <- cv.biglasso(p.bm.train, labels_train, penalty = 'enet', family = 'gaussian', alpha = alpha, seed = 1234, nfolds = 10, ncores = n_cores)
  proc.time() - ptm
  enet$lambda.min
  saveRDS(enet, paste(path_summary,version,name_model,'cvEnetModel_.Rds', sep = '_'))

  png(paste(path_summary,version,name_model,'cvEnetModel_.png', sep = '_'), height = 500, width = 500, units = 'px')
  par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
  plot(enet, type = "all")
  dev.off()
  summary(enet)
  
  # LASSO Model:
  } else if (model_type == 'RIDGE') {
  ptm <- proc.time()
  ridge <- cv.biglasso(p.bm.train, l.train, penalty = 'ridge', family = 'gaussian', seed = 1234, nfolds = 10, ncores = n_cores)
  proc.time() - ptm
  ridge$lambda.min
  saveRDS(ridge, paste(path_summary,version,name_model,'cvRidgeModel_.Rds', sep = '_'))
  
  png(paste(path_summary,version,name_model,'cvRidgeModel_.png', sep = '_'), height = 500, width = 500, units = 'px')
  par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
  plot(ridge, type = "all")
  dev.off()
  summary(ridge)
  }
  
  if (model_type == 'ENET') {
    return(enet)
  } else if (model_type == 'RIDGE') {
    return(ridge)
  }
}

Generate_Multivariate_Regression_v2 <- function(
  matrix_train, 
  labels_train, 
  model_type, 
  alpha = 0.5, 
  n_cores, 
  path_summary = './data/interim/machine_learning/multivariate/', 
  name_model = '..') {
  
  require(biglasso)
  require(tidyverse)
  require(MASS)
  require(Matrix)
  
  version <- '2'
  dir.create(path_summary)
  
  # Generate a big.matrix object for each set:
  p.bm.train <- as.big.matrix(matrix_train)
  
  # PERFORM REGRESSION ANALYSIS:
  
  # ElasticNet Model:
  if (model_type == 'ENET') {
    ptm <- proc.time()
    enet <- biglasso(p.bm.train, labels_train, penalty = 'enet', family = 'gaussian', alpha = alpha, seed = 1234, ncores = n_cores, screen = 'SSR-BEDPP')
    proc.time() - ptm
    #enet$lambda.min
    #saveRDS(enet, paste(path_summary,version,name_model,'cvEnetModel_.Rds', sep = '_'))
    
    # LASSO Model:
  } else if (model_type == 'RIDGE') {
    ptm <- proc.time()
    ridge <- biglasso(p.bm.train, l.train, penalty = 'ridge', family = 'gaussian', seed = 1234, ncores = n_cores, screen = 'SSR-BEDPP')
    proc.time() - ptm
    
  }
  
  if (model_type == 'ENET') {
    return(enet)
  } else if (model_type == 'RIDGE') {
    return(ridge)
  }
}

Generate_Multivariate_CrossValidation_v2 <- function(
  matrix_train, 
  labels_train, 
  model_type, 
  alpha = 0.5, 
  n_cores, 
  n_fold, 
  name_model = '..') {
  
  require(biglasso)
  require(tidyverse)
  require(MASS)
  require(Matrix)
  require(matrixStats)
  require(keras)
  
  version <- '2'
  
  # Generate index for tracking data in cross-validation:
  index <- c(1:nrow(matrix_train))
  indexr <- split(sample(index,length(index)), c(1:n_fold))
  
  #Generate objects to hold results:
  models <- list()
  cross_validation_predictions <- tibble(y=labels_train, index=index)
  lambda_all <- c()
  
  for (i in c(1:n_fold)) {
    
    index_train <- setdiff(index, indexr[[i]])
    
    train <- matrix_train[index_train,]
    y = labels_train[index_train]
    oob <- matrix_train[indexr[[i]],]
    
    train_data <- scale(train)
    mean_train <- colMeans(train_data)
    sd_train <- colSds(train_data)
    test_data <- scale(oob, center = mean_train, scale = sd_train)
    
    train_big <- as.big.matrix(train_data)
    test_big <- as.big.matrix(test_data)
    
    model <- biglasso(X = train_big, y = y, 
                      
                      penalty = model_type,
                      family = "gaussian", 
                      screen = "SSR-BEDPP", 
                      ncores = n_cores, 
                      alpha = alpha, 
                      nlambda = 100,
                      output.time = T, 
                      return.time = T, 
                      verbose = T)
    
    
    models[[i]] <- model
    lambda <- model$lambda[which.min(model$loss)]
    
    oob_prediction <- predict(object = model, 
                              X = test_big,
                              type = 'response',
                              lambda = lambda)
    
    output_oob <- tibble(index = indexr[[i]], 
                         one = oob_prediction@x)
    colnames(output_oob)[2] <- i
    cross_validation_predictions <- left_join(cross_validation_predictions,
                                              output_oob, 
                                              by = 'index')
    lambda_all <- c(lambda_all,lambda)
  }
  
  enet <- list(models, cross_validation_predictions, lambda_all)
 
  return(enet)
}


Evaluate_Multivariate <- function(model, matrix_test, test_labels, metadata, 
	path_summary = './data/interim/machine_learning/multivariate/', 
	name_model = '..') {
  
  require(biglasso)
  require(tidyverse)
  require(MASS)
  require(Matrix)
  require(Metrics)
  
  version <- '6.3.2'
  dir.create(path_summary)
  
  # Generate a big.matrix object for each set:
  p.bm.test <- as.big.matrix(matrix_test)
  test_index <- rownames(matrix_test)
  
  # Generate Summary Statistics:
  # Coefficients:
  coefficients <- coef(model)
  coeff_names <- coefficients@Dimnames[[1]]
  df_coef <- data.frame(gene = coeff_names, coef = coefficients)
  write.csv(df_coef, paste(path_summary,version,'_multivariate_coefficients.csv', sep = ''))
  
  # Predict the ages of cells using the chosen model:
  test_predictions <- predict(model, p.bm.test, type = 'response')
  fin_test <- metadata[test_index,]
  fin_test <- fin_test %>% mutate(tp = test_predictions[ , 1], labels = test_labels)
  write.csv(fin_test, paste(path_summary,version,'_multivariate_predictions.csv', sep = ''))
  
  fin_testt <- fin_test %>%
    group_by(Age) %>%
    filter(tp < quantile(tp, 0.995) & tp > quantile(tp, 0.005))
  
  ggplot(fin_testt,aes(x=factor(Age),y=tp)) +
    theme_minimal(base_size = 20) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    geom_violin(fill = 'gold3', scale = 'width') +
    scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
    geom_boxplot(width=0.1, fill = 'grey80') +
    geom_point(aes(x = factor(Age), y = Age), size = 3, color = 'darkblue')
  ggsave(paste(path_summary,version,'_ViolinPlot_multivariate_TestPrediction.png', sep = ''), height = 9, width = 16, units = 'in', dpi = 600)
  
  cor_pred <- cor(fin_test$Age, fin_test$tp, method = 'pearson')
  rmse <- rmse(fin_test$Age, fin_test$tp)
  print(paste('Pearson\'s correlation: ', cor_pred, ', Root Mean Square Error: ', rmse,sep=''))

}

Generate_Transfer_Data <- function(percent_split, split_type, pbmc.data, df5, genes_used, prev_donor_list = c(1:100)) {
  set.seed(42)
  
  # index the dataset:
  df5$ind <- c(1:nrow(df5))
  colnames(pbmc.data) <- df5$ind
  
  if (split_type == 'donor' & length(prev_donor_list) < 100) {
    
  	# choose donors for transfer and test data:
  	donors_young <- unique(df5[df5$AGE < 60,]$BC_ID)
  	donors_old <- unique(df5[df5$AGE >= 60,]$BC_ID)
  	
  	# If previous donor list > 1, use the list and add 2 new donors
  	if (length(prev_donor_list) > 0) {
  	  train_set <- c(prev_donor_list, sample(donors_young[!(donors_young %in% prev_donor_list)],1), sample(donors_old[!(donors_old %in% prev_donor_list)],2))
  	  
  	} else {
  	number_young <- ceiling(length(donors_young)*percent_split)
  	number_old <- ceiling(length(donors_old)*percent_split)
    train_set <- c(sample(donors_old, number_old), sample(donors_young, number_young))
  	
    }
  	meta_train <- df5 %>% filter(BC_ID %in% train_set)
  	index_train <- meta_train$ind
  	
  } else if (split_type == 'cell') {
  	train_set <- df5 %>% group_by(donor) %>% sample_frac(percent_split, replace = F)
  	index_train <- train_set$ind
  }
  #Filter transfer set donors:
  index_train <- sample(index_train,length(index_train)) 
  
  df_0 <- df5
  index_test <- df_0[!(df_0$ind %in% index_train),]$ind
  index_test <- sample(index_test,length(index_test))
  
  df5 <- df5 %>% mutate(Sex_OHE = ifelse(SEX == 'M',1,0))
  
  # Generate Training/Test Data/Labels:
  
  # One-hot encode the sex variable
  labels <- matrix(c(df5$AGE,df5$Sex_OHE),ncol=2)
  colnames(labels) <- c('Age','Sex'); rownames(labels) <- df5$ind
  # Data Cleaning:
  # Use cells as observations
  transform_data <- t(as.matrix(pbmc.data))
  transform_data <- transform_data[,head(genes_used,-1)]
  transform_data <- cbind(transform_data,Sex=as.numeric(labels[,'Sex']))
  # Training
  train_data <- transform_data[index_train,]
  train_labels <- labels[index_train,1]
  # Test
  test_data <- transform_data[index_test,]
  test_labels <- labels[index_test,1]
  
  return(list(train_data, train_labels, test_data, test_labels, train_set))
}


Transfer_DL_Model <- function(train_data, transfer_data, transfer_labels, test_data, test_labels,
	trainable = T, l1_node = 30, l2_node = 15, l1_do = 0.25, l2_do = 0.15, learning.rate = 0.0005,
	epochs = 10, patience = 10, times_run = 1, 
	model_weights = '/gpfs/gsfs5/users/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/data/interim/machine_learning/_6.0.1_/_0.2_30_0.25_15_0.15_/_7/weights.241-144.46.hdf5',
	name = '') {
	
  # train_data <- training_data
  # transfer_data <- readRDS(paste('./data/interim/machine_learning/transfer/',version,'_',percent,'_transfer_data.Rds',sep=''))
  # transfer_labels <- readRDS(paste('./data/interim/machine_learning/transfer/',version,'_',percent,'_transfer_labels.Rds',sep=''))
  # test_data <- readRDS(paste('./data/interim/machine_learning/transfer/',version,'_',percent,'_testing_data.Rds',sep=''))
  # test_labels <- readRDS(paste('./data/interim/machine_learning/transfer/',version,'_',percent,'_testing_labels.Rds',sep=''))
  
  require(keras)
  require(tensorflow)
  require(tidyverse)
  require(MASS)
  require(Matrix)
  require(gtools)
  require(MLmetrics)
  
	# Standardize Data:

	train_data <- scale(train_data)
	col_means_train <- attr(train_data, "scaled:center") 
	col_stddevs_train <- attr(train_data, "scaled:scale")
	transfer_data <- scale(transfer_data, center = col_means_train, scale = col_stddevs_train)
	test_data <- scale(test_data, center = col_means_train, scale = col_stddevs_train)

	# Setup Transfer:

	set.seed(1234)

	### Setup Network Topology

	build_model <- function(trainable) {
	  
	  model <- keras_model_sequential() %>%
	    layer_dense(units = l1_node, 
	                trainable = trainable,
	                kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01),
	                activation = "relu", 
	                input_shape = dim(train_data)[2]) %>%
	    layer_dropout(l1_do) %>%
	    layer_dense(units = l2_node,
	    		#trainable = trainable,
	                activation = "relu") %>%
	    layer_dropout(l2_do) %>%
	    layer_dense(units = 1)
	  
	  model %>% compile(
	    loss = "mse",
	    optimizer = optimizer_rmsprop(lr = learning.rate),
	    metrics = list("mean_absolute_error")
	  )
	  
	  model
	}

	model <- build_model(T)
	model %>% summary(model)

	### Evaluate Model

	# Display training progress by printing a single dot for each completed epoch.
	print_dot_callback <- callback_lambda(
	  on_epoch_end = function(epoch, logs) {
	    if (epoch %% 80 == 0) cat("\n")
	    cat(".")
	  }
	)    

	model1 <- build_model(trainable)
	summary(model1)
	model1 %>% load_model_weights_hdf5(model_weights)

	early_stop <- callback_early_stopping(monitor = "val_loss", patience = patience)
	dir.create(paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1',sep=''))
	dir.create(paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/weights',sep=''))
	weights_control <- paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/weights/control',sep='')
	weights_transfer <- paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/weights/transfer',sep='')
	dir.create(weights_control)
	dir.create(weights_transfer)

	# Control Study:
	checkpoint_dir <- weights_control
	filepath <- file.path(checkpoint_dir, paste("weights.{epoch:02d}-{val_loss:.2f}", '.hdf5', sep = '_'))
	check <- callback_model_checkpoint(filepath = filepath, monitor='val_loss',save_weights_only=T)

	history <- model %>% fit(
	  transfer_data,
	  transfer_labels,
	  epochs = epochs,
	  validation_split = 0.1,
	  #validation_data = list(test_data2,test_labels2),
	  verbose = T,
	  callbacks = list(print_dot_callback, early_stop, check)
	)

	plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
	  theme_minimal() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	  coord_cartesian(xlim = c(0, 10), ylim = c(0, 50))
	ggsave(paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/control_learning_curve.png', sep = ''), height = 5, width = 4, units = 'in')

	# Transfer Study:
	checkpoint_dir <- weights_transfer
	filepath <- file.path(checkpoint_dir, paste("weights.{epoch:02d}-{val_loss:.2f}", '.hdf5', sep = '_'))
	check <- callback_model_checkpoint(filepath = filepath, monitor='val_loss',save_weights_only=T)

	history <- model1 %>% fit(
	  transfer_data,
	  transfer_labels,
	  epochs = epochs,
	  validation_split = 0.1,
	  #validation_data = list(test_data2,test_labels2),
	  verbose = T,
	  callbacks = list(print_dot_callback, early_stop, check)
	)

	print('plotting transfer study...')
	plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
	  theme_minimal() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	  coord_cartesian(xlim = c(0, 10), ylim = c(0, 50))
	ggsave(paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/transfer_learning_curve.png', sep = ''), height = 5, width = 4, units = 'in')

	# Visualize on Test set 1:

	print('visualizing on test set...')
	c(loss, mae) %<-% (model1 %>% evaluate(test_data, test_labels, verbose = 0))
	paste0("Mean absolute error on test set: ", mae)

	test_predictions <- model1 %>% predict(test_data)
	fin_test <- tibble(Age = test_labels, prediction = test_predictions[,1])

	#fin_test_avg <- fin_test %>%
	 # group_by(Age) %>%
	  #summarise(avg = mean(tp), sd = sd(tp))

	ggplot(data = fin_test, aes(x = factor(Age), y = prediction)) +
	  theme_minimal() + theme_bw(base_size = 20) +
	  geom_violin(fill='forestgreen') +
	  geom_boxplot(width=0.1, fill = 'grey80') +
	  scale_y_continuous(limits = c(0,100)) +
	  geom_point(aes(x = factor(Age), y = Age), color = 'orange') +
	  theme(aspect.ratio = 5/3)
	  #geom_errorbar(aes(x = Age, ymin = avg - sd, ymax = avg + sd)) +
	  #scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100),limits = c(0,100)) +
	  #scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100),limits = c(0,100))

	ggsave(paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/transfer_prediction.png', sep = ''), height = 5, width = 5, units = 'in', dpi = 600)

  print('generating stats...')
	transfer.stats <- fin_test %>% mutate(error = Age - prediction) %>%
	  summarize(rmse = sqrt(mean(error^2)), corre = cor(Age,prediction,method = 'pearson'))

	write.csv(fin_test, paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/transfer_prediction.csv', sep = ''))
write.csv(transfer.stats, paste('./data/interim/machine_learning/transfer/',version,'_',name,'_transfer_1/transfer_stats.csv', sep = ''))
}

custom_predict_keras <- function(model, newdata)  {

	require(keras)

  newdata_keras <- newdata
  res <- predict(model, newdata_keras)
  return(as.numeric(res[,1]))
}


parse_single_cell_prediction <- function(data_mat, data_anno, training_data) {

	# matrix - the gene x cell matrix of counts data
	# df - the metadata for the matrix (includes 'barcode', 'sample', and 'sex' columns)
	# 		if no 'sex' column included, it will be inferred based on expression of XIST and RPS4Y1
	# training_data - original data that model was built on

	require(tidyverse)
	require(Seurat)

	# Record the positions with 0 standard deviation:
	sds <- apply(training_data, 2, sd)
	cd8_genes <- colnames(training_data)[sds != 0]
	cd8_genes <- cd8_genes[-21093]

	# Modify names to correct possible changes to cell barcode format:
	colnames(data_mat) <- paste(substr(colnames(data_mat),1,16), '-', substr(colnames(data_mat),18,18), sep = '')

	# 
	pbmc <- CreateSeuratObject(counts = data_mat)
	pbmc_genes <- rownames(pbmc@assays$RNA@counts)

	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	plot1 <- FeatureScatter(pbmc, pt.size = 0.001, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(pbmc, pt.size = 0.001, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

	pbmc_trim <- subset(pbmc, subset = nFeature_RNA > 300); print('remove Genes < 500: '); pbmc_trim
	pbmc_trim <- subset(pbmc_trim, subset = nFeature_RNA < 2300); print('remove Genes > 2300: '); pbmc_trim
	pbmc_trim <- subset(pbmc_trim, subset = percent.mt < 10); print('remove MT > 10: '); pbmc_trim
	pbmc_trim <- subset(pbmc_trim, subset = percent.mt > 1); print('remove MT < 1: '); pbmc_trim
	pbmc_trim <- subset(pbmc_trim, subset = nCount_RNA > 500); print('remove RNA < 500: '); pbmc_trim
	FeatureScatter(pbmc_trim, pt.size = 0.005, feature1 = "nCount_RNA", feature2 = "percent.mt")
	FeatureScatter(pbmc_trim, pt.size = 0.005, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

	pbmc_trim <- NormalizeData(pbmc_trim, normalization.method = "LogNormalize", scale.factor = 10000)






}