Regress_Variable <- function(matrix, metadata, variable, n_cores) {

	require(doParallel)
	require(parallel)
	require(foreach)
  require(tidyverse)

	cl <- makeCluster(n_cores)
	registerDoParallel(cl, cores = n_cores)

	df <- metadata
	data <- matrix
	var <- variable
	
	part <- df %>% filter(get(cluster_level) == l & BC_ID != 1 & BC_ID != 2)
	cell_list <- as.character(part$bc)
	p_clust <- data[,cell_list]
	t_clust <- t(as.matrix(p_clust))
	
	Age <- part[,colnames(part) == var]
	Sex <- part$Sex
	#Age[Age == 0.1] <- 0
	Age[Age == 83.1] <- 83

	listcl <- foreach (i=1:length(colnames(t_clust)), .combine='rbind') %dopar% {
      		#i = 1
      		exp <- t_clust[,i]
      		mdf <- data.frame(exp = exp, Age = Age, Sex = Sex)

	model <- lm(exp~Age, data = mdf) # This is the actual linear model generated, change this to preferred model if necessary
	model <- glm(formula = exp ~ Age, family = binomial(link = "logit"), data = mdf)
      		#summary(model)
      		age_p <- summary(model)$coefficients[,4]["Age"]
      		co <- summary(model)$coefficients[,1]["Age"]
      		sex_p <- 0 
      		co_sex <- 0 
      		intercept <- summary(model)$coefficients[,1]["(Intercept)"]
      		r2 <- summary(model)$r.squared
      		adj_r2 <- summary(model)$adj.r.squared

    df_log2_ind <- data.frame(cluster = l, gene = as.character(colnames(t_clust)[i]), age_p = age_p, age_p_adj = 0, sex_p = sex_p, sex_p_adj = 0, co = co, intercept = intercept, co_sex = co_sex, r2 = r2, adj_r2 = adj_r2,
                                age_p_3_2 = age_p_3_2, age_p_adj_3_2 = 0, age_p_3_1 = age_p_3_1, age_p_adj_3_1 = 0, sex_p_3 = sex_p_3, sex_p_adj_3 = 0, co_3_2 = co_3_2, co_3_1 = co_3_1, intercept_3 = intercept_3, co_sex_3 = co_sex_3, r2_3 = r2_3, adj_r2_3 = adj_r2_3,
                                age_p_4_3 = age_p_4_3, age_p_adj_4_3 = 0, age_p_4_2 = age_p_4_2, age_p_adj_4_2 = 0, age_p_4_1 = age_p_4_1, age_p_adj_4_1 = 0, sex_p_4 = sex_p_4, sex_p_adj_4 = 0, co_4_3 = co_4_3, co_4_2 = co_4_2, co_4_1 = co_4_1, intercept_4 = intercept_4, co_sex_4 = co_sex_4, r2_4 = r2_4, adj_r2_4 = adj_r2_4)
     }

}



Regress_Variable_From_Clusters <- function(matrix, metadata, variable, cluster_level, n_cores) {
	print('This function regresses out a variable from the dataset and returns summaries for each gene')
  # 'Age' is used as a dummy variable for the regression code
  
	require(doParallel)
	require(parallel)
	require(foreach)
  require(tidyverse)

	cl <- makeCluster(n_cores)
	registerDoParallel(cl, cores = n_cores)

  # Assign Variables:
	df <- metadata
	data <- matrix
	var <- variable
	
  	final <- data.frame(cluster = 0, gene = 0, age_p = 0, age_p_adj = 0, sex_p = 0, sex_p_adj = 0, co = 0, intercept = 0, co_sex = 0, r2 = 0, adj_r2 = 0) # ,
  	#age_p_3_2 = 0, age_p_adj_3_2 = 0, age_p_3_1 = 0, age_p_adj_3_1 = 0, sex_p_3 = 0, sex_p_adj_3 = 0, co_3_2 = 0, co_3_1 = 0, intercept_3 = 0, co_sex_3 = 0, r2_3 = 0, adj_r2_3 = 0,
  	#age_p_4_3 = 0, age_p_adj_4_3 = 0, age_p_4_2 = 0, age_p_adj_4_2 = 0, age_p_4_1 = 0, age_p_adj_4_1 = 0, sex_p_4 = 0, sex_p_adj_4 = 0, co_4_3 = 0, co_4_2 = 0, co_4_1 = 0, intercept_4 = 0, co_sex_4 = 0, r2_4 = 0, adj_r2_4 = 0)
  
  	clust <- unique(df[,colnames(df) == cluster_level])
  	clust <- clust[[1]]
  	clust <- clust[is.na(clust) == F]
  
  # Perform Regression for Each Cluster:
  	for (l in clust) {
    	#l <- 0
    	part <- df %>% filter(get(cluster_level) == l) # & BC_ID != 1 & BC_ID != 2)
    	cell_list <- as.character(part$bc)
    	p_clust <- data[,cell_list]
    	t_clust <- t(as.matrix(p_clust))
    
    	Age <- part[,colnames(part) == var]
    	Age <- Age[[1]]
    	Sex <- part$Sex
    	#Age[Age == 0.1] <- 0
    	Age[Age == 83.1] <- 83
    
    	listcl <- foreach (i=1:length(colnames(t_clust)), .combine='rbind') %dopar% {
      		#i = 1
      		exp <- t_clust[,i]
      		mdf <- data.frame(exp = exp, Age = Age, Sex = Sex)
      
      		model <- lm(exp~Age+Sex, mdf) # This is the actual linear model generated, change this to preferred model if necessary
      		#summary(model)
      		age_p <- summary(model)$coefficients[,4]["Age"]
      		co <- summary(model)$coefficients[,1]["Age"]
      		sex_p <- summary(model)$coefficients[,4]["SexM"]
      		co_sex <- summary(model)$coefficients[,1]["SexM"]
      		intercept <- summary(model)$coefficients[,1]["(Intercept)"]
      		r2 <- summary(model)$r.squared
      		adj_r2 <- summary(model)$adj.r.squared
      
      		# model <- lm(exp~poly(Age, 2, raw=T),mdf)
      		# age_p_3_2 <- summary(model)$coefficients[,4]["poly(Age, 2, raw = T)2"]
      		# co_3_2 <- summary(model)$coefficients[,1]["poly(Age, 2, raw = T)2"]
      		# age_p_3_1 <- summary(model)$coefficients[,4]["poly(Age, 2, raw = T)1"]
      		# co_3_1 <- summary(model)$coefficients[,1]["poly(Age, 2, raw = T)1"]
      		# sex_p_3 <- 0 #summary(model)$coefficients[,4]["SexM"]
      		# co_sex_3 <- 0 #summary(model)$coefficients[,1]["SexM"]
      		# intercept_3 <- summary(model)$coefficients[,1]["(Intercept)"]
      		# r2_3 <- summary(model)$r.squared
      		# adj_r2_3 <- summary(model)$adj.r.squared
      		# 
      		# model <- lm(exp~poly(Age, 3, raw=T),mdf)
      		# age_p_4_3 <- summary(model)$coefficients[,4]["poly(Age, 3, raw = T)3"]
      		# co_4_3 <- summary(model)$coefficients[,1]["poly(Age, 3, raw = T)3"]
      		# age_p_4_2 <- summary(model)$coefficients[,4]["poly(Age, 3, raw = T)2"]
      		# co_4_2 <- summary(model)$coefficients[,1]["poly(Age, 3, raw = T)2"]
      		# age_p_4_1 <- summary(model)$coefficients[,4]["poly(Age, 3, raw = T)1"]
      		# co_4_1 <- summary(model)$coefficients[,1]["poly(Age, 3, raw = T)1"]
      		# sex_p_4 <- 0 #summary(model)$coefficients[,4]["SexM"]
      		# co_sex_4 <- 0 #summary(model)$coefficients[,1]["SexM"]
      		# intercept_4 <- summary(model)$coefficients[,1]["(Intercept)"]
      		# r2_4 <- summary(model)$r.squared
      		# adj_r2_4 <- summary(model)$adj.r.squared
      		# 
      
      	df_log2_ind <- data.frame(cluster = l, gene = as.character(colnames(t_clust)[i]), age_p = age_p, age_p_adj = 0, sex_p = sex_p, sex_p_adj = 0, co = co, intercept = intercept, co_sex = co_sex, r2 = r2, adj_r2 = adj_r2) #,
                               # age_p_3_2 = age_p_3_2, age_p_adj_3_2 = 0, age_p_3_1 = age_p_3_1, age_p_adj_3_1 = 0, sex_p_3 = sex_p_3, sex_p_adj_3 = 0, co_3_2 = co_3_2, co_3_1 = co_3_1, intercept_3 = intercept_3, co_sex_3 = co_sex_3, r2_3 = r2_3, adj_r2_3 = adj_r2_3,
                               # age_p_4_3 = age_p_4_3, age_p_adj_4_3 = 0, age_p_4_2 = age_p_4_2, age_p_adj_4_2 = 0, age_p_4_1 = age_p_4_1, age_p_adj_4_1 = 0, sex_p_4 = sex_p_4, sex_p_adj_4 = 0, co_4_3 = co_4_3, co_4_2 = co_4_2, co_4_1 = co_4_1, intercept_4 = intercept_4, co_sex_4 = co_sex_4, r2_4 = r2_4, adj_r2_4 = adj_r2_4)
      
    }
    
    final <- rbind(final, listcl)
  }
  
  #stopCluster(cl)
  
  for (l in clust) {
    final[final$cluster == l,]$age_p_adj <- p.adjust(final[final$cluster == l,]$age_p, method = 'BH')
    final[final$cluster == l,]$sex_p_adj <- p.adjust(final[final$cluster == l,]$sex_p, method = 'BH')
    
  #   final[final$cluster == l,]$age_p_adj_3_2 <- p.adjust(final[final$cluster == l,]$age_p_3_2, method = 'BH')
  #   final[final$cluster == l,]$age_p_adj_3_1 <- p.adjust(final[final$cluster == l,]$age_p_3_1, method = 'BH')
  #   final[final$cluster == l,]$sex_p_adj_3 <- p.adjust(final[final$cluster == l,]$sex_p_3, method = 'BH')
  #   
  #   final[final$cluster == l,]$age_p_adj_4_3 <- p.adjust(final[final$cluster == l,]$age_p_4_3, method = 'BH')
  #   final[final$cluster == l,]$age_p_adj_4_2 <- p.adjust(final[final$cluster == l,]$age_p_4_2, method = 'BH')
  #   final[final$cluster == l,]$age_p_adj_4_1 <- p.adjust(final[final$cluster == l,]$age_p_4_1, method = 'BH')
  #   final[final$cluster == l,]$sex_p_adj_4 <- p.adjust(final[final$cluster == l,]$sex_p_4, method = 'BH')
  }

  final <- final[-1,]
	return(final)
}


Generate_Regression_Report <- function(model, gene, cluster) {
	
	age_p <- summary(model)$coefficients[,4]["Age"]
    co <- summary(model)$coefficients[,1]["Age"]
      
    intercept <- summary(model)$coefficients[,1]["(Intercept)"]
    r2 <- summary(model)$r.squared
    adj_r2 <- summary(model)$adj.r.squared
      
    regression <- data.frame(cluster = cluster, gene = gene, co = co, p = age_p, adj_r2 = adj_r2)
    return(regression)
}

Calculate_ModesOfChange <- function(matrix, metadata, cluster_level, 
                                    threshold, regression, p_value, fold_change, 
                                    regression_level, n_cores) {
	print('clusters: dataframe
		column 1 - \'bc\',
		column 2 - \'cluster\'
		column 3 - \'BC_ID\'
		column 4 - \'Age\'
		column 5 - \'Sex\'
		column 6 - \'cluster\'
	  column 7 - \'subset\'
	  column 8 - \'cell\'')
  
  require(doParallel)
  require(parallel)
  require(foreach)
  require(tidyverse)
  source('/gpfs/gsfs5/users/TCR/__SCRIPTS_Raheel/scripts-Project_Independent/regression_functions.R')
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("Generate_Regression_Report"))
  registerDoParallel(cl, cores = n_cores)
  
  #matrix <- pbmc_counts
  #metadata <- clusters
  #cluster_level <- 'cluster'
  #threshold <- 0.1
  #regression <- regression.cluster
  #p_value <- t_ptotal
  #fold_change <- t_exptotal
  
	clust <- unique(metadata[,colnames(metadata) == cluster_level])
	clust <- clust[[1]]
	clust <- clust[is.na(clust) == F]
	
	dfin <- data.frame()
	dfin2 <- data.frame()
	for (l in clust) {

		cluster.metadata <- metadata %>% filter(get(cluster_level) == l)
		cluster.cells <- as.character(cluster.metadata$bc)
		cluster.top.regression <- regression %>% filter(cluster == l & age_p_adj < p_value & abs(co) > fold_change)
		cluster.top.genes <- cluster.top.regression$gene

		cluster.exp.mat <- matrix[,cluster.cells]
		
		comb <- function(...) {
		  mapply('rbind', ..., SIMPLIFY=FALSE)
		}
		
		listcl <- foreach (i=1:length(cluster.top.genes), .packages='tidyverse', .combine='comb', .multicombine = T) %dopar% {
		  
		  gene <- cluster.top.genes[i]
			expression.vector <- as.vector(cluster.exp.mat[as.character(gene),])
			cluster.metadata$expression <- expression.vector

			if (sum(expression.vector) > 0) {

			 	  ##############################
    			### Test 1: PERCENT CHANGE ###
    			##############################

		 	    # calculate percentage of cells detectable (plot + regression)
    			
    			test1 <- cluster.metadata
    			test1.percents <- cluster.metadata %>%
      			group_by(Age) %>%
      			summarise(norm = sum(expression > threshold)/n())
    			test1.percents$norm <- test1.percents$norm*100
    			model.1 <- lm(norm ~ Age, data = test1.percents)
      		regression.1.summary <- Generate_Regression_Report(model.1, gene, l) %>% mutate(mode = 'percent')
      		
          if (regression_level == 'cell') {
            ############################################
      			### Test 2: EXPRESSION CHANGE-CELL LEVEL ###
      			############################################
        		
        		test2 <- test1 %>% filter(expression > threshold)
      
            model.2 <- lm(expression ~ Age, data = test2)
            regression.2.summary <- Generate_Regression_Report(model.2, gene, l) %>% mutate(mode = 'expression')
            
            all <- test1.percents %>% mutate(cluster = l, gene = gene)
          } else if (regression_level == 'person') {
        		##################################################
     			  ### Test 2: EXPRESSION CHANGE-PERSON LEVEL     ###
      			##################################################
  
      			test1.positive <- test1 %>% filter(expression > threshold)
      			test2 <- test1.positive %>%
        		group_by(Age) %>%
        		summarise(avg = mean(expression))
  
      			# Perform regression to determine if the expression changes in single cells:
      			model.2 <- lm(avg ~ Age, data = test2)
      			regression.2.summary <- Generate_Regression_Report(model.2, gene, l) %>% mutate(mode = 'expression')
      			all <- test1.percents %>% inner_join(test2, by = 'Age') %>% mutate(cluster = l, gene = gene)
          }
      		regression.final <- rbind(regression.1.summary, regression.2.summary)
      		
      	
			}
			
			output <- list(regression.final, all)
			
		}
		
		dfin <- rbind(dfin, listcl[[1]])
		dfin2 <- rbind(dfin2, listcl[[2]])
		
	}

	data_list <- list(dfin, dfin2)
	return(data_list)
}

Score_Modes <- function(score, thresh_p, thresh_e) {
  print(summary(score$mode=='percent'))
  score <- inner_join(score[score$mode == 'percent',], score[score$mode == 'expression',], by = c('cluster','gene'))
  score <- score %>% mutate(exp_percent = (2^(`co.y`*70)-1))
  score <- score %>% mutate(type = ifelse(abs(co.x) >  thresh_p& abs(exp_percent) > thresh_e, 'both',
                                          ifelse(abs(co.x) > thresh_p & abs(exp_percent) < thresh_e, 'percent', 
                                                 ifelse(abs(co.x) < thresh_p & abs(exp_percent) > thresh_e, 'expression', 'none'))))
}