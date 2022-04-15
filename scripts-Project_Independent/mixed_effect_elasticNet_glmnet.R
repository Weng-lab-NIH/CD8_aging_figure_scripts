#############################################################################################################
#													#####
##				Ahlem HAJJEM, François BELLAVANCE, and Denis LAROCQUE		 	 ####
###				Department of Management Sciences HEC Montréal			  	  ###
####				HEC Montréal, 3000, chemin de la Côte-Sainte-Catherine,	   	 	   ##
#####				Montréal, QC, Canada H3T 2A7

###	Revised: 27st september 2017	###
#############################################################################################################

###########################################
###	Description of MERF arguments	###
###########################################

#xnam			#A charcter vector of p columns, corresponding to the names of the p fixed effect covariates (as they appear in the learning dataset).

#MERF.lDB 		#The learning dataset: a dataframe of N rows (i.e. N level-one observations) and (p+2) columns, 
#where the column "cluster.id" corresponds to the variable uniquely identifying the n clusters (or level-two obervations), 
#the column "Y" corresponds to the name of the continuous response variable,
#and the other p columns correspond to the p fixed effects covariates.	

#ni			#A vector of n columns, where ni[i] corresponds to the size of cluster i, 
#for i = 1, ..., n (ATTENTION: should keep the same order as in MERF.lDB).

#Zi			#A list of n matrices Zi[[i]] of dimension (ni[i] X q), where  
#Zi[[i]] corresponds to the q random effects covariates values of the ni[i] observations nested within cluster i, for i= 1, ..., n. 
#Note that q=1 when there is only a random intercept, i.e. no random effect covariates.

#Yi			#A list of n vectors Yi[i] of ni[i] rows, where 
#Yi[i] corresponds to the response values of the ni[i] observations nested within cluster i, for i= 1, ..., n.

#ntree			#ntree argument of randomForest function, i.e., number of trees to grow (e.g. 300).

#mtry			#mtry argument of randomForest function, i.e., number of variables randomly sampled as candidates at each split (e.g. 3).

#nodesize		#nodesize argument of randomForest function, i.e., minimum size of terminal nodes (e.g.5).

#sigmasqzero = NULL	#Starting value of s2, where the covariance matrix of the errors Ri = s2Ini , for i = 1, ..., n. 
#Default: if( is.null(sigmasqzero) ) sigmasqzero <- 1.

#Dzero = NULL		#Starting value of the covariance matrix of the random effects. 
#Default:
#if( is.null(Dzero) ){
#	Dzero <- diag(0.01,nrow=q,ncol=q)
#}

#bizero = NULL		#Starting values of the random effects: a list of n matrices of q × 1 unknown vector of random effects. 
#Default:
#if( is.null(bizero) ){
#bizero <- list(); length(bizero)<- n
#	for(i in 1:n) bizero[[i]] <- matrix(0,nrow=q,ncol=1)}
#}

#F.niter		#The number of iterations forced to avoid early stopping (e.g. 100).

#max.niter		#Maximum number of iterations, in addition to the F.niter forced iterations (e.g. 300).

#smallest.Jump.allowed 	#A given small value (e.g. 1e-4). 

#verbose		#Logical. Should R report extra information on progress?



####################################	BEGIN THE SCRIPT		#####################################

#setwd("/home/...")
#sink("Routput.doc")
#set.seed(321)

########################################################################################################################
###############################                   ### MERF function ###  	         ##############################
##################		fits a mixed effects random forest of regression trees model		###############
########################################################################################################################


library(glmnet)
library(glmnetUtils)

MEEN_glmnet  <- function(
  
  xnam
  
  ,MERF.lDB 	
  
  ,ni
  ,Zi
  ,Yi	
  
  ,model_type = 'enet'
  ,alpha = 0.5 
  ,name_model = '..'
  ,n_fold = 10
  
  
  ,sigmasqzero = NULL
  ,Dzero = NULL
  ,bizero = NULL
  
  ,F.niter
  ,max.niter
  ,smallest.Jump.allowed
  
  ,verbose = TRUE
  ,ncore = 1
  
){
  
  ####################
  ####	STEP 0	####
  #################### 
  
  #Memory Allocation and initialization:
  
  #Parameters values 
  n <- length(ni)
  N <- sum(ni)
  q <- dim(Zi[[1]])[2] 	# q=1 in random intercept case
  
  
  #Initial values of sigmasqzero, Dzero, and bizero
  if( is.null(sigmasqzero) ) sigmasqzero <- 1
  else sigmasqzero <- sigmasqzero
  
  if( is.null(Dzero) ){
    Dzero <- diag(0.01 ,nrow=q,ncol=q)
  }
  else Dzero <- Dzero
  
  if( is.null(bizero) ){
    bizero <- list(); length(bizero)<- n
    for(i in 1:n) bizero[[i]] <- matrix(0,nrow=q,ncol=1)
  }	else { bizero <- bizero }
  
  #iter number
  r <- 1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #transformed outcome, star.Yi[[r]][[i]], initialized with the original values
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bizero[[i]]
  }
  
  #one STD random forest
  ######################
  
  MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  
  fit.rf.formula <-  as.formula(paste("star.Yi ~ ", paste(xnam, collapse= "+")))
  
  fit.enet <- cv.glmnet(formula = fit.rf.formula,
                        data = MERF.lDB, 
                       alpha = 0.5, 
                       n_folds = n_fold,
                       keep = T)
  
  #fixed part
  #as vector
  currInd <- match(fit.enet$lambda.min,fit.enet$glmnet.fit$lambda)
  
  MERF.lDB$f.pred <- as.numeric(fit.enet$fit.preval[,currInd]) # THese are the out of bag predictions!!!!
  # predict(fit.ranger, data = MERF.lDB)			#!!! use the out-of-bag predictions - these are likely not out of bag
  
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, MERF.lDB$cluster.id) 	
  for(i in 1:n) {
    fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
  }
  
  #random	part
  ############
  #random effects parameters in list format
  bi <- list(list()) ; length(bi) <- r
  for(i in 1:n)bi[[r]][[i]] <- bizero[[i]] 	
  #print("bizero");print(bizero)
  rm(bizero) ; gc(verbose=FALSE)
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
  
  sigma.sq <- vector(mode="numeric") ;length(sigma.sq) <- r
  sigma.sq[r] <- sigmasqzero
  #print("sigmasqzero") ;print(sigmasqzero)			
  rm(sigmasqzero) ; gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  D <- list() ;length(D) <- r
  D[[r]] <- Dzero		#!!!Dzero <- diag(x=0.01, nrow=q, ncol = q)
  #print("Dzero") ;print(Dzero)	
  rm(Dzero) ; gc(verbose=FALSE)
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  Vi <- list() 
  inv.Vi <- list(list()) ; length(inv.Vi) <- r
  
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]]) + sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                           -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                 , ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  
  Vi <- list(NULL) 
  #inv.Vi[[r-1]] <- list(NULL) #not to run at step 0
  
  #the generalized log-likelihood (GLL) 
  GLL <- vector(mode="numeric") ; length(GLL) <- r
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
  + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
  + log(abs(D[[r]]))
  + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #convergence criterion
  Jump <- rep(NA,r) 		#at this first iteration Jump = NA
  convergence.iter<- rep(NA,r) 	#at this first convergence.iter = NA
  
  ####################
  ####	STEP 1	####        
  ####################
  
  #update iteration number r
  r <- r+1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #update the length of the different lists
  
  length(sigma.sq) <- r
  length(D) <- r
  length(inv.Vi) <- r
  length(bi) <- r
  length(GLL) <- r
  
  length(Jump) <- r
  length(convergence.iter) <- r
  
  #update the transformed outcome, star.Yi
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
  }	
  
  #one STD random forest
  ######################
  
  MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  
  fit.enet <- cv.glmnet(formula = fit.rf.formula,
                        data = MERF.lDB, 
                        alpha = 0.5, 
                        n_folds = n_fold,
                        keep = T)
  
  #fixed part
  #as vector
  currInd <- match(fit.enet$lambda.min,fit.enet$glmnet.fit$lambda)
  
  MERF.lDB$f.pred <- as.numeric(fit.enet$fit.preval[,currInd]) # THese are the out of bag predictions!!!!
 

  
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, MERF.lDB$cluster.id) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
  
  #random	part
  ############
  for(i in 1:n)
    bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
  
  bi[r-1] <- list(NULL)		 
  
  ####################
  ####	STEP 2	####         
  ####################
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
  
  
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i] <- crossprod(epsili[[i]]) + 
    sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
  sigma.sq[r] <- (1/N)*(sum(term))
  rm(term) ;gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  term <- list()
  term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
    )
  for(i in 2:n) 
    term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
    )
  term <- term[[n]]
  D[[r]] <- (1/n)*term
  rm(term) ;gc(verbose=FALSE)	 
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  inv.Vi[[r]] <-list()
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                           -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                 , ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  Vi <- list(NULL) 
  inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
  
  
  #the generalized log-likelihood (GLL) 
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
  + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
  + log(abs(D[[r]]))
  + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #update the value of the Jump in GLL
  Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
  
  if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
    convergence.iter[r] <- r
    if (verbose) message("Converg. at iter no: ", r)
  } 
  
  
  ####################
  ####	STEP 3	####         
  ####################
  
  
  ###################################################
  #Iterating "F.niter" times to avoid early stopping#
  ###################################################
  
  for(I in 1:F.niter){#repeat step 1 and 2
    
    ####################
    ####	STEP 1	####        
    ####################
    
    #update iteration number r
    r <- r+1
    if (verbose) 	
      message("MERF iter no: ", r)
    
    #update the length of the different lists
    
    length(sigma.sq) <- r
    length(D) <- r
    length(inv.Vi) <- r
    length(bi) <- r
    length(GLL) <- r
    
    length(Jump) <- r
    length(convergence.iter) <- r
    
    #update the transformed outcome, star.Yi
    star.Yi <- list() 
    for(i in 1:n){
      star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
    }
    
    #one STD random forest
    ######################
    
    MERF.lDB$star.Yi <- unlist(star.Yi) 
    rm(star.Yi) ;  gc(verbose=FALSE)
    
    fit.enet <- cv.glmnet(formula = fit.rf.formula,
                          data = MERF.lDB, 
                          alpha = 0.5, 
                          n_folds = n_fold,
                          keep = T)
    
    #fixed part
    #as vector
    currInd <- match(fit.enet$lambda.min,fit.enet$glmnet.fit$lambda)
    
    MERF.lDB$f.pred <- as.numeric(fit.enet$fit.preval[,currInd]) # THese are the out of bag predictions!!!!
    #in matrix format
    fixed.pred <- list()	
    fixed.pred <- split(MERF.lDB, MERF.lDB$cluster.id) 	
    for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
    
    #random	part
    ############
    for(i in 1:n)
      bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
    bi[r-1] <- list(NULL)		 
    
    
    ####################
    ####	STEP 2	####        
    ####################
    
    #level-1 variance component
    #residuals
    epsili <- list()
    for(i in 1:n)
      epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
    
    
    term <- vector(mode="numeric",length=n)
    for(i in 1:n)
      term[i] <- crossprod(epsili[[i]]) + 
      sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
    sigma.sq[r] <- (1/N)*(sum(term))
    rm(term) ;gc(verbose=FALSE)
    #message("sigmasq of current micro iter", sigma.sq[r] )
    
    #level-2 variance component
    term <- list()
    term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
      (	D[[r-1]] - 
          D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
      )
    for(i in 2:n) 
      term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
      (	D[[r-1]] - 
          D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
      )
    term <- term[[n]]
    D[[r]] <- (1/n)*term
    rm(term) ;gc(verbose=FALSE)	 
    #message("D of current micro iter: ", D[[r]] )
    
    #level-1 and level-2 variance components (or typical or total variance)
    inv.Vi[[r]] <-list()
    for(i in 1:n){
      Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
      if(q==1)
        inv.Vi[[r]][[i]] <- 
          (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                             -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                   , ncol=ni[i], nrow=ni[i]) )
      else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
    }
    Vi <- list(NULL) 
    inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
    
    
    #the generalized log-likelihood (GLL) 
    term <- vector(mode="numeric",length=n)
    for(i in 1:n)
      term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
    + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
    + log(abs(D[[r]]))
    + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
    GLL[r] <- sum(term)
    rm(term)
    gc(verbose=FALSE)
    
    #update the value of the Jump in GLL
    
    Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
    
    if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
      convergence.iter[r] <- r
      if (verbose) message("Converg. at iter no: ", r)
    } 
    
  }
  #end for (I in 1: F.niter)
  ###########################
  
  ######################################
  #Iterating "max.niter" times at most #
  ######################################
  
  while( r < (2 + F.niter + max.niter) ){
    
    if(Jump[r] > smallest.Jump.allowed){ #repeat step 1 and 2
      
      ####################
      ####	STEP 1	####        
      ####################
      
      #update iteration number r
      r <- r+1
      if (verbose) 	
        message("MERF iter no: ", r)
      
      #update the length of the different lists
      
      length(sigma.sq) <- r
      length(D) <- r
      length(inv.Vi) <- r
      length(bi) <- r
      length(GLL) <- r
      
      length(Jump) <- r
      length(convergence.iter) <- r
      
      #update the transformed outcome, star.Yi
      star.Yi <- list() 
      for(i in 1:n){
        star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
      }
      
      #one STD random forest
      ######################
      
      MERF.lDB$star.Yi <- unlist(star.Yi) 
      rm(star.Yi) ;  gc(verbose=FALSE)
      
      fit.enet <- cv.glmnet(formula = fit.rf.formula,
                            data = MERF.lDB, 
                            alpha = 0.5, 
                            n_folds = n_fold,
                            keep = T)
      
      #fixed part
      #as vector
      currInd <- match(fit.enet$lambda.min,fit.enet$glmnet.fit$lambda)
      
      MERF.lDB$f.pred <- as.numeric(fit.enet$fit.preval[,currInd]) # THese are the out of bag predictions!!!!
      
      #in matrix format
      fixed.pred <- list()	
      fixed.pred <- split(MERF.lDB, MERF.lDB$cluster.id) 	
      for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
      
      #random	part
      ############
      for(i in 1:n)
        bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
      bi[r-1] <- list(NULL)		
      
      
      ####################
      ####	STEP 2	####        
      ####################
      
      #level-1 variance component
      #residuals
      epsili <- list()
      for(i in 1:n)
        epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
      
      
      term <- vector(mode="numeric",length=n)
      for(i in 1:n)
        term[i] <- crossprod(epsili[[i]]) + 
        sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
      sigma.sq[r] <- (1/N)*(sum(term))
      rm(term) ;gc(verbose=FALSE)
      #message("sigmasq of current micro iter", sigma.sq[r] )
      
      #level-2 variance component
      term <- list()
      term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
        (	D[[r-1]] - 
            D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
        )
      for(i in 2:n) 
        term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
        (	D[[r-1]] - 
            D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
        )
      term <- term[[n]]
      D[[r]] <- (1/n)*term
      rm(term) ;gc(verbose=FALSE)	 
      #message("D of current micro iter: ", D[[r]] )
      
      #level-1 and level-2 variance components (or typical or total variance)
      inv.Vi[[r]] <-list()
      for(i in 1:n){
        Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
        if(q==1)
          inv.Vi[[r]][[i]] <- 
            (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                               -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                     , ncol=ni[i], nrow=ni[i]) )
        else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
      }
      Vi <- list(NULL) 
      inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
      
      
      #the generalized log-likelihood (GLL) 
      term <- vector(mode="numeric",length=n)
      for(i in 1:n)
        term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
      + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
      + log(abs(D[[r]]))
      + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
      GLL[r] <- sum(term)
      rm(term)
      gc(verbose=FALSE)
      
      #update the value of the Jump in GLL
      
      Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
      
      if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
        convergence.iter[r] <- r
        if (verbose) message("Converg. at iter no: ", r)
      } 
      
      
    }
    #end if(Jump[r] > smallest.Jump.allowed) and STOP repeating step 1 and 2
    
    else break 
    #end while( r < (2 + F.niter + max.niter) )
    
    
    ############################
  }###	END OF STEP 3	####        
  ############################
  
  
  #output to be returned (MERF model is the one at the last iteration)
  ###################### 
  
  
  output <- list(
    Jump[r]
    ,GLL
    ,convergence.iter
    ,fit.enet 
    ,bi[[r]]
    ,sigma.sq[r]
    ,D	#,D[[r]]
  )
  
  names(output) <- c(
    "Jump[r]"
    ,"GLL"
    ,"convergence.iter"
    ,"fit.enet"
    ,"bi[[r]]"
    ,"sigma.sq[r]"
    ,"D"	#,"D[[r]]"
  )
  
  #clean memory
  #############
  rm(
    xnam
    
    ,MERF.lDB	
    
    ,ni,n,N
    ,Zi,q
    ,Yi	
    
    ,fit.rf.formula 
    ,fit.enet
    ,fixed.pred ,epsili
    
    ,sigma.sq,D,bi,Vi,inv.Vi
    
    ,F.niter ,max.niter
    ,smallest.Jump.allowed ,GLL ,Jump ,convergence.iter
    ,r,i,I
    
    ,verbose 		
  )
  gc(verbose=FALSE)
  
  
  #return
  #######
  output
  
}#end of MERF function
############################################################################################################
############################################################################################################
############################################################################################################
### end of the script	###
#sink()