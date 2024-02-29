library(dplyr)
library(tidyr)
d_convergence_threshold =  1e-10

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

fBM_ALg1_inner_loop <- function(dta_s, df_cluster, sX, sY, theta_0, alpha_0, iG, sOLS_formula,
                                time_dummies, group_dummies, max_iter=100,
                                clustered_coefs=FALSE,
                                s_interaction_clustered_coefs=NULL,
                                s_interaction_terms=NULL){
  theta_s <- theta_0
  alpha_s <- alpha_0
  old_param <- 0
  
  iT <- length(unique(dta_s$year))
  n_iter <- 0
  change_param <- 1
  while (((change_param > d_convergence_threshold) & (n_iter < max_iter)) | (n_iter < 3)){
    n_iter <- n_iter + 1
    
    # print(paste0('Change: ', as.character(change_param)))
    
    # 1st step: assignment
    # first part of the error:
    if (clustered_coefs){
      error_1 <- dta_s[sY]-0
      error_g <- list()
      for (g in 1:iG){
        # First the time dummies
        alpha_g <- alpha_s[grep(paste0('G_', g), names(alpha_s))]
        # Make sure they have the correct order
        names(alpha_g) <- gsub('^I\\((time_[0-9]{4}).*$', '\\1' , names(alpha_g))
        error_g[[g]] <- as.matrix(dta_s[time_dummies]) %*% alpha_g[time_dummies]
        # Second the clustered coefs
        theta_g <- theta_s[grep(paste0('G_', g), names(theta_s))]
        names(theta_g) <- gsub('^I\\((.*) \\* G\\_[0-9]\\)$', '\\1' , names(theta_g))
        error_g[[g]] <- error_g[[g]] + as.matrix(dta_s[sX]) %*% theta_g[sX]
      }
    }else{
      error_1 <- dta_s[sY] - (as.matrix(dta_s[sX]) %*% theta_s)
      error_g <- list()
      for (g in 1:iG){
        alpha_g <- alpha_s[(iT*(g-1) + 1):(iT*g)] # MUST BE IN ORDER
        error_g[[g]] <- as.matrix(dta_s[time_dummies]) %*% alpha_g
      }
    }
    
    old_clustering <- df_cluster['cluster']
    clustering_error <- c()
    for (i in 1:nrow(df_cluster)){
      code <- df_cluster$code[i]
      # Take the error:
      error_1_code <- error_1[dta_s$code == code,]
      error_tot_g <- c()
      for (g in 1:iG){
        error_g_code <- error_g[[g]][dta_s$code == code,]
        error_tot_g[g] <- mean((error_1_code - error_g_code)^2)
      }
      df_cluster[i, 'cluster'] <- which.min(error_tot_g) # New assginment
      clustering_error[i] <- min(error_tot_g)
    }
    clustering_change <- sum(old_clustering != df_cluster['cluster'])
    
    # Redo the cluster variable in the main dataset
    for(j in group_dummies){
      iCluster <- as.numeric(substr(j, 3, 3))
      dta_s[j] <- as.numeric(dta_s$code %in% df_cluster$code[df_cluster$cluster == iCluster])
    }
    
    # 2nd step: OLS
    model_1st_step <- lm(sOLS_formula, dta_s)
    if (clustered_coefs){
      theta_s <- coef(model_1st_step)[strsplit(s_interaction_clustered_coefs, ' \\+ ')[[1]][-1]]
      alpha_s <- coef(model_1st_step)[strsplit(s_interaction_terms, ' \\+ ')[[1]][-1]]
    } else{
      theta_s <- coef(model_1st_step)[sX]
      alpha_s <- coef(model_1st_step)[(length(sX)+1):(length(sX) + length(time_dummies)*iG)] #Takes the other ones
    }
    
    new_param <- coef(model_1st_step)
    new_param[is.na(new_param)] <- 0
    change_param <- mean((old_param - new_param)^2)
    old_param <- new_param
  }
  return(list(model_1st_step, df_cluster, clustering_error, theta_s, alpha_s))
}

fBM_ALg1 <- function(dta_s, df_cluster, sX, sY, iG){
  iT <- length(unique(dta_s$year))
  # iN <- length(unique(dta_s$code))
  
  time_dummies <- paste0('time_', sort(unique(dta_s$year)))
  group_dummies <- paste0('G_', 1:iG)
  # Time dummies are equal to 1 if the year is equal to their year
  for(i in c(time_dummies)){
    sYear <- substr(i, 6, 9)
    dta_s[i] <- as.numeric(dta_s$year == as.numeric(sYear))
  }
  for(i in group_dummies){
    iCluster <- as.numeric(substr(i, 3, 3))
    dta_s[i] <- as.numeric(dta_s$code %in% df_cluster$code[df_cluster$cluster == iCluster])
  }
  
  # Inital OLS step
  s_interaction_terms <- ''
  for (i in group_dummies){
    s_interaction_terms <- paste0(paste0(s_interaction_terms, ' + '),
                                  paste0('I(', paste(time_dummies, i, sep='*', collapse=') + I('), ')'))
  }
  sOLS_formula <- paste0(sY, ' ~ ', paste0(sX, collapse=' + 0 + '), ' + ', s_interaction_terms)
  model_1st_step <- lm(sOLS_formula, dta_s)
  
  old_param <- coef(model_1st_step)
  old_param[is.na(old_param)] <- 0
  change_param <- 1
  
  theta_0 <- coef(model_1st_step)[sX]
  alpha_0 <- coef(model_1st_step)[(length(sX)+1):(length(sX) + length(time_dummies)*iG)] #Takes the other ones
  
  lResults <- fBM_ALg1_inner_loop(dta_s, df_cluster, sX, sY, theta_0, alpha_0, iG, sOLS_formula,
                                  time_dummies, group_dummies, max_iter=100)
  
  return(lResults)
}

fBM_ALg2 <- function(dta_s, df_cluster, sX, sY, iG){
  iN <- length(unique(df_cluster$code))
  n_neigh_jumpsize <- iN/10
  n_iter_LS_max <- 20
  n_neigh_max <- iN
  j_iter_max <- 25
  iT <- length(unique(dta_s$year))
  n_neigh <- 1
  # iN <- length(unique(dta_s$code))
  
  time_dummies <- paste0('time_', sort(unique(dta_s$year)))
  group_dummies <- paste0('G_', 1:iG)
  # Time dummies are equal to 1 if the year is equal to their year
  for(i in c(time_dummies)){
    sYear <- substr(i, 6, 9)
    dta_s[i] <- as.numeric(dta_s$year == as.numeric(sYear))
  }
  for(i in group_dummies){
    iCluster <- as.numeric(substr(i, 3, 3))
    dta_s[i] <- as.numeric(dta_s$code %in% df_cluster$code[df_cluster$cluster == iCluster])
  }
  
  # Inital OLS step
  s_interaction_terms <- ''
  for (i in group_dummies){
    s_interaction_terms <- paste0(paste0(s_interaction_terms, ' + '),
                                  paste0('I(', paste(time_dummies, i, sep='*', collapse=') + I('), ')'))
  }
  sOLS_formula <- paste0(sY, ' ~ ', paste0(sX, collapse=' + 0 + '), ' + ', s_interaction_terms)
  model_1st_step <- lm(sOLS_formula, dta_s)
  
  old_param <- coef(model_1st_step)
  old_param[is.na(old_param)] <- 0
  change_param <- 1
  
  gamma_star <- df_cluster
  objective_func_star <- 1e10
  n_neigh <- 1
  j <- 1
  n_new_optimum <- 0
  
  while(j < j_iter_max){
    # Neighborhood jump
    df_cluster[sample(1:nrow(df_cluster), n_neigh, replace=FALSE), 2] <- sample(1:iG, size=n_neigh, replace=TRUE)
    
    # Starting values:
    theta_0 <- coef(model_1st_step)[sX]
    alpha_0 <- coef(model_1st_step)[(length(sX)+1):(length(sX) + length(time_dummies)*iG)] #Takes the other ones
    
    lResults_Alg1 <- fBM_ALg1_inner_loop(dta_s, df_cluster, sX, sY, theta_0, alpha_0, iG, sOLS_formula,
                                    time_dummies, group_dummies, max_iter=100)
    df_cluster <- lResults_Alg1[[2]]
    model_1st_step <- lResults_Alg1[[1]]
    
    theta_s <- coef(model_1st_step)[sX]
    alpha_s <- coef(model_1st_step)[(length(sX)+1):(length(sX) + length(time_dummies)*iG)] #Takes the other ones
    
    #############################
    # Local Search
    # first part of the error:
    error_1 <- dta_s[sY] - (as.matrix(dta_s[sX]) %*% theta_s)
    error_g <- list()
    for (g in 1:iG){
      alpha_g <- alpha_s[(iT*(g-1) + 1):(iT*g)] # MUST BE IN ORDER
      error_g[[g]] <- as.matrix(dta_s[time_dummies]) %*% alpha_g
    }
    clustering_change <- 1
    n_iter_LS <- 1
    while((clustering_change != 0) & (n_iter_LS < n_iter_LS_max)){
      old_clustering <- df_cluster['cluster']
      clustering_error <- c()
      for (i in 1:nrow(df_cluster)){
        code <- df_cluster$code[i]
        # Take the error:
        error_1_code <- error_1[dta_s$code == code,]
        error_tot_g <- c()
        for (g in 1:iG){
          error_g_code <- error_g[[g]][dta_s$code == code,]
          error_tot_g[g] <- mean((error_1_code - error_g_code)^2)
        }
        df_cluster[i, 'cluster'] <- which.min(error_tot_g) # New assginment
        clustering_error[i] <- min(error_tot_g)
      }
      clustering_change <- sum(old_clustering != df_cluster['cluster'])
    }
    gamma_primeprime <- df_cluster
    objective_func_prime <- mean(clustering_error)
    cat('\r',
        paste0('n_neigh: ', as.character(n_neigh),
               '; j: ', as.character(j),
               '; n_iter_LS: ', as.character(n_iter_LS),
               '; n_new_optimum: ', as.character(n_new_optimum),
               '; obj. diff.: ', as.character(objective_func_prime - objective_func_star)))
    flush.console() 
    # print(as.character(objective_func_prime - objective_func_star))
    
    #############################
    if (objective_func_prime < objective_func_star){
      # print('New optimum')
      n_new_optimum <- n_new_optimum + 1
      n_neigh <- 1
      objective_func_star <- objective_func_prime
      gamma_star <- gamma_primeprime
    } else {
      n_neigh <- floor(n_neigh + n_neigh_jumpsize)
      df_cluster <- gamma_star
    }
    if (n_neigh > n_neigh_max){
      n_neigh <- 1
      j <- j + 1
    }
  }
  
  return(lResults_Alg1)
}

fBM_ALg1_clustered_coefs <- function(dta_s, df_cluster, sX, sY, iG){
  iT <- length(unique(dta_s$year))
  # iN <- length(unique(dta_s$code))
  
  time_dummies <- paste0('time_', sort(unique(dta_s$year)))
  group_dummies <- paste0('G_', 1:iG)
  # Time dummies are equal to 1 if the year is equal to their year
  for(i in c(time_dummies)){
    sYear <- substr(i, 6, 9)
    dta_s[i] <- as.numeric(dta_s$year == as.numeric(sYear))
  }
  for(i in group_dummies){
    iCluster <- as.numeric(substr(i, 3, 3))
    dta_s[i] <- as.numeric(dta_s$code %in% df_cluster$code[df_cluster$cluster == iCluster])
  }
  
  # Inital OLS step
  s_interaction_terms <- ''
  for (i in group_dummies){
    s_interaction_terms <- paste0(paste0(s_interaction_terms, ' + '),
                                  paste0('I(', paste(time_dummies, i, sep=' * ', collapse=') + I('), ')'))
  }
  # Clusterd coefs:
  s_interaction_clustered_coefs <- ''
  for (i in group_dummies){
    s_interaction_clustered_coefs <- paste0(paste0(s_interaction_clustered_coefs, ' + '),
                                  paste0('I(', paste(sX, i, sep=' * ', collapse=') + I('), ')'))
  }
  sOLS_formula <- paste0(sY, ' ~ ', ' 0 ', s_interaction_clustered_coefs, ' + ', s_interaction_terms)
  model_1st_step <- lm(sOLS_formula, dta_s)
  
  old_param <- coef(model_1st_step)
  old_param[is.na(old_param)] <- 0
  change_param <- 1
  
  theta_0 <- coef(model_1st_step)[strsplit(s_interaction_clustered_coefs, ' \\+ ')[[1]][-1]]
  # alpha_0 <- coef(model_1st_step)[(length(sX)+1):(length(sX) + length(time_dummies)*iG)] #Takes the other ones
  alpha_0 <- coef(model_1st_step)[strsplit(s_interaction_terms, ' \\+ ')[[1]][-1]]
  
  lResults <- fBM_ALg1_inner_loop(dta_s, df_cluster, sX, sY, theta_0, alpha_0, iG, sOLS_formula,
                                  time_dummies, group_dummies, max_iter=100,
                                  clustered_coefs=TRUE,
                                  s_interaction_clustered_coefs=s_interaction_clustered_coefs,
                                  s_interaction_terms=s_interaction_terms)
  
  return(lResults)
}