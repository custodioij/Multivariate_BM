library(dplyr)
library(tidyr)
d_convergence_threshold =  1e-10

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

fBM_ALg1_inner_loop_timevar <- function(dta_s, sX, sY, theta_0, alpha_0,
                                       iG, nPartitions, sOLS_formula,
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
    # if (clustered_coefs){
    #   error_1 <- dta_s[sY]-0
    #   error_g <- list()
    #   for (g in 1:iG){
    #     # First the time dummies
    #     alpha_g <- alpha_s[grep(paste0('G_', g), names(alpha_s))]
    #     # Make sure they have the correct order
    #     names(alpha_g) <- gsub('^I\\((time_[0-9]{4}).*$', '\\1' , names(alpha_g))
    #     error_g[[g]] <- as.matrix(dta_s[time_dummies]) %*% alpha_g[time_dummies]
    #     # Second the clustered coefs
    #     theta_g <- theta_s[grep(paste0('G_', g), names(theta_s))]
    #     names(theta_g) <- gsub('^I\\((.*) \\* G\\_[0-9]\\)$', '\\1' , names(theta_g))
    #     error_g[[g]] <- error_g[[g]] + as.matrix(dta_s[sX]) %*% theta_g[sX]
    #   }
    # }else{
    #}
    error_1 <- dta_s[sY] - (as.matrix(dta_s[sX]) %*% theta_s)
    error_g <- list()
    for (g in 1:iG){
      # First the time dummies
      alpha_g <- alpha_s[grep(paste0('g_', g), names(alpha_s))]
      # Make sure they have the correct order
      names(alpha_g) <- gsub('^I\\((t_[0-9]{4}).*$', '\\1' , names(alpha_g))
      error_g[[g]] <- as.matrix(dta_s[time_dummies]) %*% alpha_g[time_dummies]
    }
    # To make the partition loop easier: (there should be a betetr way)
    error_1_all <- cbind.data.frame(error_1, dta_s[, 'code'])
    names(error_1_all) <- c('error', 'code')
    
    #########
    # Re shuffle cluster independently in each partition
    #########
      
    old_clustering <- dta_s['cluster']
    clustering_error <- c()
    for (p in 1:nPartitions){
      p_subset <- dta_s$partition == p
      local_cluster_p <- distinct(dta_s[p_subset, c('code', 'cluster')])
      for (i in 1:nrow(local_cluster_p)){
        code <- local_cluster_p$code[i]
        p_i_subset <- (dta_s$code == code) & p_subset
        # Take the error:
        error_1_code <- error_1[p_i_subset,]
        error_tot_g <- c()
        for (g in 1:iG){
          error_g_code <- error_g[[g]][p_i_subset,]
          error_tot_g[g] <- mean((error_1_code - error_g_code)^2)
        }
        dta_s[p_i_subset, 'cluster'] <- which.min(error_tot_g)  # New assginment
        # clustering_error[i] <- min(error_tot_g)
      }
    }
    # clustering_change <- sum(old_clustering != dta_s['cluster'])
    
    #########
    # Probably should redo OLS here!
    #########
    
    # Redo the cluster variable in the main dataset
    for(i in group_dummies){
      iCluster <- as.numeric(substr(i, 3, 3))
      dta_s[i] <- as.numeric(dta_s$cluster == iCluster)
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
    
    #########
    # Probably should redo OLS here!
    #########
    
    old_partitioning <- dta_s['partition']
    # clustering_error <- c()
    all_partitions <- distinct(dta_s[c('partition','year')])
    cluster_by_partition <- distinct(dta_s[c('code','partition','cluster')])
    # Fill the holes in cluster_by_partition so that we don't drop countries that
    # only appear in one partition.
    cluster_by_partition <- merge(expand.grid(code=unique(dta_s$code),
                                              partition=1:nPartitions),
                                  distinct(dta_s[c('code','partition','cluster')]),
                                  all.x=TRUE, sort=FALSE)
    # The NAs show where a country has no observations in that partition and so it does
    # not have an assigned cluster. 
    # Fill it by copying a cluster from another partition
    # TODO: find a more efficient way of doing this
    cluster_by_partition <- cluster_by_partition %>% group_by(code) %>%
      summarise(cl1 = first(cluster)) %>% merge(cluster_by_partition)
    cluster_by_partition[is.na(cluster_by_partition$cluster), 'cluster'] <- 
      cluster_by_partition[is.na(cluster_by_partition$cluster), 'cl1']
    cluster_by_partition <- select(cluster_by_partition, -cl1)
    # Check what would be the error if moving each t to a different partition.
    # For that we need to calculate the error based on the cluster of each country
    # At the other partitions.
    for (t in all_partitions$year){
      t_subset <- dta_s$year == t
      error_1_t <- error_1_all[t_subset,]
      # Take the error:
      error_tot_p <- c()
      # Use the fact that error_g[[g]] is constant over t
      for (p in 1:nPartitions){
        clusters_p = cluster_by_partition[cluster_by_partition$partition == p,
                                          c('cluster', 'code')]
        clusters_p$error <- NA
        for (g in 1:iG){
          alpha_tg <- alpha_s[paste0('I(t_', t, ' * g_', g, ')')]
          clusters_p[clusters_p$cluster == g, 'error'] <- alpha_tg
        }
        error_p <- merge(error_1_t, clusters_p,
                         by='code', sort=FALSE)[,c('error.x', 'error.y')]
        error_tot_p[p] <- mean((error_p$error.x - error_p$error.y)^2)
      }
      
      dta_s[t_subset, 'partition'] <- which.min(error_tot_p)  # New assginment
      # Reajust clustering
      dta_s$cluster <- NULL
      dta_s <- merge(dta_s, cluster_by_partition, sort=FALSE)
    }
    # clustering_change <- sum(old_partitioning != dta_s['partition'])
    
    #########
    # Back to OLS
    #########
    
    # Redo the cluster variable in the main dataset
    for(i in group_dummies){
      iCluster <- as.numeric(substr(i, 3, 3))
      dta_s[i] <- as.numeric(dta_s$cluster == iCluster)
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
  clustering_error <- NULL # to be removed probably
  return(list(model_1st_step, dta_s, clustering_error, theta_s, alpha_s, n_iter))
}

fBM_ALg1_timevar <- function(dta, sX, sY, iG=2, nPartitions=2, n_starts=100){
  iT <- length(unique(dta$year))
  # iN <- length(unique(dta_s$code))
  
  llResults <- list()
  total_error <- c()
  for (i_start in 1:n_starts){
    dta_s <- dta
    # Generate random time partitions and clustering
    
    df_partitions <- unique(dta_s[c('year')])
    df_partitions$partition <- sample(x=rep_len(1:nPartitions, iT), size=iT)
    df_cluster <- c()
    for (i in 1:nPartitions){
      df_cluster_p <- unique(dta_s[c('code')])
      df_cluster_p$partition <- i
      df_cluster_p$cluster <- sample(x=rep_len(1:iG, nrow(df_cluster_p)),
                                     size=nrow(df_cluster_p))
      df_cluster <- rbind(df_cluster, df_cluster_p)
    }
    df_cluster <- merge(df_cluster, df_partitions)
    dta_s <- merge(dta_s, df_cluster)
    
    
    time_dummies <- paste0('t_', sort(unique(dta_s$year)))
    group_dummies <- paste0('g_', 1:iG)
    # Time dummies are equal to 1 if the year is equal to their year
    for(i in c(time_dummies)){
      sYear <- substr(i, 3, 6)
      dta_s[i] <- as.numeric(dta_s$year == as.numeric(sYear))
    }
    for(i in group_dummies){
      iCluster <- as.numeric(substr(i, 3, 3))
      dta_s[i] <- as.numeric(dta_s$cluster == iCluster)
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
    
    llResults[[i_start]] <- list()
    # tryCatch({
    llResults[[i_start]] <- fBM_ALg1_inner_loop_timevar(dta_s, sX, sY,
                                           theta_0, alpha_0, iG, nPartitions,
                                           sOLS_formula,
                                           time_dummies, group_dummies, max_iter=100)
    total_error[i_start] <- mean(llResults[[i_start]][[1]]$residuals^2)
    # },
    # error={function(ex){
    #   print('error:')
    #   print(ex)
    #   llResults[[i_start]] <- as.list(rep(NA, 6))
    #   }}
    # )
  }
  # Keep the one with best regression fit
  lResults <- llResults[[which.min(total_error)]]
  lResults[[length(lResults)+1]] <- total_error
  
  return(lResults)
}

