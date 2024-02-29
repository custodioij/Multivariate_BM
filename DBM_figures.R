rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(countrycode)
source('B_M_functions_timevar.R')


sOutDir <- "Figures/"
cluster_palette <- brewer.pal(name="YlGnBu",n=5)[2:5]

lPart <- 1:4
lG <- 2:4
dta_bm <- read_excel('Income-and-Democracy-Data-AER-adjustment.xls', sheet='5 Year Panel')

####################
# To make the maps:
world_coordinates_base <- map_data("world")
world_coordinates_base$region <- countryname(world_coordinates_base$region)
dta_names <- dta_bm
dta_names <- distinct(dta_names[c('country', 'code')])
dta_names$country <- countryname(dta_names$country)
# save(dta_names, file='dta_names.Rdata')
# load('dta_names.Rdata')
####################

sY <- 'fhpolrigaug'
sX <- c('lag_fhpolrigaug', 'lag_lrgdpch')

iPart <- 4
iG <- 4
for (iPart in lPart){
  for (iG in lG){
    # Load results:
    load(file=paste0(sOutDir, 'results_', iG, '_cl_',iPart, '_part.Rdata'))
    model_1st_step <- res_to_save[[1]]
    df_cluster <- res_to_save[[2]]
    
    # Align clusters across partitions
    # 1. calculate cluster-partition means
    align_clusters <- df_cluster %>% group_by(cluster, partition) %>%
      summarise(Y_mean = mean(fhpolrigaug)) %>% 
      group_by(partition) %>% mutate(Y_rank = as.integer(rank(Y_mean)))
    align_clusters$new_cluster <- align_clusters$Y_rank
    align_clusters <- align_clusters[c('cluster', 'new_cluster', 'partition')]
    
    df_cluster <- merge(df_cluster, align_clusters) %>%
      select(-cluster) %>% rename(cluster = new_cluster)
    
    # 2. Align partitions by start year
    align_partition <- df_cluster[c('year', 'partition')] %>% distinct() %>% group_by(partition) %>%
      summarise(min_year=min(year)) %>% 
      mutate(new_partition = as.integer(rank(min_year))) %>% select(-min_year)
    
    df_cluster <- merge(df_cluster, align_partition) %>%
      select(-partition) %>% rename(partition = new_partition)
    
    # 3. Name partitions
    years_in_partition <- df_cluster[c('year', 'partition')] %>% distinct()
    partition_names <- data.frame(partition=1:iPart, part_name = 1:iPart)
    for (part in 1:iPart){
      partition_names[partition_names$partition == part, 'part_name'] <- 
        paste0(sort(years_in_partition[years_in_partition$partition == part, 'year']), collapse = ', ')
    }
    
    # Map it again and replace because these are better.
    map_nrows <- iPart #1 + (as.numeric(iPart) > 2)
    # map_ncols <- 2
    # Trick to have the NAs as grey in all partitions:
    world_coordinates_part <- list()
    for (part in sort(unique(df_cluster$partition))){
      world_coordinates_part[[part]] <- world_coordinates_base
      world_coordinates_part[[part]]$partition <- part
    }
    world_coordinates <- bind_rows(world_coordinates_part)
    
    world_coordinates <- merge(world_coordinates,
                               merge(dplyr::select(df_cluster,
                                                   c('cluster', 'partition', 'code')),
                                     # c('cluster', 'code')),
                                     dta_names),
                               by.x=c('region', 'partition'), by.y=c('country','partition'),
                               # by.x=c('region'), by.y=c('country'),
                               all.x=TRUE, sort=FALSE)
    world_coordinates <- merge(world_coordinates, partition_names, sort=FALSE)
    world_coordinates <- arrange(world_coordinates, order)
    world_coordinates$col <- 'white'
    ggplot() +
      geom_map(
        data = world_coordinates, map = world_coordinates,
        aes(long, lat, map_id = region, fill=factor(cluster))) +
      facet_wrap(part_name~., ncol=1) +
      ylim(-50, 70) +
      theme_void() + theme(legend.position = "bottom") + 
      scale_fill_manual(values=cluster_palette, name = "Cluster", na.value="grey")
    ggsave(paste0(sOutDir, 'map', iG, '_cl_',iPart, '_part.png'),
           width = 8, height = 4*map_nrows)
    
  }
}