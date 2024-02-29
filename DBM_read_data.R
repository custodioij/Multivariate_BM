rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)
source('B_M_functions_timevar.R')


sOutDir <- "Figures/"
year_max <- 2000
year_min <- 1970
iB <- 100
lPart <- 1:4
lG <- 2:4

dta_bm <- read_excel('Income-and-Democracy-Data-AER-adjustment.xls', sheet='5 Year Panel')

dta_s <- dta_bm %>% group_by(code) %>%
  dplyr::mutate(lag_fhpolrigaug = lag(fhpolrigaug, order_by=year),
                lag_lrgdpch = lag(lrgdpch, order_by=year)) %>% 
  subset(samplebalancefe == 1)
# write.csv(dta_s, "BM_data.csv")

dta_s <- dta_s[dta_s$year <= year_max & dta_s$year >= year_min,]
dta_s <- dta_s[c('code', 'year', 'fhpolrigaug', 'lag_fhpolrigaug', 'lag_lrgdpch')]

colSums(is.na(dta_s))

sY <- 'fhpolrigaug'
sX <- c('lag_fhpolrigaug', 'lag_lrgdpch')

iT <- length(unique(dta_s$year))
iN <- length(unique(dta_s$code))

for (iPart in lPart){
for (iG in lG){

lResults <- fBM_ALg1_timevar(dta_s, sX, sY, iG, nPartitions=iPart, n_starts=iB)

model_1st_step <- lResults[[1]]
df_cluster <- lResults[[2]]

res_to_save <- lResults
save(res_to_save, file=paste0(sOutDir, 'results_', iG, '_cl_',iPart, '_part.Rdata'))

}
}