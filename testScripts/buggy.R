rm(list = ls())
Rcpp::sourceCpp('src/countPifamsV2.1.cpp')
source('~/CompletnessTestData/R/completenessTestData.R')
data = readRDS("~/Daten/data.Rds")
catalogue = readRDS("~/Daten/bac_dt.Rds")$DD

for(i in 1:100){
    x = completenessTestData(data,catalogue,10000,25000,10000,seed = i)
}