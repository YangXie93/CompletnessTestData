rm(list = ls())
Rcpp::sourceCpp('src/countPifamsV2.1.cpp')
source('~/CompletnessTestData/R/completenessTestData.R')
data = readRDS("~/Daten/data.Rds")
catalogue = readRDS("~/Daten/bac_dt.Rds")$DD

x = completenessTestData(data,catalogue,10000,20000,1000,seed = 1)
