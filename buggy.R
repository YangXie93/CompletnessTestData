rm(list = ls())
Rcpp::sourceCpp('countPifamsV2.2.cpp')
source("Com.R")
data = readRDS("~/Daten/data.Rds")
catalogue = readRDS("~/Daten/bac_dt.Rds")$DD

x = completenessTestData(data,catalogue,10000,20000,1000)



