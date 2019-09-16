library(IRanges)
library(data.table)

Rcpp::sourceCpp('CompTestNeu.cpp')
source('~/work/CompletnessTestData/compTestNeu.R')
data = readRDS("~/work/Data/newDatafull.Rds")
cata = readRDS("~/work/Data/cata.Rds")

x = pfamCounter(data,cata,10000,20000,number = 1000,seed = 1)

