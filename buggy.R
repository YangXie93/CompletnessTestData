rm(list = ls())
library(IRanges)
library(data.table)

Rcpp::sourceCpp('CompTestNeu.cpp')
source('compTestNeu.R')
data = readRDS("~/Daten/newData1000.Rds")
cata = readRDS("~/Daten/cata.Rds")

x = getContigsAsIRanges(data,cata,10000,20000,number = 100,seed = 1,debugInfo = TRUE)
