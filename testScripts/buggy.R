rm(list = ls())
library(CompletenessTestData)
library(IRanges)
data = readRDS("~/Daten/data_red.Rds")
cat = readRDS("~/Daten/bac_dt.Rds")$DD


x = completenessTestData(data,cat,10000,15000,10000,seed =1)
print(x[[44]])
