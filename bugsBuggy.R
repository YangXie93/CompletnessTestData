rm(list = ls())
library(CompletenessTestData)
zero = readRDS("~/Daten/data.Rds")
cat = readRDS("~/Daten/bac_dt.Rds")$DD
for(i in 6:100){
    x = completenessTestData(zero,cat,10000,20000,1000,seed = i)
}
