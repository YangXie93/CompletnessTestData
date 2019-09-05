rm(list = ls())
library(CompletenessTestData)
data = readRDS("~/Daten/data.Rds")
catalogue = readRDS("~/Daten/bac_dt.Rds")$DD
for(i in 1:100){
    x = completenessTestData(data,catalogue,10000,20000,1000)
}
