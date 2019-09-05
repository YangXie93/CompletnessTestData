rm(list = ls())
library(CompletenessTestData)
zero = readRDS("~/Daten/data.Rds")
cat = readRDS("~/Daten/bac_dt.Rds")$DD

x = completenessTestData(zero,cat,10,15,10,cont = c(0,0.01),comp = c(0.1,0.2))

