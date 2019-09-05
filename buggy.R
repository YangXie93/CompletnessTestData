
library(CompletenessTestData)

zero = readRDS("~/Daten/genome_R/IDz_8_9.Rds")
cat = readRDS("~/Daten/bac_dt.Rds")$DD


t = testForSpeedDependenceOnData("~/Daten/genome_R",cat,10000,20000,1000)

