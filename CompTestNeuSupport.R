x = IRanges(start = c(25,10),end = c(36,20))
y =IRanges(start = c(17),end = c(30))


overlapsRanges(y,x)




data = readRDS("~/Daten/data_red.Rds")
cata = readRDS("~/Daten/cata_red.Rds")

n = c()
for(i in 1:length(data[[1]]$GENOME)){
    n = c(n,sapply(1:length(data[[1]]$GENOME[[i]]@lengths), function(x) sum(data[[1]]$GENOME[[i]]@lengths[1:x])))
}

m = c()
for(i in 1:length(data[[1]]$GENOME)){
    m = c(m,sapply(1:length(data[[1]]$GENOME[[i]]@lengths), function(x) sum(data[[1]]$GENOME[[i]]@lengths[1:x]) - data[[1]]$GENOME[[i]]@lengths[x]))
}

pfms = c()
for(i in 1:length(data[[1]]$GENOME)){
    pfms = c(pfms,data[[1]]$GENOME[[i]]@values)
}

t = IRanges(start = m,end = n,names = pfms)
t = subset(t,names(t) > 0)

n = c()
for(i in 1:length(data[[1]]$ORF)){
    n = c(n,sapply(1:length(data[[1]]$ORF[[i]]@lengths), function(x) sum(data[[1]]$ORF[[i]]@lengths[1:x])))
}

m = c()
for(i in 1:length(data[[1]]$ORF)){
    m = c(m,sapply(1:length(data[[1]]$ORF[[i]]@lengths), function(x) sum(data[[1]]$ORF[[i]]@lengths[1:x]) - data[[1]]$ORF[[i]]@lengths[x]))
}

orfs = c()
for(i in 1:length(data[[1]]$ORF)){
    orfs = c(orfs,data[[1]]$ORF[[i]]@values)
}

o = IRanges(start = m,end = n,names = orfs)
o = subset(o,names(o) > 0)

subsetByOverlaps(o,t[1])


ORFs = sapply(1:length(t), function(x) names(subsetByOverlaps(o,t[x])))

res = list()
res[[1]] = data.table(start = start(t),end = end(t),pfam = names(t),ORF = ORFs)






test = IRanges(start = c(100),end = c(100000))
pfamName = names(subsetByOverlaps(t,test))

which(res[[1]]$pfam == pfamName[1])

res[[1]]$ORF

unique(res[[1]]$orf[which(res[[1]]$pfam == pfamName[1])])

length(unique(res[[1]]$ORF[res[[1]]$pfam == pfamName[1]]))


testPCount = sapply(pfamName[1:length(pfamName)], function(x) length(unique(res[[1]]$ORF[res[[1]]$pfam == x])))


 
test1 = IRanges(start = c(10,15),end = c(20,25))
test2 = IRanges(start = c(15),end = c(20))

overlapsRanges(test1,test2)






overlapsRanges(t,contigs[[1]][[1]][[1]])


pp = sapply(1:length(inv$GENOME), function(x) unique(inv$GENOME[[x]]@values))
pp = unlist(pp)
pp =  unique(pp)

length(pruf[[1]]$completeness$pfamName)/length(pp[pp > 0])

length(res[[1]]$completeness$pfamName)/length(pp[pp > 0])


pc =  sapply(1:length(inv$GENOME), function(x) inv$GENOME[[x]]@values)
pc = sort(unlist(pc))
pc = rle(pc)

sum(pruf[[1]]$completeness$pfamCount)/sum(pc$lengths[pc$values >0])
sum(res[[1]]$completeness$pfamCount)/sum(pc$lengths[pc$values >0])

bc = sapply(1:length(inv$GENOME), function(x) inv$GENOME[[x]]@lengths[inv$GENOME[[x]]@values > 0 ])
bc = sum(unlist(bc))

sum(pruf[[1]]$completeness$baseCount)/sum(bc)
sum(res[[1]]$completeness$baseCount)/sum(bc)





contigs = list(list(list(IRanges(start = 1,end = max(testDat[[1]]$end),names = "478476202"),1)))


baseCount = sapply(1:length(pfamName), function(x) sum(width(overlaps[names(pfamRange) == pfamName[x] ])))



pfamCount = sapply(1:length(pfamName), function(x) length(unique(this$ORF[this$pfam == pfamName[x]])))

