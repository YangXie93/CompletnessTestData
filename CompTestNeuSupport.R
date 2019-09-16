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





#################### format data ###############

data = readRDS("~/work/Data/data.Rds")

res = list()
for(i in 1:length(data)){
    tmp1 = data.table::data.table()
    pfam = data[[i]]$GENOME
    orf = data[[i]]$ORF
    for(j in 1:length(data[[i]]$GENOME)){
        starts = sapply(1:length(pfam[[j]]@lengths),function(x) sum(pfam[[j]]@lengths[1:x]) - pfam[[j]]@lengths[x] +1)
        ends = sapply(1:length(pfam[[j]]@lengths),function(x) sum(pfam[[j]]@lengths[1:x]))
        pfams = pfam[[j]]@values
        orfStarts =  sapply(1:length(orf[[j]]@lengths),function(x) sum(orf[[j]]@lengths[1:x]) - orf[[j]]@lengths[x] +1)
        orfEnds = sapply(1:length(orf[[j]]@lengths),function(x) sum(orf[[j]]@lengths[1:x]))
        ORF = orf[[j]]@values
        pRange = IRanges(start = starts,end = ends,names = pfams)
        oRange = IRanges(start = orfStarts,end = orfEnds,names = ORF)
        
        pRange = subset(pRange,as.integer(names) > 0)
        oRange = subset(oRange,as.integer(names) > 0)
        
        tr = findOverlaps(oRange,pRange)
        ORF = names(oRange)[queryHits(tr)]
        tmp = data.table::data.table(start = start(pRange),end = end(pRange),pfam = names(pRange),ORF = ORF)
        if(length(tmp1) == 0){
            tmp1 = tmp
        }
        else{
            tmp1 = rbind(tmp1,tmp)
        }
    }
    res[[i]] = tmp1
    print(paste(i,"von:",length(data)))
}
names(res) = names(data)

saveRDS(res,"~/work/Data/newDatafull.Rds")

cata = readRDS("~/work/Data/cata.Rds")

###############################################################################

for(i in 1:length(x)){
    tmp = 0
    tmp1 = 0
    tmp2 = 0
    pCount = 0
    pCountOrg = 0
    wrongPfams = list()
    for(j in 1:length(x[[i]]$completeness$chromID)){
        tmp = tmp + sum(data[[ x[[i]]$completeness$chromID[j] ]]$end - data[[ x[[i]]$completeness$chromID[j] ]]$start +1)
        tmp1 = tmp1 + max(data[[ x[[i]]$completeness$chromID[j] ]]$end)
        tmp2 = tmp2 + sum(width(conts[[i]][[1]][[j]]))
        pCountOrg = pCountOrg + sum(rle(sort(as.integer(data[[ x[[i]]$completeness$chromID[j] ]]$pfam)))$lengths)
        pCount = pCount + sum(x[[i]]$completeness$pfamCount)
        wrongPfams[[j]] = data[[ x[[i]]$completeness$chromID[j] ]]$pfam
    }
    wrongPfams = unique(unlist(wrongPfams))
    if(sum(x[[i]]$completeness$baseCount)/tmp -x[[i]]$completeness$completeness > 0.1 || sum(x[[i]]$completeness$baseCount)/tmp -x[[i]]$completeness$completeness < -0.1){
        print(paste("bCount",i))
        print(length(x[[i]]$completeness$chromID))
        print(paste(sum(x[[i]]$completeness$baseCount)/tmp,tmp2/tmp1,x[[i]]$completeness$completeness))
    }
    if(pCount/pCountOrg - x[[i]]$completeness$completeness > 0.1 || pCount/pCountOrg - x[[i]]$completeness$completeness < -0.1 ){
        print(paste("pCount",i))
        print(length(x[[i]]$completeness$chromID))
        print(paste(pCount/pCountOrg,x[[i]]$completeness$completeness))
    }
    if(length(x[[i]]$completeness$pfamName)/length(wrongPfams) -x[[i]]$completeness$completeness > 0.1 ||length(x[[i]]$completeness$pfamName)/length(wrongPfams) -x[[i]]$completeness$completeness < -0.1 ){
        print(paste("Nr Unique ps",i))
        print(length(x[[i]]$completeness$chromID))
        print(paste(length(x[[i]]$completeness$pfamName)/length(wrongPfams),x[[i]]$completeness$completeness))
    }
}

# 391
# 972

