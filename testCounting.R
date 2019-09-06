Rcpp::sourceCpp("testcountPifamsV2.1.cpp")

data = readRDS("~/Daten/IDz_Zero_Zero.Rds")
catalogue = readRDS("~/Daten/bac_dt.Rds")$DD

minContigLength = 10000
meanContigLength = 20000
number = 1
comp = c(0.2)
cont = c(0.1)
distr = ""
seed = 0
data = data[unique(names(data))]
cat = subset(catalogue,GI.Vec %in% names(data))
nms = names(data) %in% cat$GI.Vec
IDs = names(data)[nms]
tmp1 = list()
isUsed = c()
for(i in 1:length(cat$comb)){
    if(!(cat$comb[i] %in% isUsed)){
        tmp1[[length(tmp1)+1]] = cat$GI.Vec[which(cat$comb == cat$comb[i])]
        isUsed[length(isUsed)+1] = cat$comb[i]
        
    }
}
lengths = list()
lengthSums = c()
accession = list()
for(i in 1:length(tmp1)){
    tmp2 = c()
    tmp5 = c()
    lengthSums[i] = 0
    for(j in 1:length(tmp1[[i]])){
        tmp5[j] = which(IDs == tmp1[[i]][j]) -1
        tmp2[length(tmp2)+1] = as.numeric(tmp1[[i]][j])
        tmp2[length(tmp2)+1] = sum(data[[tmp1[[i]][j]]]$GENOME[[1]]@lengths)
        lengthSums[i] = lengthSums[i] + tmp2[length(tmp2)]
    }
    lengths[[i]] = tmp2
    accession[[i]] = tmp5
} 
pifams = list()
Orfs = list()
for(i in 1:length(data)){
    if(nms[i]){
        tmp3 = list()
        tmp4 = list()
        for(j in 1:length(data[[i]]$GENOME)){
            if((!is.null(data[[i]]$GENOME[[j]])) && (!is.null(data[[i]]$ORF[[j]]))){
                tmp3[[length(tmp3)+1]] = data[[i]]$GENOME[[j]]@lengths
                tmp3[[length(tmp3)+1]] = data[[i]]$GENOME[[j]]@values
                tmp4[[length(tmp4)+1]] = data[[i]]$ORF[[j]]@lengths
                tmp4[[length(tmp4)+1]] = data[[i]]$ORF[[j]]@values
            }
        }
        pifams[[length(pifams)+1]] = tmp3
        Orfs[[length(Orfs)+1]] = tmp4
    }
}
print(lengths)
if(length(lengths) > 1){
    res = compTestData(pifams,Orfs,lengths,lengthSums,minContigLength,meanContigLength,number,comp,cont,accession,seed,distr)
}
print(res)