

test = data[1:100]

################################################## get all single copy pifams #################################################
res = list()
swtch = FALSE
nms = c()
nms1 = c()
for(i in 1:length(test)){
    tmp = list()
    for(j in 1:length(test[[i]]$GENOME)){
        x = rle(sort(test[[i]]$GENOME[[j]]@values))
        x = x$values[x$lengths == 1]
        if(length(x) > 0){
            for(n in 1:length(x)){
                isIn = FALSE
                for (v in 1:(length(test[[i]]$GENOME))) {
                    if(v != j){
                        if(x[n] %in% test[[i]]$GENOME[[v]]@values){
                            isIn = TRUE
                        }
                    }
                }
                if(!isIn){
                    tmp2 = list()
                    tmp2[[length(tmp2) +1]] = j
                    tmp2[[length(tmp2) +1]] = test[[i]]$GENOME[[j]]@lengths[which(test[[i]]$GENOME[[j]]@values == x[n])]
                    names(tmp2) = c("GENOMEvector#","baseNr")
                    tmp [[length(tmp) +1]] = tmp2
                    nms1 = c(nms1,as.character(x[n]))
                    swtch = TRUE
                }
            }
        }
    }
    if(swtch){
        names(tmp) = nms1
        nms1 = c()
        res[[length(res) +1]] = tmp
        nms = c(nms,names(test)[i])
        swtch = FALSE
    }
}
names(res) = nms

###########################################################################################################################


chromIDs = c()
j =1
for(i in 1:length(bins)){
    chromIDs[j] = bins[[i]][[1]]$compChromID
    j = j+1
    if(length(bins[[i]]) > 1){
        chromIDs[j] = bins[[i]][[2]]$contChromID
        j = j+1
    }
}

bins = completenessTestData(test,catalogue,15000,25000,1000)

moreInvest = data.table()
for(j in 1:length(bins)){
    w = which(bins[[j]][[1]]$compPifamNames %in%  names(res[as.character(bins[[j]][[1]]$compChromID)][[1]]) )
    
    toInvestigate = data.table(pifams = bins[[j]][[1]]$compPifamNames[w],count = bins[[j]][[1]]$compPifamCount[w],bases = bins[[j]][[1]]$compBaseCount[w])
    
    for(i in 1:length(toInvestigate$pifams)){
        y = which(names(res[as.character(bins[[j]][[1]]$compChromID)][[1]]) == as.character(toInvestigate$pifams[i]))
        if(length(y) > 0){
            if(toInvestigate$count[i] > 1){
                print("Fehler!")
                print(i)
                moreInvest = rbind(moreInvest,toInvestigate[i,])
            }
            if(toInvestigate$bases[i] > res[as.character(bins[[j]][[1]]$compChromID)][[1]][y][[1]]$baseNr){
                print("Fehler!")
                print(i)
                moreInvest = rbind(moreInvest ,toInvestigate[i,])
            }
        }
    }
}

buildUp = c()
x = 0
for(i in 1:length(data$`763146576`$ORF[[6]]@lengths)){
    x = x + data$`763146576`$ORF[[6]]@lengths[i]
    buildUp[i] = x
}

for(i in 1:100){
    conts = getContigsAsIRanges(test,catalogue,15000,25000,1000,i)
    for(j in 1:length(conts)){
        print(conts[j][(end(conts[[j]]) -start(conts[[j]]) < 0)])
    }
}


conts = list()
for(i in 1:length(data)){
    conts[[length(conts)+1]] = vector(mode = "integer")
    conts[[length(conts)+1]] = vector(mode = "integer")
    conts[[length(conts)+1]] = as.integer(names(data)[i])
    conts[[length(conts)+1]] = c(1)
    conts[[length(conts)+1]] = c(sum(data[[i]]$GENOME[[1]]@lengths))
}

acc = c(0:(length(data)-1))

x = countPifams(pifams,Orfs,conts,acc)

unq = list()
for(i in 1:length(data[[1]]$GENOME)){
     unq[[i]] = data[[1]]$GENOME[[i]]@values
}
unq = unique(unlist(unq))
 
sum(unq %in% x[[1]][[1]]$compPifamNames) == sum(x[[1]][[1]]$compPifamNames %in% unq)
bn = c
for(j in 1:length(unq)){
    b = list()
    for(i in 1:length(data[[1]]$GENOME)){
        b[[i]] = data[[1]]$GENOME[[i]]@lengths[data[[1]]$GENOME[[i]]@values == unq[j]]
    }
    bn[j] = sum(unlist(b))
}

for(i in 1:length(x[[1]][[1]]$compBaseCount)){
    n = which(unq == x[[1]][[1]]$compPifamNames[i])
    if(bn[n] != x[[1]][[1]]$compBaseCount[i]){
        print(paste(bn,x[[1]][[1]]$compBaseCount[i]))
    }
}


for(i in 1:length(c)){
    for(j in 1:length(c[[i]])){
        for(n in 1:length(c[[i]][[j]])){
            for(x in 1:length(c[[i]][[j]][[n]])){
                print(i)
            }
        }
    }
}

for(i in 1:100){
    print(i)
    c = getContigsAsIRanges(data,catalogue,10000,15000,44,1)
}



############################################ see weather any completeness or contamination is not in intervall #########################################

compInt = c(0.7,0.9)
contInt = c(0.1,0.4)

x = completenessTestData(data,catalogue,10000,15000,10000,seed = 2,comp = comtInt,cont = contInt)

for(i in 1:length(x)){
    if(x[[i]][[1]]$completness < compInt[1] || x[[i]][[1]]$completness > compInt[2]){
        print("completeness: ")
        print(i)
        print(x[[i]][[1]]$completness)
    }
    if(length(x[[i]]) > 1){
        if(x[[i]][[2]]$contamination > contInt[2] ||x[[i]][[2]]$contamination < contInt[1] ){
            print("contamination:")
            print(i)
            print(x[[i]][[2]]$contamination)
        }
    }
}


