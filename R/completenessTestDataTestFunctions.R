#' @export
testBaseCountcomp<- function(x,data){
    res = c()
    for(i in 1:length(x)){
        y = sum(x[[i]][[1]]$compBaseCount)
        z = sum(c(data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[2]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@values >= 0]))
        res[i] = y/z-x[[i]][[1]]$completness
    }
    return(res)
}
#' @export
testBaseCountcont<- function(x,data){
    res = c()
    for(i in 1:length(x)){
        if(length(x[[i]]) > 1){
            y = sum(x[[i]][[2]]$contBaseCount)
            z = sum(c(data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@values >= 0]))
            res[i] = y/z-x[[i]][[2]]$contamination
        }
    }
    return(res)
}

#' @export
testPifamCountComp <- function(x,data){
    res = c()
    for(i in 1:length(x)){
        y = sum(x[[i]][[1]]$compPifamCount)
        z = sum(rle(sort(c(data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@values >= 0])))$lengths)
        res[i] = y/z-x[[i]][[1]]$completness
    }
    return(res)
}
#' @export
testPifamCountcont <- function(x,data){
    res = c()
    for(i in 1:length(x)){
        if(length(x[[i]]) > 1){
            y = sum(x[[i]][[2]]$contPifamCount)
            z = sum(rle(sort(c(data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@values >= 0])))$lengths)
            res[i] = y/z-x[[i]][[2]]$contamination
        }
    }
    return(res)
}


# die funktion guckt ob die completeness und contamination werte inerhalb ihres intervalls liegen
#' @export
contAndCompInIntervall <- function(x,compInt,contInt){
    for(i in 1:length(x)){
        if(x[[i]][[1]]$completness < compInt[1] || x[[i]][[1]]$completness > compInt[2]){
            print(paste("completness out of bounds at: ",i,"value: ",x[[i]][[1]]$completness))
        }
        if(x[[i]][[2]]$contamination < contInt[1] || x[[i]][[2]]$contamination > contInt[2]){
            print(paste("contamination out of bounds at: ",i))
        }
    }
}

#' @export
meanDiffPiCountAndBaseCount <- function(data,cat,minContigLength,meanContigLength,nrExmlps,startingSeed = 0,times){
    tmpBase = c()
    tmpPifam = c()
    tmpBaseDiv = c()
    tmpPifamDiv = c()
    for(i in 1:times){
        x = completenessTestData(data,cat,minContigLength,meanContigLength,nrExmlps,seed = startingSeed+i)
        tmpBase[i] = mean(c(testBaseCountcomp(x,data),testBaseCountcont(x,data)))
        tmpBaseDiv = append(tmpBaseDiv,c(testBaseCountcomp(x,data),testBaseCountcont(x,data)))
        tmpPifam[i] = mean(c(testPifamCountComp(x,data),testPifamCountcont(x,data)))
        tmpPifamDiv = append(tmpPifamDiv,c(testPifamCountComp(x,data),testPifamCountcont(x,data)))
    }
    print(paste("mean baseCount diff: ",mean(tmpBase)))
    print(paste("sd baseCount diff: ",sd(tmpBaseDiv)))
    print(paste("max base diff: " ,max(tmpBaseDiv),"min base diff: ",min(tmpBaseDiv)))
    print(paste("meen pifamCount diff: ",mean(tmpPifam)))
    print(paste("sd pifamCount diff: ",sd(tmpPifamDiv)))
    print(paste("max pifam diff: " ,max(tmpPifamDiv),"min pifam diff: ",min(tmpPifamDiv)))

}

#' @export
testForSpeedDependenceOnData <- function(datadir,cat,minContigLength,meanContigLength,nrExmpls){
    data = dir(datadir)
    data = data[substr(data,1,3) == "IDz"]
    t = c()
    l = c()
    for(i in 1:length(data)){
        dt = readRDS(paste(datadir,"/",data[i],sep = ""))
        l[i] = length(dt)
        print(paste(data[i]," length: ",l[i]))
        tmp = Sys.time()
        x = completenessTestData(dt,cat,minContigLength,meanContigLength,nrExmpls)
        t[i] = Sys.time() - tmp
    }
    return(t)
}


#' @export
plotIRanges <-function(ra){
    dat = as.data.table(ra)
    dat = cbind(dat,bin = bins)
    ggplot(dat) + 
        geom_rect(aes(xmin = start, xmax = end,
                      ymin = 1, ymax = 1.9)) +
        theme_bw()
}

#' @export
getContigsAsIRanges <- function(data,catalogue,minContigLength,meanContigLength,number,seed){
    
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
    
    
    n = vector(mode = "integer")
    c = CompletenessTestData::mkContigs(lengths,lengthSums,minContigLength,meanContigLength,number,c(0.6,1),c(0,0.4),accession,n,seed = seed)
    
    cc = list()
    i = 1
    while(i <= length(c)){
        if(length(c[[i]]) > 0){
            if(length(c[[i +1]]) > 0){
                tc = list()
                while(i <= length(c) && length(c[[i]]) > 0){
                    tc[[length(cc)+1]] = IRanges(start = c[[i+1]],end = c[[i+2]],names = rep(c[[i]],length(c[[i+1]])))
                    i = i +3
                }
                cc[[length(cc) +1]] = tc
            }
        }
        i = i +1
    }
    return(cc)
    
}
#' @export
chooseGenomes <- function(data,catalogue,minContigLength,times,seed){
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
    
    return(chooseGenomes(lengthSums,minContigLength,seed = seed,times = times))
}
