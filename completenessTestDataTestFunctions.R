testBaseCountcomp<- function(x,data){
    res = c()
    for(i in 1:length(x)){
        y = sum(x[[i]][[1]]$compBaseCount)
        z = sum(c(data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[2]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@lengths[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@values >= 0]))
        res[i] = y/z-x[[i]][[1]]$completness
    }
    return(res)
}

testBaseCountcont<- function(x,data){
    res = c()
    for(i in 1:length(x)){
        y = sum(x[[i]][[2]]$contBaseCount)
        z = sum(c(data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@lengths[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@values >= 0]))
        res[i] = y/z-x[[i]][[2]]$contamination
    }
    return(res)
}


testPifamCountComp <- function(x,data){
    res = c()
    for(i in 1:length(x)){
        y = sum(x[[i]][[1]]$compPifamCount)
        z = sum(rle(sort(c(data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[1]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@values[data[[which(names(data) == x[[i]][[1]]$compChromID)]]$GENOME[[6]]@values >= 0])))$lengths)
        res[i] = y/z-x[[i]][[1]]$completness
    }
    return(res)
}

testPifamCountcont <- function(x,data){
    res = c()
    for(i in 1:length(x)){
        y = sum(x[[i]][[2]]$contPifamCount)
        z = sum(rle(sort(c(data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[2]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[3]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[4]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[5]]@values >= 0],data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@values[data[[which(names(data) == x[[i]][[2]]$contChromID)]]$GENOME[[6]]@values >= 0])))$lengths)
        res[i] = y/z-x[[i]][[2]]$contamination
    }
    return(res)
}


# die funktion guckt ob die completeness und contamination werte inerhalb ihres intervalls liegen

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


testForSpeedDependenceOnData <- function(datadir,cat,minContigLength,meanContigLength,nrExmpls,seed){
    data = dir(datadir)
    data = data[substr(data,1,3) == "IDz"]
    t = c()
    l = c()
    for(i in 1:length(data)){
        dt = readRDS(paste(datadir,"/",data[i],sep = ""))
        l[i] = length(dt)
        print(paste(data[i]," length: ",l[i]))
        tmp = Sys.time()
        x = completenessTestData(dt,cat,minContigLength,meanContigLength,nrExmpls,seed = seed)
        t[i] = Sys.time() - tmp
    }
    return(t)
}







