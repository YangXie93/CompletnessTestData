completenessTestData <- function(data,catalogue,minContigLength,meanContigLength,number,seed = 0,distr = "normal",comp = c(0.6,1.0),cont = c(0.0,0.4)){
    x = Sys.time()
    cat = subset(catalogue,GI.Vec %in% names(data))
    
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
    for(i in 1:length(tmp1)){
        tmp2 = c()
        lengthSums[i] = 0
        for(j in 1:length(tmp1[[i]])){
            tmp2[length(tmp2)+1] = as.integer(tmp1[[i]][j])
            tmp2[length(tmp2)+1] = sum(data[[tmp1[[i]][j]]]$GENOME[[1]]@lengths)
            lengthSums[i] = lengthSums[i] + tmp2[length(tmp2)]
        }
        lengths[[i]] = tmp2
    } 
    pifams = list()
    Orfs = list()
    for(i in 1:length(data)){
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
    res = compTestData(pifams,Orfs,lengths,lengthSums,minContigLength,meanContigLength,number,comp,cont,seed,distr)
    print(Sys.time() -x)
    return(res)
}




timeCompletenessTestData <- function(times,data,catalogue,minContigLength,meanContigLength,seed = 0){
    time = Sys.time()

    x = completenessTestData(data,catalogue,minContigLength,meanContigLength,times,seed)

    print(Sys.time() - time)
}



meanTimeCompletenessTestData <- function(meanTime,times,data,catalogue,minContigLength,meanContigLength,seed = 0){
    mn = c()
    for(i in 1:meanTime){    
        time = Sys.time()
        
        x = completenessTestData(data,catalogue,minContigLength,meanContigLength,times,seed)
        
        mn[i] = Sys.time() - time
    }
    res = list()
    res[[1]] = mean(mn)
    res[[2]] = sd(mn)
    return(res)
}