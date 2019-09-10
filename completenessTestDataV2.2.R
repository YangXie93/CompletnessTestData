#' @useDynLib CompletenessTestData
#' @importFrom Rcpp sourceCpp
#' @export
completenessTestData <- function(data,catalogue,minContigLength,meanContigLength,number,seed = 0,distr = "normal",comp = c(0.6,1.0),cont = c(0.0,0.4)){
    y = Sys.time()
    data = data[unique(names(data))]
    cat = subset(catalogue,GI.Vec %in% names(data))
    nms = names(data) %in% cat$GI.Vec
    IDs = names(data)[nms]
    tmp1 = list()
    isUsed = c()
    
    n = 1
    for(i in 1:length(cat$comb)){
        if(!(cat$comb[i] %in% isUsed)){
            tmp1[[n]] = cat$GI.Vec[which(cat$comb == cat$comb[i])]
            isUsed[n] = cat$comb[i]
            n = n+1
        }
    }
    
    lengths = list()
    lengthSums = c()
    accession = list()
    n = 1
    for(i in 1:length(tmp1)){
        tmp2 = c()
        tmp5 = c()
        lengthSums[i] = 0
        n = 1
        for(j in 1:length(tmp1[[i]])){
            tmp5[j] = which(IDs == tmp1[[i]][j]) -1
            tmp2[n] = as.numeric(tmp1[[i]][j])
            tmp2[n+1] = sum(data[[tmp1[[i]][j]]]$GENOME[[1]]@lengths)
            lengthSums[i] = lengthSums[i] + tmp2[length(tmp2)]
            n = n +2
        }
        lengths[[i]] = tmp2
        accession[[i]] = tmp5
    } 
    
    pifams = list()
    Orfs = list()
    m = 1
    for(i in 1:length(data)){
        if(nms[i]){
            tmp3 = list()
            tmp4 = list()
            n = 1
            for(j in 1:length(data[[i]]$GENOME)){
                if((!is.null(data[[i]]$GENOME[[j]])) && (!is.null(data[[i]]$ORF[[j]]))){
                    tmp3[[n]] = data[[i]]$GENOME[[j]]@lengths
                    tmp3[[n+1]] = data[[i]]$GENOME[[j]]@values
                    tmp4[[n]] = data[[i]]$ORF[[j]]@lengths
                    tmp4[[n+1]] = data[[i]]$ORF[[j]]@values
                    n = n +2
                }
            }
            pifams[[m]] = tmp3
            Orfs[[m]] = tmp4
            m = m +1
        }
    }
    
    if(length(lengths) > 1){
        print(Sys.time() -y)
        x = Sys.time()
        res = compTestData(pifams,Orfs,lengths,lengthSums,minContigLength,meanContigLength,number,comp,cont,accession,seed,distr)
    }
    else{
        res = list()
        print("zu wenig Genome in dem Daten Satz")
    }
    
    print(Sys.time() -x)
    print(Sys.time() -y)
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