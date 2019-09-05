#' @useDynLib CompletenessTestData
#' @importFrom Rcpp sourceCpp
#' @export
completenessTestData <- function(data,catalogue,minContigLength,meanContigLength,number,seed = 0,distr = "normal",comp = c(0.6,1.0),cont = c(0.0,0.4)){
    x = Sys.time()
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
    unqs = list()
    for(i in 1:length(data)){
        if(nms[i]){
            tmp3 = list()
            tmp4 = list()
            tmp6 = c()
            for(j in 1:length(data[[i]]$GENOME)){
                if((!is.null(data[[i]]$GENOME[[j]])) && (!is.null(data[[i]]$ORF[[j]]))){
                    tmp3[[length(tmp3)+1]] = data[[i]]$GENOME[[j]]@lengths
                    tmp3[[length(tmp3)+1]] = data[[i]]$GENOME[[j]]@values
                    tmp4[[length(tmp4)+1]] = data[[i]]$ORF[[j]]@lengths
                    tmp4[[length(tmp4)+1]] = data[[i]]$ORF[[j]]@values
                    tmp6 = append(tmp6,data[[i]]$GENOME[[j]]@values)
                }
            }
            pifams[[i]] = tmp3
            Orfs[[i]] = tmp4
            tmp6 = sort(unique(tmp6))
            unqs[[i]] = tmp6[2:length(tmp6)]
            typeof(tmp6)
        }
    }
    if(length(lengths) > 1){
        print(Sys.time() -x)
        res = compTestData(pifams,Orfs,lengths,lengthSums,minContigLength,meanContigLength,number,comp,cont,accession,unqs,seed,distr)
    }
    else{
        res = list()
        print("zu wenig Genome in dem Daten Satz")
    }
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