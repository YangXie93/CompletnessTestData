completenessTestData <- function(data,catalogue,minContigLength,meanContigLength,number,seed = 0){
    ID = names(data)
    name = catalogue$GI.Vec
    connec = catalogue$comb
    
    tmp1 = list()
    for(i in 1:length(ID)){
        tmp = c()
        n = which(name == ID[i])
        tmp[length(tmp)+1] = ID[i]
        flag = TRUE
        while(flag){
            n = which(ID == name[connec[n]])
            if(length(n) > 0 && n != 0){
                if(name[n] != tmp[length(tmp)]){
                    tmp[length(tmp)+1] = name[n]
                    print(tmp[length(tmp)])
                }
                else{
                    flag = FALSE
                }
            }
            else{
                flag = FALSE
            }
        }
        tmp1[[length(tmp1)+1]] = tmp
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
            tmp3[[length(tmp3)+1]] = data[[i]]$GENOME[[j]]@lengths
            tmp3[[length(tmp3)+1]] = data[[i]]$GENOME[[j]]@values
            tmp4[[length(tmp4)+1]] = data[[i]]$ORF[[j]]@lengths
            tmp4[[length(tmp4)+1]] = data[[i]]$ORF[[j]]@values
        }
        pifams[[length(pifams)+1]] = tmp3
        Orfs[[length(Orfs)+1]] = tmp4
    }
    conts = mkContigs(lengths,lengthSums,minContigLength,meanContigLength,number,seed)
    res = countPifams(pifams,Orfs,conts,as.integer(ID))
    return(res)
}

timeCompletenessTestData <- function(times,data,catalogue,minContigLength,meanContigLength,seed = 0){
    time = Sys.time()

    x = completenessTestData(data,catalogue,minContigLength,meanContigLength,times,seed)

    print(Sys.time() - time)
}
