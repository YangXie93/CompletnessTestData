getContigsAsIRanges <- function(data,catalogue,minContigLength,meanContigLength,comp = c(0.6,1),cont = c(0,0.4),number,seed,distr = "normal",debugInfo = FALSE){
    
    data = data[unique(names(data))]
    cat = subset(catalogue,GI.Vec %in% names(data))
    nms = names(data) %in% cat$GI.Vec
    
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
    IDs = list()
    for(i in 1:length(tmp1)){
        tmp2 = c()
        tmp3 = c()
        lengthSums[i] = 0
        for(j in 1:length(tmp1[[i]])){
            tmp3[j] = as.numeric(tmp1[[i]][j])
            tmp2[j] = max(data[tmp1[[i]][j]][[1]]$end)
            lengthSums[i] = lengthSums[i] + tmp2[length(tmp2)]
        }
        lengths[[i]] = tmp2
        IDs[[i]] = tmp3
        
    } 
    
    contigs = mkContigs(lengths,IDs,lengthSums,minContigLength,meanContigLength,number,comp,cont,seed,distr,debugInfo)
    contranges = list()
    
    for(i in 1:length(contigs)){
        j = 1
        p = 1
        tmpc = list()
        comp = list()
        cont = list()
        
        scndLayr = list()
        for(n in seq(1,(length(contigs[[i]][[1]])-2),3)){
            for(m in 1:length(contigs[[i]][[1]][[n+1]])){
                scndLayr[[m]] = IRanges(start = contigs[[i]][[1]][[n+1]][m],end = contigs[[i]][[1]][[n+2]][m],names = contigs[[i]][[1]][[n]])
                
            }
            comp[[j]] = scndLayr
            j = j +1
        }
        
        
        comp[[j]] = contigs[[i]][[1]][[n+3]]/contigs[[i]][[1]][[n+4]]
        tmpc[[p]] = comp
        p = p +1
        tmpNms = c("completeness")
        if(length(contigs[[i]]) > 1){
            j = 1
            scndLayr = list()
            for(n in seq(1,(length(contigs[[i]][[2]])-2),3)){
                for(m in 1:length(contigs[[i]][[2]][[n+1]])){
                    scndLayr[[m]] = IRanges(start = contigs[[i]][[2]][[n+1]][m],end = contigs[[i]][[2]][[n+2]][m],names = contigs[[i]][[2]][[n]])
                }
                cont[[j]] = scndLayr
                j = j +1
            }
            cont[[j]] = contigs[[i]][[2]][[n +3]]/contigs[[i]][[2]][[n+4]]
            tmpc[[p]] = cont
            tmpNms[[p]] = "contamination"
            p = p+1
            
        }
        names(tmpc) = tmpNms
        contranges[[i]] = tmpc
    }
    
    return(contranges)
    
}