#' @export
getContigsAsIRanges <- function(data,catalogue,minContigLength,meanContigLength,comp = c(0.6,1),cont = c(0,0.4),number,seed,distr = "normal"){
    
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
    
    contigs = mkContigs(lengths,IDs,lengthSums,minContigLength,meanContigLength,number,comp,cont,seed,distr)
    contranges = list()
    print(length(contigs))
    
    for(i in 1:length(contigs)){
        j = 1
        p = 1
        tmpc = list()
        comp = list()
        cont = list()
        for(n in seq(1,(length(contigs[[i]][[1]])-2),3)){
            comp[[j]] = IRanges(start = contigs[[i]][[1]][[n+1]],end = contigs[[i]][[1]][[n+2]],names = rep(contigs[[i]][[1]][[n]],length(contigs[[i]][[1]][[n+1]])))
            j = j +1
        }
        comp[[j]] = contigs[[i]][[1]][[n+3]]/contigs[[i]][[1]][[n+4]]
        tmpc[[p]] = comp
        p = p +1
        tmpNms = c("completeness")
        if(length(contigs[[i]]) > 1){
            j = 1
            for(n in seq(1,(length(contigs[[i]][[2]])-2),3)){
                cont[[j]] = IRanges(start = contigs[[i]][[2]][[n+1]],end = contigs[[i]][[2]][[n+2]],names = rep(contigs[[i]][[2]][[n]],length(contigs[[i]][[2]][[n+1]])))
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



pfamCounter <- function(pfams,catalogue,minContigLength,meanContigLength,number,comp = c(0.6,1),cont = c(0,0.4),seed =1){
    x = Sys.time()
    res = list()
    
    contigs = getContigsAsIRanges(pfams,catalogue,minContigLength ,meanContigLength,comp,cont,number,seed)
    print(length(contigs))
    print("start counting")
    
    for(i in 1:length(contigs)){    # fuer alle bins
            bin = list()
            print(i)
            for(j in 1:(length(contigs[[i]]))){       # fuer completeness und contamination
                compCont = list()
                chromID =c()
                pfamName =c()
                pCount = c()
                bCount = c()
                nms = list()

                
                for(n in 1:(length(contigs[[i]][[j]]) -1)){
                    
                    chromID[n] = names(contigs[[i]][[j]][[n]])[1]
                    
                    this = pfams[[chromID[n]]]
                    
                    nms[[n]] = this$pfam
                }
                
                nms = unique(unlist(nms))
                dat = data.table(pfam = as.integer(nms),pfamCount = rep(0,length(nms)),baseCount = rep(0,length(nms)),key = "pfam")

                for(n in 1:(length(contigs[[i]][[j]]) -1)){     #fuer alle chromosomen
                    
                    chromID[n] = names(contigs[[i]][[j]][[n]])[1]
                    
                    this = pfams[[chromID[n]]]
                    pfamRange = IRanges(start = this$start ,end = this$end,names = this$pfam)

                    ol = findOverlaps(pfamRange,contigs[[i]][[j]][[n]])
                    overlaps = overlapsRanges(pfamRange,contigs[[i]][[j]][[n]],hits = ol)
                    this = this[queryHits(ol)]
                    dt = data.table(pfam = as.integer(this$pfam), width = width(overlaps),contigORFComb = paste(this$ORF,subjectHits(ol)),key = "pfam")
                    
                    bCount = dt[,sum(width),pfam]
                    pCount = dt[,uniqueN(contigORFComb),pfam]

                    dat[pfam %in% bCount$pfam,baseCount := baseCount + bCount$V1]
                    dat[pfam %in% pCount$pfam,pfamCount := pfamCount + pCount$V1]


                    
                }   #ende alle chromosomen
                    
                last = length(contigs[[i]][[j]])
                
                compCont= list(chromID,dat[pfamCount > 0]$pfam,dat[pfamCount > 0]$pfamCount,dat[pfamCount > 0]$baseCount,contigs[[i]][[j]][[last]])
                if(j == 1){
                    names(compCont) = c("chromID","pfamName","pfamCount","baseCount","completeness")
                }
                else{
                    names(compCont) = c("chromID","pfamName","pfamCount","baseCount","contamination")
                }
                    
                bin[[j]] = compCont
            
            }       # ende  contamination und completeness
            if(length(contigs[[i]]) > 1){
                names(bin) = c("completeness","contamination")
            }
            else{
                names(bin) = c("completeness")
            }
            res[[i]] = bin
    }       # ende alle bins
    print("end counting")
    print(Sys.time() -x)
    return(res)
}




