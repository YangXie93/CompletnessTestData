getContigsAsIRanges3 <- function(data,catalogue,minContigLength,meanContigLength,comp = c(0.6,1),cont = c(0,0.4),number,seed,distr = "normal",debugInfo = FALSE){
    
    library("data.table")
    library("IRanges")
    
    data = data[unique(names(data))]
    cat = subset(catalogue,GI.Vec %in% names(data))
    nms = names(data) %in% cat$GI.Vec
    
    cat.tmp = cat[, .(MaxVal = max(as.numeric(length.Vec))), by = comb]
    cat.dt = cat[,.(LV = unique(length.Vec)),by = .(comb,GI.Vec)]
    
    IDs = split(as.numeric(cat.dt$GI),cat.dt$comb)
    names(IDs) = NULL
    
    lengths = split(as.numeric(cat.dt$LV),cat.dt$comb)
    names(lengths) = NULL
    
    lengthSums = sapply(1:length(lengths), function(x) sum(lengths[[x]]))
    
    
    cValVec = sapply(1:length(IDs), function(x) min(IDs[[x]]))
    SortIdx = order(cValVec)
    IDs = IDs[SortIdx]
    lengths = lengths[SortIdx]
    lengthSums = lengthSums[SortIdx]
    
    contigs = mkContigs(lengths,IDs,lengthSums,minContigLength,meanContigLength,number,comp,cont,seed,distr,debugInfo)
    contranges = list()
    
    
    for(i in 1:length(contigs)){
        
        q = 1
        p = 1
        tmpc = list()
        comp = list()
        cont = list()
        for(j in 1:length(contigs[[i]])){
            for(n in seq(1,(length(contigs[[i]][[j]])-2),3)){
                for(m in 1:length(contigs[[i]][[j]][[n+1]])){
                    comp[[1]] = list(IRanges(start = contigs[[i]][[j]][[n+1]][m],end = contigs[[i]][[j]][[n+2]][m]),contigs[[i]][[j]][[n+3]]/contigs[[i]][[j]][[n+4]],contigs[[i]][[j]][[n]])
                    contranges[[q]] = comp
                    q = q +1

                }
            }
 
        }
        
    }

    return(contranges)
    
}



pfamsPerContig <- function(pfams,catalogue,minContigLength,meanContigLength,comp = c(0.6,1),cont = c(0,0.4),seed =1){
    x = Sys.time()
    res = list()
    
    contigs = getContigsAsIRanges3(pfams,catalogue,minContigLength ,meanContigLength,comp,cont,1,seed)
    print(head(contigs))
    print("start overlap detection")
    
    xL = list(); yL = list(); zL = list(); nL = list(); iL = list()
    
    for (i in 1:2){
        xL[[i]] = list();      yL[[i]] = list();      zL[[i]] = list();      nL[[i]] = list();      iL[[i]] = list()
    }
    
    #do all Completness
    j = 1
    LL = contigs.to.list(contigs,pfams,j)
    str(LL)
    xL[[j]] = LL[[1]]; yL[[j]] = LL[[2]]; zL[[j]] = LL[[3]]; nL[[j]] = LL[[4]]; iL[[j]] = LL[[5]];
    j = 2
    LL = contigs.to.list(contigs,pfams,j)
    xL[[j]] = LL[[1]]; yL[[j]] = LL[[2]]; zL[[j]] = LL[[3]]; nL[[j]] = LL[[4]]; iL[[j]] = LL[[5]];
    
    
    
    print("end overlap detection")
    print(Sys.time() -x)
    
    print("building dt")
    comp.dt = data.table(x = unlist(xL[[1]]), y = unlist(yL[[1]]), z = unlist(zL[[1]]), n = unlist(nL[[1]]), i = unlist(iL[[1]]))
    cont.dt = data.table(x = unlist(xL[[2]]), y = unlist(yL[[2]]), z = unlist(zL[[2]]), n = unlist(nL[[2]]), i = unlist(iL[[2]]))
    
    
    print("end building dt")
    print(Sys.time() -x)
    
    print("start counting")
    
    print(comp.dt)
    print(cont.dt)
    
    comp.cdt = comp.dt[,.(N = .N, s = sum(z)), by = .(x,y,n,i)]
    # cont.cdt = cont.dt[,.(N = .N, s = sum(z)), by = .(x,y,n,i)]
    comp.cdt2 = comp.cdt[,.(pCount = .N, pLen = sum(s)),by = .(i,y)]
    # cont.cdt2 = cont.cdt[,.(pCount = .N, pLen = sum(s)),by = .(i,y)]

    print(comp.cdt)
    # print(cont.cdt)
    
    setkey(comp.cdt2,i)
    # setkey(cont.cdt2,i)
    
    
    #return(comp.cdt)
    
    #print(comp.cdt2)
    
    print("end counting")
    print(Sys.time() -x)
    
    print("start reformating")
    L.Comp = cdt.to.list(comp.cdt2)
    # L.Cont = cdt.to.list(cont.cdt2)
    
    AAV =
        sapply(1:length(contigs), function(i) {
            gi = as.character(contigs[[i]][[1]][[length(contigs[[i]][[1]])]][[1]])
            cat(gi,'\n')
            return(gi[1])
            #catalogue[GI.Vec == gi,comb]
        })
    
    ABV =
        sapply(1:length(contigs), function(i) {
            return(contigs[[i]][[1]][[(length(contigs[[i]][[1]])-1)]])
        })
    
    ACV =
        sapply(1:length(contigs), function(i) {
            if (length(contigs[[i]]) == 1){
                gi = as.character(contigs[[i]][[1]][[length(contigs[[i]][[1]])]][[1]])
                return(gi)
                #return(catalogue[GI.Vec == gi,comb])
            }
            
            return(NA)
        })
    
    ADV =
        sapply(1:length(contigs), function(i) {
            if (length(contigs[[i]]) == 2){
                return(contigs[[i]][[2]][[(length(contigs[[i]][[2]])-1)]])
            }
            
            return(0)
        })
    print(AAV)
    AAV = match(AAV,catalogue$GI.Vec)
    
    print(AAV)
    
    AAV = catalogue$comb[AAV]
    AAV[is.na(AAV)] = -1
    
    # ACV = match(ACV,catalogue$GI.Vec)
    # ACV = catalogue$comb[ACV]
    ACV[is.na(ACV)] = -1
    
    print("end reformating")
    print(Sys.time() -x)
    
    print(Sys.time() -x)
    return(list(values = L.Comp,compComb.Nr = AAV, comp = ABV,comp.GI = ACV,cont = ADV))
}



obtain.full.genome.as.dt <- function(some.list){
    for (i in 1:6){
        GVec = some.list$GENOME[[i]]@values
        GEVec = cumsum(some.list$GENOME[[i]]@lengths)
        GSVec = c(1,(GEVec+1)[1:(length(GEVec)-1)])
        IG = GVec != -1
        GVec = GVec[IG]
        GEVec = GEVec[IG]
        GSVec = GSVec[IG]
        
        OVec = some.list$ORF[[i]]@values
        OEVec = cumsum(some.list$ORF[[i]]@lengths)
        OSVec = c(1,(OEVec+1)[1:(length(OEVec)-1)])
        IO = OVec != -1
        OVec = OVec[IO]
        OEVec = OEVec[IO]
        OSVec = OSVec[IO]
        
        irG = IRanges(start = GSVec, end = GEVec, names = GVec)
        irO = IRanges(start = OSVec, end = OEVec, names = OVec)
        ol = findOverlaps(irO,irG)
        dt.tmp = data.table(start = GSVec[ol@to], end = GEVec[ol@to], pfam = GVec[ol@to], ORF = OVec[ol@from])
        
        
        #x = OVec[ol@from], y = GVec[ol@to], s = GSVec[ol@to], e = GEVec[ol@to])
        
        if (i == 1){
            dt.tmp2 = dt.tmp
        }
        
        else{
            dt.tmp2 = rbindlist(list(dt.tmp2,dt.tmp))
            
        }
    }
    return(dt.tmp2)
}

contigs.to.list <- function(contigs,pfams,j){
    kkk = 1
    
    xL = list(); yL = list(); zL = list(); nL = list(); iL = list()
    
    for (i in 1:2){
        xL[[i]] = list();      yL[[i]] = list();      zL[[i]] = list();      nL[[i]] = list();      iL[[i]] = list()
    }
    print(j)
    if (j == 1){
        for(i in 1:length(contigs)){
            chromID = as.character(contigs[[i]][[j]][[length(contigs[[i]][[j]])]])
            print(chromID)
            for(n in 1:(length(contigs[[i]][[j]]) -2)){
                this = pfams[[chromID[n]]]
                pfamRange = IRanges(start = this$start, end = this$end)
                
                ol = findOverlaps(contigs[[i]][[j]][[n]], pfamRange)
                overlaps = overlapsRanges(contigs[[i]][[j]][[n]], pfamRange,hits = ol)

                
                idx = subjectHits(ol)
                
                ORF.ids = this$ORF[idx]
                Pfam.ids = this$pfam[idx]
                
                xL[[j]][[kkk]] = ORF.ids
                yL[[j]][[kkk]] = Pfam.ids
                zL[[j]][[kkk]] = width(overlaps)
                nL[[j]][[kkk]] = rep(n,length(idx))
                iL[[j]][[kkk]] = rep(i,length(idx))
                kkk = kkk + 1
            }
        }
    }
    
    else{
        for(i in 1:length(contigs)){
            if (length(contigs[[i]]) == 2){
                chromID = as.character(contigs[[i]][[j]][[length(contigs[[i]][[j]])]])
                print(chromID)
                for(n in 1:(length(contigs[[i]][[j]]) -2)){
                    this = pfams[[chromID[n]]]
                    pfamRange = IRanges(start = this$start, end = this$end)
                    
                    ol = findOverlaps(contigs[[i]][[j]][[n]], pfamRange)
                    overlaps = overlapsRanges(contigs[[i]][[j]][[n]], pfamRange,hits = ol)
                    
                    idx = subjectHits(ol)
                    
                    ORF.ids = this$ORF[idx]
                    Pfam.ids = this$pfam[idx]
                    
                    xL[[j]][[kkk]] = ORF.ids
                    yL[[j]][[kkk]] = Pfam.ids
                    zL[[j]][[kkk]] = width(overlaps)
                    nL[[j]][[kkk]] = rep(n,length(idx))
                    iL[[j]][[kkk]] = rep(i,length(idx))
                    kkk = kkk + 1
                }
            }
            
        }
    }
    return(list(xL,yL,zL,nL,iL))
}

cdt.to.list <- function(cdt){
    diffVec = diff(cdt$i)
    EndsVec = which(diffVec != 0)
    StartsVec = c(1,EndsVec+1)
    EndsVec = c(EndsVec,length(cdt$i))
    
    maxi = max(cdt$i)
    
    L <- vector(mode = 'list', length = maxi)
    
    for (iter in 1:length(StartsVec)){
        jter = cdt$i[ StartsVec[iter] ]
        L[[jter]] <- list()
        L[[jter]]$pID = cdt$y[ StartsVec[iter] : EndsVec[iter] ]
        L[[jter]]$pCount = cdt$pCount[ StartsVec[iter] : EndsVec[iter] ]
        L[[jter]]$pLen = cdt$pLen[ StartsVec[iter] : EndsVec[iter] ]
    }
    
    return(L)
}