getSingleGemomeContigsAsIRanges <- function(lengths,minContigLength,meanContigLength,comp = 0.1,seed,distr = "normal",debugInfo = FALSE){
    

    IDs = c(1:length(lengths))
    
    contigs = singleGenomeMkContigs(lengths,IDs,comp,minContigLength,meanContigLength,seed,distr,debugInfo)
    contranges = list()

    n = 1
    j = 2
    p = 3
    for(i in 1:((length(contigs[[1]])-2)/3)){
        contranges[[i]] = IRanges(start = contigs[[1]][[j]],end = contigs[[1]][[p]],names = rep(contigs[[1]][[n]],length(contigs[[1]][[p]])))
        p = p +3
        j = j +3
        n = n + 3
    }
    
    print(paste("tatsaechliche comp: ",contigs[[1]][[length(contigs[[1]])-1]]/contigs[[1]][[length(contigs[[1]])]]))
    
    return(contranges)
    
}