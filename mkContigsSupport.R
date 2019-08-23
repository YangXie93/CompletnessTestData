con = readRDS("~/Downloads/IDz_1_0.Rds")

timeMkContigs <- function(times,y,z){
    t = Sys.time()
    for(i in 1:times){
        x = mkContigs(lengths = y,lengthSums = z,minContigLength = 1000,meanContigLength = 2000)
    }
    print(Sys.time() -t)
}

checkForZeroLength <- function(y){
    for(i in 1:length(y)){
        for(j in seq(1,length(y[[i]]),2)){
            if( length(y[[i]][[j]][[2]]) == 0){
                print(paste("bei:",i,"bei",j))
            }
        }
    }
}

y = list()
z = c()
n = 1
for(i in seq(1,length(con),2)){
    y[[length(y) +1]] = c(as.integer(names(con)[i]),sum(con[[i]]$GENOME[[1]]@lengths),as.integer(names(con)[i+1]),sum(con[[i+1]]$GENOME[[1]]@lengths))
    z[n] = sum(con[[i]]$GENOME[[1]]@lengths) + sum(con[[i+1]]$GENOME[[1]]@lengths)
    n = n +1
}

#-------------- count pifams -----------------------------------

con = readRDS("~/Downloads/IDz_1_0.Rds")

pifams = list()
ORF = list()
for(j in 1:length(con)){
    pis = list()
    orf = list()
    for(i in 1:length(con[[j]]$GENOME)){
        pis[[length(pis)+1]] = con[[j]]$GENOME[[i]]@lengths
        pis[[length(pis)+1]] = con[[j]]$GENOME[[i]]@values
        orf[[length(orf)+1]] = con[[j]]$ORF[[i]]@lengths
        orf[[length(orf)+1]] = con[[j]]$ORF[[i]]@values
    }
    pifams[[length(pifams)+1]] = pis
    ORF[[length(ORF)+1]] = orf
}



name = as.integer(names(con))

contigs1 = list()
for(i in seq(1,length(con)-1,2)){
    cons = list()
    conts = list()
    conts[[length(conts)+1]] = as.integer(name[i])
    conts[[length(conts)+1]] = c(1)
    conts[[length(conts)+1]] = c(sum(con[[i]]$GENOME[[1]]@lengths))
    conts[[length(conts)+1]] = as.integer(name[i+1])
    conts[[length(conts)+1]] = c(1)
    conts[[length(conts)+1]] = c(sum(con[[i+1]]$GENOME[[1]]@lengths))
    cons[[length(cons)+1]] = conts
    contigs1[[length(contigs1)+1]] = cons
}

names = c()
for(i in 1:length(con)){
    names[i] = i
}


countPifams(pifams,contigs1)

timeCountPifams <- function(listA,listB,names,times){
    x = Sys.time()
    for(i in 1:times){
        y = countPifams(listA,listB,names)
    }
    print(Sys.time() -x)
}

#---------- test completenessTestData --------------------------



testRes <- function(data,res){
    link = as.integer(names(data))
    r = c()
    tmp = c()
    for(i in 1:length(res)){
        lnk = which(link == res[[i]][[1]]$compChromID)
        for(j in 1:length(data[[lnk]]$GENOME)){
            tmp = append(tmp,data[[lnk]]$GENOME[[j]]@values)
        }
        print(length(res[[i]][[1]]$compPifamNames) < (length(unique(tmp)) -1))
        r[i] = (length(res[[i]][[1]]$compPifamNames) == (length(unique(tmp)) -1))
    }
    return(r)
}




for(j in 1:length(data[[1]]$GENOME)){
    tmp = append(tmp,data[[1]]$GENOME[[j]]@values)
}


zero = readRDS("~/Daten/IDz_Zero_Zero.Rds")


library(data.table)
cat = data.table(GI.Vec = names(con),comb = c(1:length(con)))

x = completenessTestData(zero,cat,10000,20000,1000)


#------------------------------------------------------

l = c(1,20,2,10,3,15,4,12)
t = 47


testForToHighCont <- function(x){
    for(i in 1:length(x)){
        if(x[[i]][[2]]$contamination > 0.4){
            print(paste("cont at:",i))
            print(x[[i]][[2]]$contamination)
        }
        if(x[[i]][[1]]$completeness < 0.6){
            print(paste("comp at: ",i))
            print(x[[i]][[1]]$completeness)
        }
    }
    
}






763166122
13401


