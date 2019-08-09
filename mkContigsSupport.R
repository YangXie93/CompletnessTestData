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
        for(j in 1:length(y[[i]])){
            if( length(y[[i]][[j]]) == 0){
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
for(j in 1:length(con)){
    pis = list()
    for(i in 1:length(con[[j]]$GENOME)){
        pis[[length(pis)+1]] = con[[1]]$GENOME[[i]]@lengths
        pis[[length(pis)+1]] = con[[1]]$GENOME[[i]]@values
    }
    pifams[[length(pifams)+1]] = pis
}

contigs1 = list()
for(i in 1:length(con)){
    cons = list()
    conts = list()
    for(j in 1:1){
        conts[[j]] = i +j -1
        conts[[j+1]] = c(1)
        conts[[j+2]] = c(sum(con[[i +j -1]]$GENOME[[1]]@lengths))
    }
    cons[[1]] = conts
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

cat = data.table::data.table(GI.Vec = names(con)[seq(1,length(con),2)],comb = names(con)[seq(2,length(con),2)])
