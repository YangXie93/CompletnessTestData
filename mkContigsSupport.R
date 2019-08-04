timeMkContigs <- function(times,sameCLen,y){
    t = Sys.time()
    for(i in 1:times){
        x = mkContigs(lengths = y,minContigLength = 1000,meanContigLength = 2000,sameContigLength = sameCLen)
    }
    print(Sys.time() -t)
}


y = list()
l = c()
for(i in 1:length(con)){
    y = append(y,sum(con[[i]]$GENOME[[1]]@lengths))
    l[i] = sum(con[[i]]$GENOME[[1]]@lengths)
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


testBinComplietnessTester <- function(bins){
    complete = list()
    for(i in length(bins)){
        tmp = list()
        for(j in 1:length(bins[[i]])){
            for(k in 1:length(bins[[i]][[j]])){
                
            }
        }
        complete[[i]] = tmp
    }
}

y = list()
z = c()
n = 0
for(i in seq(1,length(con),2)){
    y[[n]] = c(sum(con[[i]]$GENOME[[1]]@lengths),sum(con[[i+1]]$GENOME[[1]]@lengths))
    z[n] = sum(con[[i]]$GENOME[[1]]@lengths) + sum(con[[i+1]]$GENOME[[1]]@lengths)
    n = n +1
}

#-------------- count pifams -----------------------------------


listA = list(list(),list())
listA[[1]][[1]] = c(10,15,15,10,10)
listA[[1]][[2]] = c(-1,3,-1,5,-1)
listA[[2]][[3]] = c(10,15,15,10,10)
listA[[2]][[4]] = c(-1,3,-1,5,-1)

listB = list(list(),list())
listB[[1]][[1]] = c(0)
listB[[1]][[2]] = c(20,40)
listB[[1]][[3]] = c(35,60)
listB[[2]][[1]] = c(1)
listB[[2]][[2]] = c(20,40)
listB[[2]][[3]] = c(35,60)

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
    conts = list()
    for(j in 1:1){
        conts[[j]] = i +j -1
        conts[[j+1]] = c(1)
        conts[[j+2]] = c(sum(con[[i +j -1]]$GENOME[[1]]@lengths))
    }
    contigs1[[length(contigs1)+1]] = conts
}

countPifams(pifams,contigs1)

timeCountPifams <- function(listA,listB,times){
    x = Sys.time()
    for(i in 1:times){
        y = countPifams(listA,listB)
    }
    print(Sys.time() -x)
}

#---------- test contamination --------------------------

res = vector(mode = "integer")
x = c(1000,100,100,100)
contamination(x,300,0,res)
