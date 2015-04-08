#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("depmixS4")
#install.packages("foreach")
#install.packages("doParallel")
#install.packages("parallel")

setwd("/home/mitchell/IW/hmm")
library(depmixS4)
library(ggplot2)
library(reshape2)
library(gtools)
library("parallel")
library("foreach")
library("doParallel")

hmmfiles = list.files("./freebayes", pattern="*.hmm")
plots = list()
input = list()
for( i in 1:length(hmmfiles)){
  file <- paste0("freebayes/", hmmfiles[i])
  input[[i]] <- as.data.frame( read.table( file ) )
  colnames(input[[i]]) <- c("chr", "pos", "state", "ao", "dp", "frac")
  p <-  ggplot(input[[i]], aes(x=chr, y=pos, color=state)) + geom_line(size=2)
  p <- p + scale_colour_gradientn(colours=c("blue","red")) 
  plots[[i]] <-p
}

# plotshmm = list()
# for( i in 1:length(hmmfiles)){
#   mod <- depmix(state ~ 1,data=input[[i]],nstates=2,family=binomial())
#   done <- fit(mod, verbose = FALSE)
#   probs <- posterior(done)  
#   input[[i]]["hmm"] <- probs[1] - 1
#   p <-  ggplot(input[[i]], aes(x=chr, y=pos, color=hmm)) + geom_line(size=2)
#   p <- p + scale_colour_gradientn(colours=c("blue","red","green")) 
#   plotshmm[[i]] <-p
# }
# 
# save(plotshmm, file="plotshmm.data")
# plotshmm[[1]]






# multiple runs
runhmm <- function(){
  setwd("/home/mitchell/IW/hmm")
  hmmfiles = list.files("./freebayes", pattern="*.hmm")
  input = list()
  for( i in 1:length(hmmfiles)){
    file <- paste0("freebayes/", hmmfiles[i])
    input[[i]] <- as.data.frame( read.table( file ) )
    colnames(input[[i]]) <- c("chr", "pos", "state", "ao", "dp", "frac")
  }
  
  for( i in 1:length(hmmfiles)){
    mod <- depmix(state ~ 1,data=input[[i]],nstates=2,family=binomial())
    done <- fit(mod, verbose = FALSE)
    probs <- posterior(done)  
    input[[i]]["hmm"] <- probs[1] - 1
  }
  
  num = 1:length(hmmfiles)
  diffs <- list()
  for(j in 1:length(hmmfiles)){
    i = 1
    num[j] <- 0
    while( i <= length(input[[j]]$hmm)){
      if( input[[j]]$hmm[i] == 1 ){    
        num[j] = num[j] + 1  
        pre <- input[[j]]$chr[i]
        prepos <- input[[j]]$pos[i]
        while( input[[j]]$hmm[i] == 1 & input[[j]]$chr[i] == pre ){
          i = i + 1
          if( i > length(input[[j]]$hmm)){
              break
          }
        }
        postpos <- input[[j]]$pos[i]
        diffs <- c(diffs, postpos - prepos) 
      }
      else{
        i = i + 1
      }  
    }
  }
  counts <- data.frame(hmmfiles, num)
  counts <- counts[ mixedsort(counts$hmmfiles),]
  return( counts$num )
}


cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)
data <- foreach(i = 1:2, .packages = c("depmixS4", "gtools"), .combine = cbind) %dopar% {
                 try({
                   x <- runhmm()
                   x
                 })
}

rownames(data) <- mixedsort(hmmfiles)
data <- as.data.frame(data)
mm <- as.matrix(data)

v = list()
m = list()
rtn = list()
for(i in 1:17){
    v[i] <- var(mm[i,])
    m[i] <- mean(mm[i,])
    temp <- table(as.vector(mm[i,]))
    rtn[i] <- names(temp)[temp == max(temp)][1]
}

data$mean <- m
data$var <- v
data$mode <- as.numeric(rtn)
save(data, file="hmmcounts.data")
stopCluster(cl)



