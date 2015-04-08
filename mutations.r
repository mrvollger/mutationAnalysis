#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("depmixS4")
#install.packages("devtools")  # so we can install from github
#library(devtools)
#install_github("ropensci/plotly")  # plotly is part of ropensci

library(depmixS4)
library(ggplot2)
library(reshape2)
library(gtools)
library(scales)
library(dplyr)
library(devtools)
library(plotly)

setwd("/home/mitchell/IW/hmm")

mutationfiles =list.files("./mutations", pattern="*.mutations")
plots = list()
input = list()
for( i in 1:length(mutationfiles)){
  file <- paste0("mutations/", mutationfiles[i])
  input[[i]] <- as.data.frame( read.table( file ) )
  colnames(input[[i]]) <- c("chr", "pos", "change", "ref", "alt", "dp", "ao")
  strain <- mutationfiles[i]
  strain <- substr(strain, 0, nchar(strain)-10)
  input[[i]]$strain <- rep(strain, length(input[[i]]$change))
  input[[i]]$Percentage <- rep(1, length(input[[i]]$change))
  
  #A:T
  #G:C
  fix <- gsub("A>C|T>G", "A>C:T>G", input[[i]]$change)
  fix <- gsub("A>G|T>C", "A>G:T>C", fix)
  fix <- gsub("A>T|T>A", "A>T:T>A", fix)
  fix <- gsub("G>A|C>T", "G>A:C>T", fix)
  fix <- gsub("G>T|C>A", "G>T:C>A", fix)
  fix <- gsub("G>C|C>G", "G>C:C>G", fix)
  input[[i]]$Mutation <- fix
  
  input[[i]] <- input[[i]][order(input[[i]]$Mutation),]
  
}

mutations <- rbind(input[[1]], input[[2]], input[[3]], input[[4]], input[[5]], input[[6]])

#remove systematic mutaiton calls
mutations$id <- with(mutations, paste0(chr, pos, Mutation))
index1 <- duplicated(mutations$id) 
index2 <- duplicated(mutations$id, fromLast=TRUE) 
index <- (index1+index2)==0
mutations<- mutations[index,]

# make cis and uv data frames to mkae it easer to compare rad14 changes
cis <- mutations[grep("cis*", mutations$strain),]
uv <- mutations[grep("uv*", mutations$strain),]

pcount <- qplot(strain, data=mutations, geom="bar", fill=Mutation)
pcountcis <-qplot(strain, data=cis, geom="bar", fill=Mutation)
pcountuv <- qplot(strain, data=uv, geom="bar", fill=Mutation)

pcis <- ggplot(cis, aes(x = strain, fill=Mutation ))
pcis <- pcis + geom_bar(stat="identity", position="fill", aes(y=Percentage ))
pcis <- pcis + scale_y_continuous(labels = percent_format())

puv <- ggplot(uv, aes(x = strain, fill=Mutation ))
puv <- puv + geom_bar(stat="identity", position="fill", aes(y=Percentage ))
puv <- puv + scale_y_continuous(labels = percent_format())

p <- ggplot(mutations, aes(x = strain, fill=Mutation ))
p <- p + geom_bar(stat="identity", position="fill", aes(y=Percentage ))
p <- p + scale_y_continuous(labels = percent_format())

pcount
pcountcis
pcountuv

p
pcis
puv

save(mutations, file="mutations.data")


