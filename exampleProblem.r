#install.packages("devtools")  # so we can install from github
#library(devtools)
#install_github("ropensci/plotly")  # plotly is part of ropensci

library(ggplot2)
library(reshape2)
library(gtools)
library(scales)
library(dplyr)

library(devtools)
library(plotly)



setwd("/home/mitchell/IW/hmm")

load("mutations.data")
head(mutations)

p <- ggplot(mutations, aes(x = strain, fill=Mutation ))
p <- p + geom_bar(stat="identity", position="fill", aes(y=Percentage ))
p <- p + scale_y_continuous(labels = percent_format())

p

py <- plotly(username="mrvollger", key="aboy5gyo0b")
r<-py$ggplotly(p)
r$response$url
