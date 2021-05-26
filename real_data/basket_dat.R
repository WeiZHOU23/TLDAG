
#########  Basketball data ############

install.packages("SportsAnalytics")
library(SportsAnalytics)
data("NBAPlayerStatistics0910")
mydata <- NBAPlayerStatistics0910

## data processing
deleteCol <- c(1,2,3,4,5,23,24)
mydataPro <- mydata[, -deleteCol]
mydataMat <- as.matrix(mydataPro)

## lambda chosen for the tuning parameter in glmnet
grid.tc <- 10^(-2+0.15*seq(0,60,by=1))

## epsilon_t for reconstructing the layer 
grid.t <- grid.tc[grid.tc < 1]



getwd()
source("tldag_basket.r")
## the estimated result
EDAG(mydataMat, grid.t, 1, 0, tuning = "adaptive", c0=2.5, lambda.index = 27)



