install.packages("referenceIntervals", dependencies=TRUE)
library(referenceIntervals)
install.packages("zoo", dependencies=TRUE)
library(zoo)




#IMPORT AND ARRANGE DATA: Example model: Age ~ Hormone level [e.g. total testosterone]
Data <- read.csv(file.choose(), header=T, sep=";", dec=",")
plot.default(Data$Age, Data$Testosterone)
y <- data.frame(Data$Age, Data$Testosterone)
colnames(y) <- c("Age","Hormone")
y <- y[complete.cases(y$Hormone), ]
y <- y[order(y$Age, decreasing = FALSE),]  
y$Age
y$Hormone


#Basic quantile function for whole column
quantile(y$Hormone, c(.5, .025, .975), na.rm=TRUE)


#Establish a series of n=120 'moving window' 95% reference intervals (p2.5 & p97.5): calculate Basic quantiles
quantiles <- c()
for(i in c(1:361)) {
  quantiles <- c(quantiles, as.numeric(quantile(y[c(seq(i,i+120,1)),], c(.5, .025, .975), na.rm=TRUE)))
  }
write.csv(quantiles, file = "quantiles.csv")
head(quantiles)


###################
#referenceIntervals package workflow
###################

#ReferenceIntervals package basic function: calculate non-parametric 95% CI (and 90% CI for both centiles) for a column
ref <- refLimit(y$Hormone, out.method = "horn", out.rm = FALSE, RI = "n", CI = "n",
                refConf = 0.95, limitConf = 0.9, bootStat = "basic")
ref
ref$Conf_Int[1:4]


#Establish a series of n=120 'moving window' 95% reference intervals: p2.5[90%CI] & p97.5[90%CI]
#Preload vectors to be filled with centile data 
LowerLim <- c()
UpperLim <- c()
CI_LL_L <- c()
CI_LL_U <- c()
CI_UL_L <- c()
CI_UL_U <- c()
test <- c()

#Centile calculation !NB heavy computation loop
for(i in seq(1, length(y$Hormone)-120, 20)) {
  
  LowerLim <- c(LowerLim, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                                        out.method = "horn", 
                                                        out.rm = F, RI = "n", 
                                                        CI = "n", refConf = 0.95, 
                                                        limitConf = 0.9, 
                                                        bootStat = "basic")$Hormone$Ref_Int[1]))
  
  UpperLim <- c(UpperLim, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                                        out.method = "horn", 
                                                        out.rm = F, RI = "n", 
                                                        CI = "n", refConf = 0.95, 
                                                        limitConf = 0.9, 
                                                        bootStat = "basic")$Hormone$Ref_Int[2]))
  
  CI_LL_L <- c(CI_LL_L, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                                        out.method = "horn", 
                                                        out.rm = F, RI = "n", 
                                                        CI = "n", refConf = 0.95, 
                                                        limitConf = 0.9, 
                                                        bootStat = "basic")$Hormone$Conf_Int[1]))
  
  CI_LL_U <- c(CI_LL_U, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                                       out.method = "horn", 
                                                       out.rm = F, RI = "n", 
                                                       CI = "n", refConf = 0.95, 
                                                       limitConf = 0.9, 
                                                       bootStat = "basic")$Hormone$Conf_Int[2]))
  
  CI_UL_L <- c(CI_UL_L, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                                       out.method = "horn", 
                                                       out.rm = F, RI = "n", 
                                                       CI = "n", refConf = 0.95, 
                                                       limitConf = 0.9, 
                                                       bootStat = "basic")$Hormone$Conf_Int[3]))
  
  CI_UL_U <- c(CI_UL_U, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                                       out.method = "horn", 
                                                       out.rm = F, RI = "n", 
                                                       CI = "n", refConf = 0.95, 
                                                       limitConf = 0.9, 
                                                       bootStat = "basic")$Hormone$Conf_Int[4]))
  
  
  
  test <- c(test,i)
}

LowerLim
UpperLim


#MEDIAN percentile 
data <- median(data)
medianRefLimit <- function(data, limitConf = 0.9){
  data <- data[order(data)]
  alpha <- 1 - limitConf
  n <- length(data)
  k <- c(qbinom(alpha/2, n, 0.5), qbinom(1-alpha/2, n, 0.5)+1)
  data[k]
}



#MEDIAN
Median <- c()
Median_L <- c()
Median_U <- c()
test2 <- c()

for(i in seq(1, length(y$Hormone)-120, 20)) {
  
  Median <- c(Median, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                              out.method = "horn", 
                                              out.rm = F, RI = "n", 
                                              CI = "n", refConf = 0, 
                                              limitConf = 0.9, 
                                              bootStat = "basic")$Hormone$Ref_Int[1]))
  
  Median_L <- c(Median_L, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                            out.method = "horn", 
                                            out.rm = F, RI = "n", 
                                            CI = "n", refConf = 0, 
                                            limitConf = 0.9, 
                                            bootStat = "basic")$Hormone$Conf_Int[1]))
  
  Median_U <- c(Median_U, as.numeric(refLimit(y[c(seq(i,i+120,1)),], 
                                            out.method = "horn", 
                                            out.rm = F, RI = "n", 
                                            CI = "n", refConf = 0, 
                                            limitConf = 0.9, 
                                            bootStat = "basic")$Hormone$Conf_Int[3]))
  
  test2 <- c(test2,i)
}

Median

#MERGE ORIGINAL DATA WITH SYNCED N=120 WINDOW CENTILE VECTORS
Data2 <- cbind(y[121:481,],LowerLim, UpperLim, CI_LL_L, CI_LL_U, CI_UL_L, CI_UL_U, Median, Median_L, Median_U)


#ORGANIZE DATAFRAME WITH CENTILE VALUES
Data3 <- data.frame(cbind(LowerLim, UpperLim, CI_LL_L, CI_LL_U, CI_UL_L, CI_UL_U, Median, Median_L, Median_U))
#DETERMINE HOW MANY WINDOWS [OF 120 OBSERVATIONS x SPACING BETEWEEN] THEM ARE BETWEEN Y OBS 481-120
(481-120)/20 # = 18.05 --> round up: 19 windows (intervals start at y obs = 120 and end at y=601)
Data4 <- Data3[1:19,]
Data4$order <- c(1:19)
#replicate Data4 to get X replications of each window
Data5 <- data.frame(Data4,i=rep(1:20,ea=NROW(Data4)))
Data6 <- Data5[order(Data5$order, decreasing = FALSE),]  
(481-120)/20  #=361 so remove tail end of centile dataset and then merge with y
Data7 <- Data6[1:362,]
Data8 <- cbind(y[120:481,], Data7)


#PLOT THE WHOLE THING
ggplot(y, aes(x=Age, y=Hormone)) + geom_point() + theme_bw() + 
  labs(title="ReferenceIntervals test", x="Age, y (girls)", y="Serum SHBG (nmol/L)") +
  geom_line(data=Data8, aes(x=Age,y=LowerLim), inherit.aes = FALSE, linetype = "solid", lwd=1.2) +
  geom_line(data=Data8, aes(x=Age,y=UpperLim), inherit.aes = FALSE, linetype = "solid", lwd=1.2) +
  geom_ribbon(data=Data8, aes(x=Age,ymin=CI_LL_L,ymax=CI_LL_U), inherit.aes = FALSE, fill = 'grey63', alpha = 0.5) +
  geom_ribbon(data=Data8, aes(x=Age,ymin=CI_UL_L,ymax=CI_UL_U), inherit.aes = FALSE, fill = 'grey63', alpha = 0.5) +
  geom_line(data=Data8, aes(x=Age,y=Median), inherit.aes = FALSE, linetype = "solid", color = 'red', lwd=1.2)


  
ggsave("ReferenceIntervalTest.pdf", dpi = 600)


#ReferenceIntervals example non-parametric
ref <- refLimit(Data2$Hormone, out.method = "horn", out.rm = FALSE, RI = "n", CI = "n",
                  refConf = NULL, limitConf = 0.9, bootStat = "basic")
ref



ReadOut <- data.frame(testVec_upper, testVec_lower, Confidence)
write.csv(testVec_upper, file = "quantiles.csv")
plot(y$Age[121:481],testVec_upper, type = "l")
+ plot(y$Age[121:481],testVec_lower, type = "l")




#ZOO PACKAGE equivalent function
t <- rollapply(y$Hormone, width = 120, by = 1, FUN = function(x) refLimit(x, out.method = "horn", 
                                                                               out.rm = F, RI = "n", 
                                                                               CI = "n", refConf = 0.95, 
                                                                               limitConf = 0.95, 
bootStat = "basic"))

write.csv(t, file = "quantiles.csv")
plot(y[119:601,]$age,t1)

plot(y$Data.Age[121:481],testVec_upper, type = "l")

Zoo <- read.csv(file.choose(), header=T, sep=";", dec=",")
plot(Zoo$LowerLim_Low)





#ALTERNATIVE APPROACH FOR BETTER SMOOTHING (separate n=120 'moving windows' by n=20)
LowerLim <- c()
UpperLim <- c()
CI_LL_L <- c()
CI_LL_U <- c()
CI_UL_L <- c()
CI_UL_U <- c()
  
  
for(i in seq(1, length(y$Hormone)-120, 20)) {
  ri  <- refLimit(y[c(seq(i,i+120,1)),], out.method = "horn", out.rm = F, RI = "n", CI = "n", refConf = 0.95, limitConf = 0.9, bootStat = "basic")

  LowerLim <- c(LowerLim, as.numeric(ri$Hormone$Ref_Int[1]))
  UpperLim <- c(UpperLim, as.numeric(ri$Hormone$Ref_Int[2]))
  CI_LL_L <- c(CI_LL_L, as.numeric(ri$Hormone$Conf_Int[1]))
  CI_LL_U <- c(CI_LL_U, as.numeric(ri$Hormone$Conf_Int[2]))
  CI_UL_L <- c(CI_UL_L, as.numeric(ri$Hormone$Conf_Int[3]))
  CI_UL_U <- c(CI_UL_U, as.numeric(ri$Hormone$Conf_Int[4]))
}

Data3 <- data.frame(cbind(LowerLim, UpperLim, CI_LL_L, CI_LL_U, CI_UL_L, CI_UL_U))

