### R code from vignette source 'quantsmooth.Rnw'

###################################################
### code chunk number 1: SimplePlot
###################################################
library(quantsmooth)
data(chr14)
plot(affy.cn[,1],pch=".")
lines(quantsmooth(affy.cn[,1]),lwd=2)


###################################################
### code chunk number 2: FirstComparison
###################################################
plot(affy.pos,affy.cn[,1],ylab="copy number",xlab="position",pch=".")
lines(affy.pos,quantsmooth(affy.cn[,1]),lwd=2)
points(bac.pos,bac.cn[,1],col="red",pch=".")
lines(bac.pos,quantsmooth(bac.cn[,1]),col="red",lwd=2)
points(ill.pos,ill.cn[,1],col="blue",pch=".")
lines(ill.pos,quantsmooth(ill.cn[,1]),col="blue",lwd=2)
legend("topleft",legend=c("affymetrix","illumina","BAC"),col=c("black","red","blue"),lty=1)


###################################################
### code chunk number 3: LambdaLengthDependent
###################################################
lambda.divisor<-50
plot(affy.pos,affy.cn[,1],ylab="copy number",xlab="position",pch=".")
lines(affy.pos,quantsmooth(affy.cn[,1],smooth.lambda=length(affy.pos)/lambda.divisor),lwd=2)
points(bac.pos,bac.cn[,1],col="red",pch=".")
lines(bac.pos,quantsmooth(bac.cn[,1],smooth.lambda=length(bac.pos)/lambda.divisor),col="red",lwd=2)
points(ill.pos,ill.cn[,1],col="blue",pch=".")
lines(ill.pos,quantsmooth(ill.cn[,1],smooth.lambda=length(ill.pos)/lambda.divisor),col="blue",lwd=2)
legend("topleft",legend=c("affymetrix","illumina","BAC"),col=c("black","red","blue"),lty=1)


###################################################
### code chunk number 4: Crossvalidation
###################################################
lambdas<-2^seq(from=-2,to=5,by=0.25)
lambda.res <- rep(NA, length(lambdas))
for (lambda in 1:length(lambdas)) lambda.res[lambda] <- quantsmooth.cv(bac.cn[,1], lambdas[lambda])
plot(lambdas,lambda.res,type="l")
abline(v=lambdas[which.min(lambda.res)])


###################################################
### code chunk number 5: PercentileSmoothing
###################################################
plot(bac.pos,quantsmooth(bac.cn[,1],smooth.lambda=length(bac.pos)/lambda.divisor),col="red",type="l",lwd=2)
lines(bac.pos,quantsmooth(bac.cn[,1],smooth.lambda=length(bac.pos)/lambda.divisor,tau=0.25),col="red",lty=2)
lines(bac.pos,quantsmooth(bac.cn[,1],smooth.lambda=length(bac.pos)/lambda.divisor,tau=0.75),col="red",lty=2)


###################################################
### code chunk number 6: plotSmoothed
###################################################
plotSmoothed(bac.cn,bac.pos,ylim=c(1,2.5),normalized.to=2,smooth.lambda=length(bac.pos)/lambda.divisor)


###################################################
### code chunk number 7: ChangedRegions
###################################################
plotSmoothed(ill.cn[,1],ill.pos,ylim=c(1,2.5),normalized.to=2,smooth.lambda=length(ill.pos)/lambda.divisor)
res<-getChangedRegions(ill.cn[,1],ill.pos,normalized.to=2,interval=0.5)
segments(res[,"start"],1.0,res[,"end"],1.0,col=2,lwd=2)


###################################################
### code chunk number 8: quantsmooth.Rnw:123-124
###################################################
toLatex(sessionInfo())


