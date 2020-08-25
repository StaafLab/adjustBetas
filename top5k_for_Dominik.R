#####======================================================================#####
### Test flexmix on top 5K CpGs  + random CpGs
#####======================================================================#####

##https://cran.r-project.org/web/packages/flexmix/vignettes/flexmix-intro.pdf
library(flexmix)

load(file="200518_betaTop5K_withPurity.RData")

ls()
#[1] "fracTum"  "testDat2" "testDat3"

str(testDat2)
# num [1:5000, 1:236] 0.008 0.044 0.163 0.705 0.65 0.067 0.057 0.704 0.06 0 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:5000] "cg06712559" "cg17928920" "cg09248054" "cg27541454" ...
#  ..$ : chr [1:236] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...

pdf("testFlexmix_100fromTop5k.pdf",width=12,height=12)
par(mfrow=c(2,1))

set.seed(12345)
for(i in sample(1:nrow(testDat2),100)) {

  y<-testDat2[i,] + c(.001,-.001)  ##flexmix throws error if one pop has e.g all zeros, add small fudgefactor
  x<-fracTum

  model <- stepFlexmix(y ~ x,k = 1:2, nrep = 7,verbose = FALSE)  ##try 1 vs 2 group solution
  model <- getModel(model, "BIC")
  #model <- flexmix(y ~ x, k = 2)#,control = list(iter.max=25)) ##straight 2-group
  plot(x, y, col = clusters(model),pch=16,main="flexmix",xlab="tum %",ylab="beta",sub=paste0("row ",i),xlim=c(0,1),ylim=c(0,1))
  if(ncol(parameters(model))>1) {
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3)
    abline(parameters(model)[1:2, 2], col = "green", lwd = 3)
  } else {
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3)
  }

  ##function for reassigning based on distance to line - ugly..
  dist2d <- function(a,b,c) {
    ya1 <- a[2]
    yb1 <- a[1]
    ya2 <- b[2]
    yb2 <- b[1]
    ##
    d<-which.min( c( (abs(yb1+ya1*c[1]-c[2])/sqrt(ya1^2+1)),
      (abs(yb2+ya2*c[1]-c[2])/sqrt(ya2^2+1)) ) )
    return(d)
  }

  if(ncol(parameters(model))>1) {
    dd<-apply(cbind(x,y),1,function(zz) { dist2d(a=parameters(model)[1:2, 1],b=parameters(model)[1:2, 2],c=zz) })
    plot(x, y, col = dd,pch=c(15,16)[1 + (dd==clusters(model))],main="reassigned d2line",,xlab="tum %",ylab="beta",xlim=c(0,1),ylim=c(0,1))
    abline(lm(y[dd==1]~x[dd==1]),col=4,lwd=3)
    abline(lm(y[dd==2]~x[dd==2]),col=3,lwd=3)
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3,lty=2)
    abline(parameters(model)[1:2, 2], col = "green", lwd = 3,lty=2)
    legend("right",legend=c("unchanged","changed"),col=1,pch=c(16,15),bty="n")
  } else {
    dd<-clusters(model)
    plot(x, y, col = dd,pch=c(15,16)[1 + (dd==clusters(model))],main="reassigned d2line",,xlab="tum %",ylab="beta",xlim=c(0,1),ylim=c(0,1))
    abline(lm(y[dd==1]~x[dd==1]),col=4,lwd=3)
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3,lty=2)
  }
  rm(x,y,model,dd)

}

dev.off()

##get just random CpGs
set.seed(20200518)
testDat3<-betaData[sample(1:nrow(betaData),5000),]

str(testDat3)
# num [1:5000, 1:236] 0.379 0.06 0.15 0.962 0.012 0.539 0.97 0.909 1 0.457 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:5000] "cg15986668" "cg10065883" "cg24015249" "cg26857712" ...
#  ..$ : chr [1:236] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...

pdf("testFlexmix_100fromRnd5k.pdf",width=12,height=12)
par(mfrow=c(2,1))

set.seed(12345)
for(i in sample(1:nrow(testDat3),100)) {

  y<-testDat3[i,] + c(.001,-.001)  ##flexmix throws error if one pop has e.g all zeros, add small fudgefactor
  x<-fracTum

  model <- stepFlexmix(y ~ x,k = 1:2, nrep = 7,verbose = FALSE)   ##try 1 vs 2 group solution
  model <- getModel(model, "BIC")
  #model <- flexmix(y ~ x, k = 2)#,control = list(iter.max=25))   ##straight 2-group
  plot(x, y, col = clusters(model),pch=16,main="flexmix",xlab="tum %",ylab="beta",sub=paste0("row ",i),xlim=c(0,1),ylim=c(0,1))
  if(ncol(parameters(model))>1) {
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3)
    abline(parameters(model)[1:2, 2], col = "green", lwd = 3)
  } else {
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3)
  }

  ##function for reassigning based on distance to line - ugly..
  dist2d <- function(a,b,c) {
    ya1 <- a[2]
    yb1 <- a[1]
    ya2 <- b[2]
    yb2 <- b[1]
    ##
    d<-which.min( c( (abs(yb1+ya1*c[1]-c[2])/sqrt(ya1^2+1)),
      (abs(yb2+ya2*c[1]-c[2])/sqrt(ya2^2+1)) ) )
    return(d)
  }

  if(ncol(parameters(model))>1) {
    dd<-apply(cbind(x,y),1,function(zz) { dist2d(a=parameters(model)[1:2, 1],b=parameters(model)[1:2, 2],c=zz) })
    plot(x, y, col = dd,pch=c(15,16)[1 + (dd==clusters(model))],main="reassigned d2line",,xlab="tum %",ylab="beta",xlim=c(0,1),ylim=c(0,1))
    abline(lm(y[dd==1]~x[dd==1]),col=4,lwd=3)
    abline(lm(y[dd==2]~x[dd==2]),col=3,lwd=3)
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3,lty=2)
    abline(parameters(model)[1:2, 2], col = "green", lwd = 3,lty=2)
    legend("right",legend=c("unchanged","changed"),col=1,pch=c(16,15),bty="n")
  } else {
    dd<-clusters(model)
    plot(x, y, col = dd,pch=c(15,16)[1 + (dd==clusters(model))],main="reassigned d2line",,xlab="tum %",ylab="beta",xlim=c(0,1),ylim=c(0,1))
    abline(lm(y[dd==1]~x[dd==1]),col=4,lwd=3)
    abline(parameters(model)[1:2, 1], col = "blue", lwd = 3,lty=2)
  }
  rm(x,y,model,dd)
}

dev.off()

#save(testDat2,testDat3,fracTum,file="200518_betaTop5K_withPurity.RData")



