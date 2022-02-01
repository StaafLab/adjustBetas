
##define function
doLmTests<-function(y=NULL) {
  useALL<-TRUE
  ##function for angle between lines
  myFunct<-function(s1=1,s2=1) { abs( atan( (s1-s2)/(1+s1*s2) ) ) * (180/pi) }
  ##define purity vector
  x<-fracTum
  ##criteria to catch uncorrectable values
  if( sd(y) > 0 & sd(y[x >=.5]) > 0 ) {
    ##check global correlation
    corALL<-(cor(x,y))
    ##Define 2 pops and check separation metrics
    ##only define two pops if purest tumors form 2 well separated clusters
    yUse<- x >=.5
    y2<-y[yUse]
    x2<-x[yUse]
    ##ySplit<-kmeans(y[yUse],2)$cluster
     kSolutions<-lapply(1:50,function(x) {
      s1<-kmeans(y2,2)
      s2<-which.min(s1$centers)
      names(s1$cluster)[s1$cluster==s2]
    })
    s1<-unique(unlist(kSolutions))
    s1<-s1[ colSums(do.call("rbind",lapply(kSolutions,function(x) { s1 %in% x }))) >= 30]
    ySplit<-vector("integer",length=length(y2))+2
    names(ySplit)<-names(y2)
   ##criteria if not at least 10 in minor cluster
    if( length(s1) < 10 | (length(y2)-length(s1)) < 10 ) {
        mu1<-NA
        mu2<-NA
        sd1<-NA
        sd2<-NA
        aOL<-NA
        dINT<-NA
        dMEAN<-NA
        dSD<-NA
        co1<-NA
        co2<-NA
        lm1c<-c(NA,NA)
        lm2c<-c(NA,NA)
        useALL<-TRUE
        yUse1<-rep(TRUE,length(x))
        yUse2<-rep(FALSE,length(x))
    } else {
      ySplit[s1]<-1
      yUse1<-ySplit == 1
      yUse2<-ySplit == 2
      ##calculate some metrics
      mu1<-mean( y2[yUse1] )
      mu2<-mean( y2[yUse2] )
      sd1<-sd( y2[yUse1] )
      sd2<-sd( y2[yUse2] )
      med1<-median(y2[yUse1])
      med2<-median(y2[yUse2])
      pMed<-pmin( pmin(abs(1-med1),abs(0-med1)),pmin(abs(1-med2),abs(0-med2)) )
      ##get normal pop intercept
      lm1<-lm(y2[yUse1]~x2[yUse1])
      lm2<-lm(y2[yUse2]~x2[yUse2])
      lm1c<-lm1$coefficients
      lm2c<-lm2$coefficients
      ##get correlations for minor pops
      co1<-summary(lm1)$r.squared
      co2<-summary(lm2)$r.squared
      ##get diff in intercept
      dINT<-(lm1$coeff[1]-lm2$coeff[1])
      ##get diff in means
      dMEAN<-(mu2-mu1)
      ##get diff in SDs
      dSD<-mean(abs(c((mu2-mu1)/sd1,(mu1-mu2)/sd2)))
      ##new aOL -> angle between lines
      aOL<-myFunct(lm1$coeff[2],lm2$coeff[2])
#      if(  aOL > 5 & abs(dINT) < .75 & abs(dMEAN) > .25 & dSD > 3 & pmax(co1,co2) > .0225 ) {   
      if(   abs(dINT) < .75 & abs(dMEAN) > .25 & dSD > 3 & (pmax(co1,co2) > .01 | pMed < .05)  ) {   
        yUse<- x >= .3
        y2<-y[yUse]
        x2<-x[yUse]
        ##ySplit<-kmeans(y[yUse],2)$cluster
        kSolutions<-lapply(1:50,function(x) {
          s1<-kmeans(y2,2)
          s2<-which.min(s1$centers)
          names(s1$cluster)[s1$cluster==s2]
        })
        s1<-unique(unlist(kSolutions))
        s1<-s1[ colSums(do.call("rbind",lapply(kSolutions,function(x) { s1 %in% x }))) >= 25]
        ySplit<-vector("integer",length=length(y2))+2
        names(ySplit)<-names(y2)
        ySplit[s1]<-1
        yUse1<-ySplit == 1
        yUse2<-ySplit == 2
        useALL<-FALSE
      } else {
        useALL<-TRUE
        yUse1<-rep(TRUE,length(x))
        yUse2<-rep(FALSE,length(x))
      }
    }
    #gather intercept+slope
    iOneT<-NA
    iOneN<-NA
    iTwoT<-NA
    iTwoN<-NA
    ww<-rep(1,length(x))
    if(useALL) {
      lmFull1<-lm(y~x)
      lmFull2<-lm(y~c(1-x))
      meNorm<-lmFull1$coefficients[1]+lmFull1$residuals
      meNorm[meNorm < 0]<-0
      meNorm[meNorm > 1]<-1
      meTum<-lmFull2$coefficients[1]+lmFull2$residuals
      meTum[meTum < 0]<-0
      meTum[meTum > 1]<-1
      iOneN<-as.numeric(lmFull1$coefficients)
      iOneT<-as.numeric(lmFull2$coefficients)
    }
    if(!useALL) {
      ##free intercept
      lmFull1.1<-lm(y2[yUse1]~x2[yUse1])
      lmFull1.2<-lm(y2[yUse1]~c(1-x2[yUse1]))

      lmFull2.1<-lm(y2[yUse2]~x2[yUse2])
      lmFull2.2<-lm(y2[yUse2]~c(1-x2[yUse2]))

      ##force common intercept for normal line = average of two estimates
        ##use this line estimate to do proximity assignment of remaning low Tumor% samples
        ##if common intercept is reasonably close
      if( abs( lmFull1.1$coefficients[1]-lmFull2.1$coefficients[1] ) < .5 ) {
        int1<-lmFull1.1$coefficients[1]-mean( c(lmFull1.1$coefficients[1],lmFull2.1$coefficients[1]) )
        int2<-lmFull2.1$coefficients[1]-mean( c(lmFull1.1$coefficients[1],lmFull2.1$coefficients[1]) )
        int1<-rep(int1,sum(yUse1))
        int2<-rep(int2,sum(yUse2))
        lmFull1.1<-lm(y2[yUse1]~x2[yUse1]+offset(int1))
        lmFull2.1<-lm(y2[yUse2]~x2[yUse2]+offset(int2))
      }
#      ##function for shortest distance to line
#       dist2d <- function(a,b,c) {
#       ya1 <- a[2]
#       yb1 <- a[1]
#       ya2 <- b[2]
#       yb2 <- b[1]
#       ##do not correct above .975 or 0.025
#       d<-which.min( c( (abs(yb1+ya1*c[1]-c[2])/sqrt(ya1^2+1)),
#        (abs(yb2+ya2*c[1]-c[2])/sqrt(ya2^2+1)) ) )
#       ##do not correct below intercept of other line 
#       if( d==1 & b[1] <= c[2] & c[1] <= 1 & c[1] < .3 ) d<-2
#       if( d==2 & a[1] >= c[2] & c[1] <= 1 & c[1] < .3 ) d<-1       
#       return(d)
#      }
      ##function for shortest distance to line
       dist2d <- function(a,b,c) {
       ya1 <- a[2]
       yb1 <- a[1]
       ya2 <- b[2]
       yb2 <- b[1]
       ##do not correct above .975 or 0.025
       d<-which.min( c( (abs(yb1+ya1*c[1]-c[2])/sqrt(ya1^2+1)),
        (abs(yb2+ya2*c[1]-c[2])/sqrt(ya2^2+1)) ) )
       ##do not correct below intercept of other line
       frx<-function(x,cc) x*cc[2]+cc[1]  
       if( d==1  & frx(c[1],b) < c[2] ) d<-2
       if( d==2  & frx(c[1],a) > c[2] ) d<-1
       return(d)
      }
      ww<-apply(cbind(x,y),1,function(zz) {
        dist2d(a=lmFull1.1$coefficients,b=lmFull2.1$coefficients,c=zz)
      })

      yUse1<- ww == 1
      yUse2<- ww == 2
      ##recalculate lm's
      lmFull1.1<-lm(y[yUse1]~x[yUse1])
      lmFull1.2<-lm(y[yUse1]~c(1-x[yUse1]))

      lmFull2.1<-lm(y[yUse2]~x[yUse2])
      lmFull2.2<-lm(y[yUse2]~c(1-x[yUse2]))


      meNorm<-rep(NA,length(x))
      meNorm[yUse1]<-lmFull1.1$coefficients[1]+lmFull1.1$residuals
      meNorm[yUse2]<-lmFull2.1$coefficients[1]+lmFull2.1$residuals
      meNorm[meNorm < 0]<-0
      meNorm[meNorm > 1]<-1

      meTum<-rep(NA,length(x))
      meTum[yUse1]<-lmFull1.2$coefficients[1]+lmFull1.2$residuals
      meTum[yUse2]<-lmFull2.2$coefficients[1]+lmFull2.2$residuals
      meTum[meTum < 0]<-0
      meTum[meTum > 1]<-1

      ##gather intercepts
      iOneN<-as.numeric(lmFull1.1$coefficients)
      iOneT<-as.numeric(lmFull1.2$coefficients)
      iTwoN<-as.numeric(lmFull2.1$coefficients)
      iTwoT<-as.numeric(lmFull2.2$coefficients)

    }
    return( list(oneGroup=useALL,
      globalCorr=corALL,
      line1N=iOneN,
      line2N=iTwoN,
      line1T=iOneT,
      line2T=iTwoT,
      mu=c(mu1,mu2),
      sd=c(sd1,sd2),
      dDiff=aOL,
      cDiff=dINT,
      mDiff=dMEAN,
      sDiff=dSD,
      corr1=co1,
      corr2=co2,
      line1=lm1c,
      line2=lm2c,
      methTum=y,
      methCalTum=round(meTum,3),
      methCalAvgDelta=mean( abs(y-meTum) ),
      methCalNorm=round(meNorm,3),
      groups=ww
      )
    )
  } else {
    return( list(oneGroup="NA",
      globalCorr=NA,
      line1N=rep(NA,2),
      line2N=rep(NA,2),
      line1T=rep(NA,2),
      line2T=rep(NA,2),
      mu=rep(NA,2),
      sd=rep(NA,2),
      dDiff=NA,
      cDiff=NA,
      mDiff=NA,
      sDiff=NA,
      corr1=NA,
      corr2=NA,
      line1=rep(NA,2),
      line2=rep(NA,2),
      methTum=y,
      methCalTum=y,
      methCalAvgDelta=NA,
      methCalNorm=NA,
      groups=NA
      )
    )
  }
}

#plot(doLmTests2(y)$methTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=doLmTests2(y)$groups)
#plot(gg[[217]]$methCalTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=gg[[217]]$groups)
#plot(doLmTests2(y)$methCalTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=doLmTests2(y)$groups)
#plot(doLmTests2(y)$methTum,x,col=doLmTests2(y)$groups)
#
#par(mfrow=c(3,2))
#plot(x,doLmTests2(y)$methTum,col=doLmTests2(y)$groups,pch=16)
#plot(x,gg[[217]]$methTum,col=gg[[217]]$groups,pch=16)
#plot(x,doLmTests2(y)$methCalTum,col=doLmTests2(y)$groups,pch=16)
#plot(x,gg[[217]]$methCalTum,col=gg[[217]]$groups,pch=16)
#plot(gg[[217]]$methCalTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=gg[[217]]$groups,pch=16)
#plot(doLmTests2(y)$methCalTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=doLmTests2(y)$groups,pch=16)
#
#table(doLmTests2(y)$groups,gg[[217]]$groups)
#
#table(doLmTests2(y)$groups, clinAnno[samples_use,"BRCA1_PromMetPc_Class"])
#
#      ww<-apply(cbind(x,y),1,function(zz) {
#        a<-dist2d(c(0,lmFull1.1$coefficients[1]),c(1,sum(lmFull1.1$coefficients)),zz)
#        b<-dist2d(c(0,lmFull2.1$coefficients[1]),c(1,sum(lmFull2.1$coefficients)),zz)
#        (c(a,b))
#      })
#


