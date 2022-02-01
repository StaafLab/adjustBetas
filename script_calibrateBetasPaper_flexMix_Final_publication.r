#####======================================================================#####
### Short paper on beta calibration based on tumor content
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncology & Pathology

################################################################################
##create work directories in default location

##work
HOME<-"~/hdd1"
DATA<-"~/hdd1/tcgaBrca"
GIT<-"~/Documents/adjustBetas"

##home
HOME<-"I:/data"
DATA<-"I:/data/tcgaBrca"
GIT<-"F:/gitProjects/adjustBetas"

##create data directories
if( !file.exists( paste0(HOME,"/","adjustBetas_paper") )   ) {
  dir.create(paste0(HOME,"/","adjustBetas_paper/me"),recursive=TRUE)
}

HOME<-"~/hdd1/adjustBetas_paper"
HOME<-"I:/data/adjustBetas_paper"

################################################################################
##load required packages

if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager") }
 
if(!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges") }

library(GenomicRanges)

if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer") }

library(RColorBrewer)

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap") }

library(pheatmap)

if(!requireNamespace("magick", quietly = TRUE)) {
  install.packages("magick") }

library(magick)

if(!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr") }

library(httr)

if(!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl") }

library(readxl)

if(!requireNamespace("ggalluvial", quietly = TRUE)) {
  install.packages("ggalluvial") }

library(ggalluvial)

if(!requireNamespace("flexmix", quietly = TRUE)) {
  install.packages("flexmix") }

library("flexmix")

if(!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret") }

library("caret")

if(!requireNamespace("doParallel", quietly = TRUE)) {
	install.packages("doParallel") }

library(doParallel)

if(!requireNamespace("parallel", quietly = TRUE)) {
	install.packages("parallel") }

library(parallel)

##source - flexmix loaded on source
source(paste0(GIT,"/function_correctBetas.r"))

################################################################################
##Get TCGA data object

load(file=paste0(DATA,"/","finalWorkspace_atacCnGexMeWes_withAnnotations.RData"))

ls()
#  [1] "adjustBeta" "annoObj"    "betaAdj"    "betaNorm"   "betaOrig"   "DATA"       "dataAtac"   "dataCn"     "dataMut"    "dataSeg"   
# [11] "gexCounts"  "gexFpkm"    "gexUq"      "GIT"        "HOME"       "MANIFEST"   "sampleAnno" "TUMOR_TYPE"

str(annoObj)
# 'data.frame':   421368 obs. of  65 variables:
#  $ illuminaID                   : chr  "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
#  $ chr                          : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start                        : num  69591 864703 870161 877159 898803 ...
#  $ end                          : num  69592 864704 870162 877160 898804 ...
#  $ hasUCSCknownGeneOverlap      : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ nameUCSCknownGeneOverlap     : chr  "OR4F5" "" "FAM41C" "FAM41C" ...
#  $ numberUCSCknownGeneOverlap   : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownGenePromOverlap  : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneUp5kbOverlap : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneDn5kbOverlap : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxPromOverlap    : num  0 0 1 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxUp5kbOverlap   : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ hasUCSCknownTxDn5kbOverlap   : num  1 1 0 1 0 0 0 0 0 0 ...
#  $ ucscKnownGeneIsDistal        : num  0 1 1 0 1 1 1 1 1 1 ...
#  $ ucscKnownGeneIsGenic         : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ ucscKnownGeneIsDistalNonGenic: num  0 1 0 0 1 1 1 1 1 1 ...
#  $ ucscKnownGeneIsGeneBody      : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ hasGeneOverlap               : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ nameGeneOverlap              : chr  "ENSG00000186092|OR4F5" "" "ENSG00000230368|FAM41C" "" ...
#  $ numberGeneOverlap            : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ hasUp5kbOverlap              : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ nameUp5kbOverlap             : chr  "" "" "ENSG00000234711|TUBB8P11" "" ...
#  $ numberUp5kbOverlap           : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ hasDn5kbOverlap              : num  1 0 0 1 0 0 0 0 0 0 ...
#  $ nameDn5kbOverlap             : chr  "ENSG00000186092|OR4F5" "" "" "ENSG00000234711|TUBB8P11" ...
#  $ numberDn5kbOverlap           : num  1 0 0 1 0 0 0 0 0 0 ...
#  $ hasPromOverlap               : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ namePromOverlap              : chr  "" "" "" "ENSG00000230368|FAM41C" ...
#  $ numberPromOverlap            : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ isPromMostVariable           : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ namePromMostVariable         : chr  "" "" "" "ENSG00000230368|FAM41C" ...
#  $ numberPromMostVariable       : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ hasAtacOverlap               : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ nameAtacOverlap              : chr  "" "" "chr1:869670-870169" "" ...
#  $ numberAtacOverlap            : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ isAtacMostVariable           : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ isDistal                     : num  0 1 0 0 1 1 1 0 0 0 ...
#  $ isGenic                      : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ isDistalNonGenic             : num  0 1 0 0 1 1 1 0 0 0 ...
#  $ isGeneBody                   : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ weberOE                      : num  0.248 0.249 0.637 0.163 0.565 ...
#  $ weberClass                   : chr  "LCP" "LCP" "HCP" "LCP" ...
#  $ saxonovOE                    : num  0.238 0.268 0.473 0.202 0.441 ...
#  $ saxonovClass                 : chr  "LCG" "LCG" "HCG" "LCG" ...
#  $ cgCount300Bp                 : int  6 2 27 4 9 13 16 8 18 16 ...
#  $ cgPer100Bp                   : num  2 0.667 9 1.333 3 ...
#  $ isCgi                        : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ isShore                      : num  0 1 0 0 0 0 0 0 1 1 ...
#  $ isOcean                      : num  1 0 0 1 1 1 1 1 0 0 ...
#  $ cgiClass                     : chr  "ocean" "shore" "cgi" "ocean" ...
#  $ featureClass                 : chr  "proximal dn" "distal" "proximal up" "promoter" ...
#  $ featureClassUcsc             : chr  "proximal dn" "distal" "distal body" "promoter" ...
#  $ illuminaCpG_CpH_Probe        : chr  "cg" "cg" "cg" "cg" ...
#  $ encodeCre                    : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ encodeCreCategory            : chr  "" "" "PLS" "" ...
#  $ encodeCreCtcf                : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ encodeChipMeanNcellPeaks     : num  0 0 4.73 0 1 ...
#  $ encodeChipUniqueTfPeaks      : num  0 0 30 0 9 7 7 1 1 0 ...
#  $ hasAnyRepeatOverlap          : num  0 0 0 0 0 0 0 0 1 1 ...
#  $ hasMultiRepeatOverlap        : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasOneRepeatOverlap          : num  0 0 0 0 0 0 0 0 1 1 ...
#  $ hasRepeatOverlap             : num  0 0 0 0 0 0 0 0 1 1 ...
#  $ repeatClass                  : chr  "" "" "" "" ...
#  $ repeatFamily                 : chr  "" "" "" "" ...
#  $ repeatName                   : chr  "" "" "" "" ...

##get tumor fraction
#fracTum<-read.table(file=paste0(DATA,"/me/samplePurityVector.txt"),header=TRUE,as.is=T,sep="\t")
fracTum<-read.table(file=paste0(DATA,"/me/sampleTumorDnaFractionVector.txt"),header=TRUE,as.is=T,sep="\t")

str(fracTum)
# 'data.frame':   630 obs. of  2 variables:
#  $ sample            : chr  "TCGA-3C-AAAU-01" "TCGA-3C-AALI-01" "TCGA-3C-AALJ-01" "TCGA-3C-AALK-01" ...
#  $ tumor_dna_fraction: num  0.85 0.73 0.69 0.62 0.64 0.53 0.89 0.47 0.56 0.65 ...

rownames(fracTum)<-fracTum$sample

identical(rownames(fracTum),sub(".$","",colnames(betaOrig)))
#[1] TRUE

rownames(fracTum)<-colnames(betaOrig)

HOME<-"~/hdd1/adjustBetas_paper"
HOME<-"I:/data/adjustBetas_paper"

################################################################################
##load data

##Load GSE67919 - 96 Normal breast samples
load("I:/data/luTnbc/data/files_js/GSE67919_Beta.RData")
load("I:/data/luTnbc/data/files_js/GSE67919_Annotations.RData")

annotations_norm<-annotations
beta_norm<-beta
rm(annotations,beta)

##do Ringer-correction

##get updated annotations
hm.manifest<-readRDS( paste0(DATA,"/me/","HM450.hg38.manifest.rds") )

##remove "unreliable" probes
probesKeep<-cbind(id=names(hm.manifest),
  chr=as.character(seqnames(hm.manifest)),
  start=as.integer(start(hm.manifest)),
  end=as.integer(end(hm.manifest)),
  isBad=hm.manifest$MASK_general,
  probeType=hm.manifest$probeType,
  designType=hm.manifest$designType)

head(probesKeep)
#      id           chr    start   end     isBad  probeType designType
# [1,] "cg13869341" "chr1" "15865" "15866" "TRUE" "cg"      "I"       
# [2,] "cg14008030" "chr1" "18827" "18828" "TRUE" "cg"      "II"      
# [3,] "cg12045430" "chr1" "29407" "29408" "TRUE" "cg"      "I"       
# [4,] "cg20826792" "chr1" "29425" "29426" "TRUE" "cg"      "I"       
# [5,] "cg00381604" "chr1" "29435" "29436" "TRUE" "cg"      "I"       
# [6,] "cg20253340" "chr1" "68849" "68850" "TRUE" "cg"      "II"      

probesKeep<-as.data.frame(probesKeep,stringsAsFactors=FALSE)
rownames(probesKeep)<-probesKeep$id

str(probesKeep)
# 'data.frame': 485577 obs. of  7 variables:
#  $ id        : chr  "cg13869341" "cg14008030" "cg12045430" "cg20826792" ...
#  $ chr       : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start     : chr  "15865" "18827" "29407" "29425" ...
#  $ end       : chr  "15866" "18828" "29408" "29426" ...
#  $ isBad     : chr  "TRUE" "TRUE" "TRUE" "TRUE" ...
#  $ probeType : chr  "cg" "cg" "cg" "cg" ...
#  $ designType: chr  "I" "II" "I" "I" ...

table(probesKeep$isBad,probesKeep$designType)
  #            I     II
  # FALSE 115147 306286
  # TRUE   20354  43790

probesKeep$start<-as.integer(probesKeep$start)
probesKeep$end<-as.integer(probesKeep$end)
probesKeep$isBad<-as.logical(probesKeep$isBad)

isect<-intersect( rownames(betaOrig),rownames(probesKeep) )

beta_norm<-beta_norm[isect,]

probesKeep<-probesKeep[isect,]

all.equal( rownames(beta_norm) , rownames(probesKeep) )
#[1] TRUE
all.equal( rownames(betaOrig) , rownames(probesKeep) )
#[1] TRUE

rm(isect)

##beta object
dim(beta_norm)
#[1] 421368     96

table(probesKeep$designType)
#      I     II 
# 115122 306246 

##define probe subsets
probes_I<- probesKeep$designType == "I"
probes_II<- probesKeep$designType == "II"

table(probes_I,probes_II)
#         probes_II
# probes_I  FALSE   TRUE
#    FALSE      0 306246
#    TRUE  115122      0

##identify zero and one peak for respective probe sets and calibrate range of one to other

##updated function for handling samples with non-canonical beta distributions
findMaxima<-function(x,granularity=512,bw=.02,adjust=1,from=0,to=1,kernel="epanechnikov",verbose=TRUE) {
  aa<-density(x,kernel=kernel,n=granularity,bw=bw,adjust=adjust,from=from,to=to,na.rm=T)
  res<-vector(length=length(aa$x))
  for(i in 2:length(res[-1])) {
    res[i]<-aa$y[i]>aa$y[i-1] & aa$y[i]>aa$y[i+1]
  }
  y<-aa$x[res][order(aa$y[res],decreasing=T)] 
  if(verbose) {
    cat("found peaks at:",y,"\n")
    cat("choosing min/max:",range(y),"\n") ###CHANGED here to pick out min-max peaks as these are typically the correct ones
  }
  range(y) 
}

##get original matrix
betaData<-beta_norm

##define results matrix
meth.cal<-matrix(0,ncol=ncol(betaData),nrow=nrow(betaData))
rownames(meth.cal)<-rownames(betaData)
colnames(meth.cal)<-colnames(betaData)

##define vector for flagging IDs
flaggedIDs<-NULL

pdf(file=paste0(HOME,"/me/","infiniumAdjustmentPlots.pdf"),width=8,height=8)
par(mfrow=c(2,2),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
sink(paste0(HOME,"/me/","infiniumAdjustmentStats.txt"))

for ( i in 1:ncol(betaData)) {
  cat("#####Sample ",i,":",colnames(betaData)[i],"\n")
  ##do pre-calibaration I vs II plot
  plot(density(betaData[probes_I,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),main=colnames(betaData)[i],xlim=c(-.05,1.05),col="orange",lwd=3,xlab="pre-adjustment beta",cex.main=.9,ylim=c(0,8))
  lines(density(betaData[probes_II,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),col="darkgreen",lwd=3)

  ##define maxima
  tempMax_I<-findMaxima(betaData[probes_I,i])
  tempMax_II<-findMaxima(betaData[probes_II,i])
  ##function not oriented so peak near zero can be in vector slot 2 -> fix
  tempMax_I<-tempMax_I[order(tempMax_I)]
  tempMax_II<-tempMax_II[order(tempMax_II)]
  ##add plot "FLAG" if peaks not well separated and close to right place..
  FLAG<-FALSE
  if( tempMax_I[1] >.2 | tempMax_I[2] <.8 | tempMax_II[1] >.2 | tempMax_II[2] <.8 ) {
    FLAG<-TRUE
    flaggedIDs<-c(flaggedIDs,colnames(betaData)[i])
  }

  ##add to plot
  abline(v=tempMax_I,col="orange",lwd=2,lty=2)
  abline(v=tempMax_II,col="darkgreen",lwd=2,lty=2)
  text(x=.5,y=c(3.5,3),labels=c("Infinium_I ","Infinium_II"),col=c("orange","darkgreen"),font=2,cex=1)
  if(FLAG)   text(x=.5,y=.5,labels=c("FLAGGED"),col=2,font=2,cex=2)

  if(FLAG) {
    cat("Something funky with this sample!!!\n")
    cat("infinium_I peaks:",tempMax_I,"\n",sep="\t")
    cat("infinium_II peaks:",tempMax_II,"\n",sep="\t")
  } else {
    cat("infinium_I peaks:",tempMax_I,"\n",sep="\t")
    cat("infinium_II peaks:",tempMax_II,"\n",sep="\t")
  }

  ##do calibration scale both to same range..
  meth.cal[probes_I,i]<-(betaData[probes_I,i]-tempMax_I[1])/(tempMax_I[2]-tempMax_I[1])
  meth.cal[probes_II,i]<-(betaData[probes_II,i]-tempMax_II[1])/(tempMax_II[2]-tempMax_II[1])

  ##cap ends so that nothing over 1 or below zero..
  meth.cal[probes_I,i][meth.cal[probes_I,i] < 0 ]<-0
  meth.cal[probes_I,i][meth.cal[probes_I,i] > 1 ]<-1
  meth.cal[probes_II,i][meth.cal[probes_II,i] < 0 ]<-0
  meth.cal[probes_II,i][meth.cal[probes_II,i] > 1 ]<-1

  ##do post-calibration I vs II plot
  plot(density(meth.cal[probes_I,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),main=colnames(meth.cal)[i],xlim=c(-.05,1.05),col="orange",lwd=3,xlab="post-adjustment beta",cex.main=.9,ylim=c(0,8))
  lines(density(meth.cal[probes_II,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),col="darkgreen",lwd=3)
  text(x=.5,y=c(3.5,3),labels=c("Infinium_I ","Infinium_II"),col=c("orange","darkgreen"),font=2,cex=1)
  abline(v=findMaxima(meth.cal[probes_I,i],verbose=F),col="orange",lwd=2,lty=2)
  abline(v=findMaxima(meth.cal[probes_II,i],verbose=F),col="darkgreen",lwd=2,lty=2)

  cat("\n")
} ; rm( list = c("i","tempMax_I","tempMax_II","FLAG") )
cat("##FlaggedIDs:\n")
cat(flaggedIDs,sep="\n")
sink()
dev.off()

##some peaks not where they should be??
flaggedIDs
NULL

writeLines(flaggedIDs,con=paste0(HOME,"/me/","infiniumAdjustmentFlaggedIDs.txt"))

if(!is.null(flaggedIDs)) {
  write(paste(sampleMat[flaggedIDs,"Sample_ID"],flaggedIDs,sep="\t"),file=paste0(HOME,"/me/","methylationSampleBlacklist.txt"),append=TRUE)
}

rm(probes_I)
rm(probes_II)
rm(findMaxima)

##write adjusted table to file
con<-gzfile(paste0(HOME,"/me/","infiniumAdjustedNormalizedData.txt.gz"),"w")
write.table(meth.cal,file=con,sep="\t",quote=F)
close(con)
rm(con)

beta_norm_adj<-meth.cal
rm(meth.cal)

################################################################################
###Beta correction has been done already on all CpGs -> see how results change now..

##filter X and get top5k
betaData<-betaOrig

tmp<-fracTum$tumor_dna_fraction
names(tmp)<-rownames(fracTum)

all.equal(names(tmp),colnames(betaData))
#[1] TRUE

fracTum<-tmp
rm(tmp)

all.equal(rownames(annoObj),rownames(betaData))
#[1] TRUE

table(annoObj$chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
# 41098 21206 25633 21624 10540 13204 13119 18627 24356  5280 21984 30550  9354 
# chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY 
#  3378  7222 22315 17646 21279 31170 25024 18194  8605  9718   242 

table(annoObj$illuminaCpG_CpH_Probe)
#     cg     ch 
# 418858   2510 

betaData<-betaData[ !annoObj$chr %in% c("chrX","chrY") &
   annoObj$illuminaCpG_CpH_Probe == "cg"
,]

varF<-apply(betaData,1,sd)

varF<-varF >= quantile(varF,1-(5000/nrow(betaData)))

table(varF)
# varF
#  FALSE   TRUE 
# 403966   5000 

varF<-rownames(betaData)[varF]

length(intersect(varF,rownames(beta_norm)))
#[1] 5000

##create data set and run function
testDat2<-betaData[varF,]

str(testDat2)
 # num [1:5000, 1:630] 0.341 0.057 0.012 0.022 0 0.004 0.021 0.892 0.022 0.163 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg09248054" "cg25340711" "cg06443533" "cg16601494" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

################################################################################
###correct betas using multicore

no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

##estimated runtime is ~500 sec on 7 cores

clusterEvalQ(cl, {
  library("flexmix")
})

#clusterExport(cl,list("adjustBeta","betaData","fracTum"))

betaRun<-cbind(seed=1:nrow(testDat2),testDat2)
betaNames<-colnames(testDat2)

#clusterSetRNGStream(cl, 20200918) ##will not make exactly replicable..
res<-parRapply(cl = cl, betaRun, adjustBeta,purity=fracTum,snames=betaNames,seed=TRUE)

res2<-parRapply(cl = cl, betaRun, adjustBeta,purity=fracTum,snames=betaNames,seed=TRUE)

table(unlist(lapply(res,function(x) x$n.groups)))
  #  2    3 
  # 20 4980 

table(unlist(lapply(res,function(x) x$n.groups)),unlist(lapply(res2,function(x) x$n.groups)))
  #      2    3
  # 2   20    0
  # 3    0 4980

rm(betaRun,betaNames)

table( unlist(lapply(1:length(res),function(x) { all( unlist(res[[x]],use.names=FALSE) == unlist(res[[x]],use.names=FALSE) ) }) ) )
# TRUE 
# 5000 
table( unlist(lapply(1:length(res),function(x) { all( unlist(res[[x]],use.names=FALSE) == unlist(res2[[x]],use.names=FALSE) ) }) ) )
# TRUE 
# 5000 

save(res,file=paste0(HOME,"/results_object_adjustedBetas.RData"))

################################################################################
### Also do infiniumPurify

##InfiniumPurify functions - Adapted from https://github.com/Xiaoqizheng/InfiniumPurify
myasin <- function(x) asin(2*x-1)

InfiniumPurify<-function (tumor.data, normal.data, purity) 
{
    if (missing(tumor.data) | missing(normal.data) | missing(purity)) {
        stop("'tumor.data', 'normal.data' and 'purity' are required.")
    }
    probes = intersect(rownames(tumor.data), rownames(normal.data))
    tumor.sample = intersect(colnames(tumor.data), names(purity))
    normal.sample = colnames(normal.data)
    purity = purity[tumor.sample]
    if (length(normal.sample) < 20 | length(tumor.sample) < 20) {
        stop("tumor and normal samples should be more than 20!")
    }
    .get_corrBeta <- function(input) {
        x = as.numeric(input[tumor.sample])
        y = as.numeric(input[normal.sample])
        type = c(rep("Tumor", length(x)), rep("Normal", length(y)))
        data = data.frame(beta = c(x, y), type = type)
        Y = myasin(data$beta)
        X = c(1 - purity, rep(0, length(y)))
        fit = lm(Y ~ X)
        tmp = resid(fit) + coef(fit)[1]
        beta.pred = (sin(tmp) + 1)/2
        beta.pred
    }
    all.data = cbind(tumor.data[probes, tumor.sample], normal.data[probes, 
        ])
    probes.rmna = probes[rowSums(is.na(all.data)) == 0]
    all.data.corr = t(apply(all.data[probes.rmna, ], 1, .get_corrBeta))
    tumor.data.corr = all.data.corr[, 1:length(tumor.sample)]
    colnames(tumor.data.corr) = tumor.sample
    tumor.data.corr
}

temp4<-beta_norm[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))

beta.infPur<-InfiniumPurify(tumor.data=testDat2,normal.data=temp4,purity=fracTum)

str(beta.infPur)
 # num [1:5000, 1:630] 0.35269 0.03614 0.01306 0.01478 0.00128 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg09248054" "cg25340711" "cg06443533" "cg16601494" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

################################################################################
###gather adjusted data and calculate new data from it..

##check stats
temp4<-do.call("rbind",lapply(res,function(x) x$y.tum))
#rownames(temp4)<-rownames(testDat2)
temp5<-do.call("rbind",lapply(res,function(x) x$y.orig))
#rownames(temp5)<-rownames(testDat2)

table(apply(temp4,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp5,1,function(x) sum(is.na(x))))
#   0
#5000

quantile(testDat2)
#    0%   25%   50%   75%  100% 
# 0.000 0.155 0.470 0.743 1.000 

quantile(temp4)
#    0%   25%   50%   75%  100% 
# 0.000 0.053 0.526 0.924 1.000 

##Plot histogram of betas before and after correction
pdf(paste0(HOME,"/top5k_betaDistribution_tumors_beforeAfterCorrection.pdf"),width=8,height=8,useDingbats=F)

par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot(1,xlim=range(round(density(temp4)$x,1)),ylim=range(round(density(temp4)$y,1))+c(0,.5),type="n",las=1,axes=F,
  xlab="beta",ylab="density"
)
lines(density(temp4),col=2,lwd=2)
lines(density(testDat2),col=1,lwd=2)
lines(density(beta.infPur),col=3,lwd=2)
axis(1,lwd=2,las=1,at=seq(0,1,by=.2))
axis(2,lwd=2,las=1)
legend("topright",legend=c("unadjusted beta","adjusted beta","infiniumPurify"),col=c(1,2,3),lwd=2,bty="n")
dev.off()

##Plot histogram of betas before and after correction
pdf(paste0(HOME,"/top5k_correlationTumorFrac_beforeCorrection.pdf"),width=8,height=8,useDingbats=F)
par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot(unlist(lapply(res,function(x) abs(x$glob.cor) )),
  unlist(lapply(res,function(x) abs(x$avg.betaDiff) )),
  pch=16,cex=.5,xlim=c(0,1),ylim=c(0,1),
  main="Average beta correction vs global correlation",
  xlab="absolute global correlation unadjusted beta - tumor fraction",
  ylab="mean absolute beta difference pre-post correction",
  axes=F#,type="n"
)
abline(a=0,b=1,lwd=3,col=2,lty=2)
axis(1,lwd=2,las=1)
axis(2,lwd=2,las=1)
dev.off()

rm(temp4,temp5)

################################################################################
###plot top5k clusters - unadjuster order

##check stats
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))
temp4<-beta_norm[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))

table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp3,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp4,1,function(x) sum(is.na(x))))
#   0
#5000

table(annoObj[rownames(temp1),"chr"])
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4 
 #  557   272   213   255   156   144    98   169   272    69   205   365   130    43    54   248   192 
 # chr5  chr6  chr7  chr8  chr9 
 #  370   392   386   353    57 

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp3) ),method="ward.D"),5)
c2<-unique(c1[hclust( as.dist( 1-cor(temp3) ),method="ward.D")$order])
r1<-hclust( dist(temp3),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp3) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)

sample_anno<-data.frame(unadj5000=as.character(c1),
  adj5000=as.character(c4),
  ER=sampleAnno$ER,
  PR=sampleAnno$PR,
  HER2=sampleAnno$HER2,
  TNBC=as.character(as.integer(sampleAnno$TNBC)),
  PAM50=sampleAnno$pam50.full,stringsAsFactors=FALSE
  )
rownames(sample_anno)<-colnames(temp1)
sample_anno<-sample_anno[,ncol(sample_anno):1]

# my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00","NA"="white"),
#     adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00","NA"="white"),
#     ER = c("NA"="#FF7F00",
#       "Negative"="#E41A1C",
#       "Positive"="#4DAF4A"
#       ),
#     PR = c("NA"="#FF7F00",
#       "Negative"="#E41A1C",
#       "Positive"="#4DAF4A"
#       ),
#     HER2 = c("NA"="#FF7F00",
#       "Negative"="#E41A1C",
#       "Positive"="#4DAF4A"
#       ),
#     TNBC = c("1"="black","0"="lightgrey","NA"="white"),
#     PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#377EB8","NA"="white"),stringsAsFactors=FALSE
#   )

my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    ER = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    PR = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    HER2 = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A",
      "NA"="white"
      ),
    TNBC = c("1"="black","0"="lightgrey"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white"),stringsAsFactors=FALSE
  )

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_unadjClust_unadjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp3,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_unadjClust_adjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_unadjClust_InferredNormalBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_unadjClust_infiniumPurifyBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(beta.infPur,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

sample_anno<-data.frame("unadj5000"=rep("NA",ncol(temp4)),
  "adj5000"=rep("NA",ncol(temp4)),
  "ER"=rep("NA",ncol(temp4)),
  "PR"=rep("NA",ncol(temp4)),
  "HER2"=rep("NA",ncol(temp4)),
  "TNBC"=rep("NA",ncol(temp4)),
  "PAM50"=rep("NA",ncol(temp4)),stringsAsFactors=FALSE
)
rownames(sample_anno)<-colnames(temp4)
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour<-list("unadj5000"=c("NA"="white"),
  "adj5000"=c("NA"="white"),
  "ER"=c("NA"="white"),
  "PR"=c("NA"="white"),
  "HER2"=c("NA"="white"),
  "TNBC"=c("NA"="white"),
  "PAM50"=c("NA"="white"))

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_unadjClust_normalBeta.tiff"),width=9*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp4,cluster_rows = r1, 
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, GSE67919 external normal, default clust , pearC/euclR"
  ,annotation_col=sample_anno,annotation_colors=my_colour
  
)
dev.off()

################################################################################
###plot top5k clusters - adjuster order

##check stats
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))
temp4<-beta_norm[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))

table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp3,1,function(x) sum(is.na(x))))
#   0
#5000

table(annoObj[rownames(temp1),"chr"])
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4 
 #  557   272   213   255   156   144    98   169   272    69   205   365   130    43    54   248   192 
 # chr5  chr6  chr7  chr8  chr9 
 #  370   392   386   353    57 

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-unique(c1[hclust( as.dist( 1-cor(temp1) ),method="ward.D")$order])
r1<-hclust( dist(temp1),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp3) ),method="ward.D"),5)

sample_anno<-data.frame(adj5000=as.character(c1),
  unadj5000=as.character(c4),
  ER=sampleAnno$ER,
  PR=sampleAnno$PR,
  HER2=sampleAnno$HER2,
  TNBC=as.character(as.integer(sampleAnno$TNBC)),
  PAM50=sampleAnno$pam50.full,stringsAsFactors=FALSE
  )
rownames(sample_anno)<-colnames(temp1)
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    ER = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    PR = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    HER2 = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A",
      "NA"="white"
      ),
    TNBC = c("1"="black","0"="lightgrey"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white"),stringsAsFactors=FALSE
  )

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_adjClust_unadjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp3,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_adjClust_InferredNormalBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_adjClust_infiniumPurifyBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(beta.infPur,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"infiniumPurify\", adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

sample_anno<-data.frame("unadj5000"=rep("NA",ncol(temp4)),
  "adj5000"=rep("NA",ncol(temp4)),
  "ER"=rep("NA",ncol(temp4)),
  "PR"=rep("NA",ncol(temp4)),
  "HER2"=rep("NA",ncol(temp4)),
  "TNBC"=rep("NA",ncol(temp4)),
  "PAM50"=rep("NA",ncol(temp4)),stringsAsFactors=FALSE
)
rownames(sample_anno)<-colnames(temp4)
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour<-list("unadj5000"=c("NA"="white"),
  "adj5000"=c("NA"="white"),
  "ER"=c("NA"="white"),
  "PR"=c("NA"="white"),
  "HER2"=c("NA"="white"),
  "TNBC"=c("NA"="white"),
  "PAM50"=c("NA"="white"))

tiff(paste0(HOME,"/top5k_heatmap_pear_eucl_adjClust_normalBeta.tiff"),width=9*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp4,cluster_rows = r1, 
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, GSE67919 external normal, default clust , pearC/euclR"
  ,annotation_col=sample_anno,annotation_colors=my_colour
  
)
dev.off()

################################################################################
################################################################################

##save.image(paste0(HOME,"/tempWorkspace_210608.RData")) ##remove later
#load(paste0(HOME,"/tempWorkspace_210608.RData"))

################################################################################
################################################################################
##Create Fig 1 panels and image

##Panel 1 - BRCA1 vs Infiltration

##Panel 2 - Schematic flexmix identify populations - infer pure tumor + infer normal


################################################################################
################################################################################
##Create Fig 2 panels and image

##Panel 1 - Adj clust 500 p x 100 rnd iter - 2 group solution - 1 example

##Panel 2 - Discrimination Basal-Luminal | Adj vs dichotomized vs none


##Do basal vs luminal split in 100 * 500 random CpG sets
  ##extract -logP(fisher) hclust 2-split vs PAM50 basal-luminal
  ##plot distribution of p-values for corrected vs uncorrected comparison

all(rownames(sampleAnno)==colnames(betaData))
#[1] TRUE
all(rownames(sampleAnno)==names(fracTum))
#[1] TRUE

##number of lists
N<-100
resMat<-matrix(ncol=12,nrow=N,dimnames=list(1:N,paste(rep(c("acc","sen","spe"),4),rep(c(".raw",".adj",".dic",".ip"),each=3),sep="")) )

##generate 500 random CpG sets
set.seed(20210705)
varF<-apply(betaData,1,sd)
varF<-varF > 0
p_list<-lapply(1:N,function(x) sample(rownames(betaData)[varF],500) )
s_list<-lapply(1:N,function(x) sample(1:1e6,500) )

##do multicore
no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

clusterEvalQ(cl, {
  library("flexmix")
})

##main loop

##metric
getAcc<-function(x,ref) {   
  a1<-sum(x==ref) / length(x)
  a1<-c( a1 , sum( abs(3-x) == ref ) / length(x) )
  a1[which.max(a1)]
}
getSens<-function(x,ref) {   
  a1<-sum(x==2 & ref==2 ) / sum( ref==2 )
  a1<-c( a1 , sum(abs(3-x)==2 & ref==2 ) / sum( ref==2 ) )
  a1[which.max(a1)]
}
getSpec<-function(x,ref) {   
  a1<-sum(x!=2 & ref!=2 ) / sum( ref!=2 )
  a1<-c( a1 , sum(abs(3-x)!=2 & ref!=2 ) / sum( ref!=2 ) )
  a1[which.max(a1)]
}

###RUN
set.seed(20210705)
for (i in 1:N) {
  cat(i," of ",N,"\n")
  ##betaData filtered for chrX/Y
  betaRun<-cbind(seed=s_list[[i]],betaData[p_list[[i]],])
  betaNames<-colnames(betaData)
  ##own method
  b<-parRapply(cl = cl, betaRun, adjustBeta,purity=fracTum,snames=betaNames,seed=TRUE)
  ##infiniumPurify
  temp5<-t( apply(beta_norm[p_list[[i]],],1,function(x) { x[is.na(x)]<-median(x[!is.na(x)]) ; return(x) } ) )
  colnames(temp5)<-colnames(beta_norm)
  c<-InfiniumPurify(tumor.data=betaData[p_list[[i]],],normal.data=temp5,purity=fracTum)

  ##adjusted
  b1<-do.call("rbind",lapply(b,function(x) x$y.tum))
    b1<-b1[!apply(b1,1,function(x) any(is.na(x))),]
    b1<-cutree( hclust( as.dist(1-cor(b1)),"ward.D"), k=2)
    ##unadjusted
    b2<-do.call("rbind",lapply(b,function(x) x$y.orig))
    b2<-b2[!apply(b2,1,function(x) any(is.na(x))),]
    b2<-cutree( hclust( as.dist(1-cor(b2)),"ward.D"), k=2)
    ##dichotomized
    b3<-do.call("rbind",lapply(b,function(x) x$y.orig>.3))
    b3<-b3[!apply(b3,1,function(x) any(is.na(x))),]
    b3<-cutree( hclust( as.dist(1-cor(b3)),"ward.D"), k=2)
    ##infPur
    b4<-c
    b4<-b4[!apply(b4,1,function(x) any(is.na(x))),]
    b4<-cutree( hclust( as.dist(1-cor(b4)),"ward.D"), k=2)
#     resMat[i,1]<- -log10(fisher.test(table(b2,sampleAnno$pam50.full != "Basal"))$p.value)
#     resMat[i,2]<- -log10(fisher.test(table(b1,sampleAnno$pam50.full != "Basal"))$p.value)
#     resMat[i,3]<- -log10(fisher.test(table(b3,sampleAnno$pam50.full != "Basal"))$p.value)
  resMat[i,"acc.raw"]<-getAcc(b2,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"sen.raw"]<-getSens(b2,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"spe.raw"]<-getSpec(b2,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"acc.adj"]<-getAcc(b1,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"sen.adj"]<-getSens(b1,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"spe.adj"]<-getSpec(b1,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"acc.dic"]<-getAcc(b3,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"sen.dic"]<-getSens(b3,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"spe.dic"]<-getSpec(b3,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"acc.ip"]<-getAcc(b4,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"sen.ip"]<-getSens(b4,1+as.integer(sampleAnno$pam50.full == "Basal") )
  resMat[i,"spe.ip"]<-getSpec(b4,1+as.integer(sampleAnno$pam50.full == "Basal") )
}
rm(b,b1,b2,b3,i,N,cl,betaRun,betaNames,getAcc,getSens,getSpec)

save(p_list,file=paste0(HOME,"/top100percentBySd_basalVsLuminalSplitIn100randomSets_usedProbeSets.RData"))
rm(p_list)
save(s_list,file=paste0(HOME,"/top100percentBySd_basalVsLuminalSplitIn100randomSets_usedRngSeeds.RData"))
rm(s_list)
save(resMat,file=paste0(HOME,"/top100percentBySd_basalVsLuminalSplitIn100randomSets_runResults.RData"))
#rm(resMat)

rm(varF)

################################################################################
################################################################################
###Do plots

resMat2<-resMat

pdf(paste0(HOME,"/fig2_top100percentBySd_confusionStats_basalVsLuminalSplitIn100randomSets.pdf"),width=10,height=10,useDingbats=F)
par(fig=c(0,.5,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=F)

boxplot(resMat2[,4]-resMat2[,1],at=1,xlim=c(0.5,11.5),ylim=c(-1.,1.1),width=2,axes=F,lwd=2,
  main="Discrimination of PAM50 Basal vs Luminal split"
)
abline(h=0,lwd=3,col="lightgrey",lty=2)
boxplot(resMat2[,7]-resMat2[,1],at=2,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,10]-resMat2[,1],at=3,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,5]-resMat2[,2],at=5,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,8]-resMat2[,2],at=6,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,11]-resMat2[,2],at=7,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,6]-resMat2[,3],at=9,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,9]-resMat2[,3],at=10,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,12]-resMat2[,3],at=11,add=T,width=2,axes=F,lwd=2)

axis(1,at=c(1:3,5:7,9:11),lwd=2,las=2,cex=1.2,
  labels=c("Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify")
)
axis(2,at=round(seq(-1,1,length.out=5),2),lwd=2,las=1,cex=1.2)
mtext(side=2, "Relative to unadjusted data",font=2,line=2.5,cex=1.2)
lines(x=c(1,3),y=c(.85,.85),lwd=3)
text(2,.85,labels="Accuracy",pos=3)

lines(x=c(5,7),y=c(.85,.85),lwd=3)
text(6,.85,labels="Sensitivity",pos=3)

lines(x=c(9,11),y=c(.85,.85),lwd=3)
text(10,.85,labels="Specificity",pos=3)

##absolute terms
par(fig=c(.5,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)

boxplot(resMat2[,4],at=1,xlim=c(0.5,11.5),ylim=c(0,1),width=2,axes=F,lwd=2,
  main="Discrimination of PAM50 Basal vs Luminal split"
)
#abline(h=0,lwd=3,col="lightgrey",lty=2)
boxplot(resMat2[,7],at=2,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,10],at=3,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,5],at=5,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,8],at=6,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,11],at=7,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,6],at=9,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,9],at=10,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,12],at=11,add=T,width=2,axes=F,lwd=2)

axis(1,at=c(1:3,5:7,9:11),lwd=2,las=2,cex=1.2,
  labels=c("Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify")
)
axis(2,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
mtext(side=2, "Absolute level",font=2,line=2.5,cex=1.2)
lines(x=c(1,3),y=c(.05,.05),lwd=3)
text(2,.05,labels="Accuracy",pos=3)

lines(x=c(5,7),y=c(.05,.05),lwd=3)
text(6,.05,labels="Sensitivity",pos=3)

lines(x=c(9,11),y=c(.05,.05),lwd=3)
text(10,.05,labels="Specificity",pos=3)

dev.off()

##as tiff 
tiff(paste0(HOME,"/fig2_top100percentBySd_confusionStats_basalVsLuminalSplitIn100randomSets.tiff"),width=12*500,height=12*500,units="px",res=500,compression="lzw")
par(fig=c(0,.5,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=F)

boxplot(resMat2[,4]-resMat2[,1],at=1,xlim=c(0.5,11.5),ylim=c(-1.,1.1),width=2,axes=F,lwd=2,
  main="Discrimination of PAM50 Basal vs Luminal split"
)
abline(h=0,lwd=3,col="lightgrey",lty=2)
boxplot(resMat2[,7]-resMat2[,1],at=2,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,10]-resMat2[,1],at=3,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,5]-resMat2[,2],at=5,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,8]-resMat2[,2],at=6,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,11]-resMat2[,2],at=7,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,6]-resMat2[,3],at=9,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,9]-resMat2[,3],at=10,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,12]-resMat2[,3],at=11,add=T,width=2,axes=F,lwd=2)

axis(1,at=c(1:3,5:7,9:11),lwd=2,las=2,cex=1.2,
  labels=c("Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify")
)
axis(2,at=round(seq(-1,1,length.out=5),2),lwd=2,las=1,cex=1.2)
mtext(side=2, "Relative to unadjusted data",font=2,line=2.5,cex=1.2)
lines(x=c(1,3),y=c(.85,.85),lwd=3)
text(2,.85,labels="Accuracy",pos=3)

lines(x=c(5,7),y=c(.85,.85),lwd=3)
text(6,.85,labels="Sensitivity",pos=3)

lines(x=c(9,11),y=c(.85,.85),lwd=3)
text(10,.85,labels="Specificity",pos=3)

##absolute terms
par(fig=c(.5,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)

boxplot(resMat2[,4],at=1,xlim=c(0.5,11.5),ylim=c(0,1),width=2,axes=F,lwd=2,
  main="Discrimination of PAM50 Basal vs Luminal split"
)
#abline(h=0,lwd=3,col="lightgrey",lty=2)
boxplot(resMat2[,7],at=2,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,10],at=3,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,5],at=5,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,8],at=6,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,11],at=7,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,6],at=9,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,9],at=10,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,12],at=11,add=T,width=2,axes=F,lwd=2)

axis(1,at=c(1:3,5:7,9:11),lwd=2,las=2,cex=1.2,
  labels=c("Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify",
  "Adjusted",
  "Beta>0.3",
  "InfiniumPurify")
)
axis(2,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
mtext(side=2, "Absolute level",font=2,line=2.5,cex=1.2)
lines(x=c(1,3),y=c(.05,.05),lwd=3)
text(2,.05,labels="Accuracy",pos=3)

lines(x=c(5,7),y=c(.05,.05),lwd=3)
text(6,.05,labels="Sensitivity",pos=3)

lines(x=c(9,11),y=c(.05,.05),lwd=3)
text(10,.05,labels="Specificity",pos=3)

dev.off()

rm(resMat2)

################################################################################  
###Do plot with one of the random 500 iterations

##2 figure panels
  ##1. heatmap uncorrected
  ##2. heatmap corrected
  ##3. basal/luminal v 2-group with 1 placeholder

##choose one of the random iters by accuracy
  ##best change pre-post?!
load(file=paste0(HOME,"/top100percentBySd_basalVsLuminalSplitIn100randomSets_usedProbeSets.RData"))

load(file=paste0(HOME,"/top100percentBySd_basalVsLuminalSplitIn100randomSets_usedRngSeeds.RData"))

iii<-which.max( resMat[,4]-resMat[,1] )
s_list<-s_list[[iii]]
iii<-p_list[[iii]]

##do multicore
no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

clusterEvalQ(cl, {
  library("flexmix")
})

betaRun<-cbind(seed=s_list,betaData[iii,])
betaNames<-colnames(betaData)
b<-parRapply(cl = cl, betaRun, adjustBeta,purity=fracTum,snames=betaNames,seed=TRUE)

rm(iii,p_list,s_list,betaRun,betaNames)

##do heatmap

##adjusted
b1<-do.call("rbind",lapply(b,function(x) x$y.tum))
##unadjusted
b2<-do.call("rbind",lapply(b,function(x) x$y.orig))

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(b1) ),method="ward.D"),2)
c2<-unique(c1[hclust( as.dist( 1-cor(b1) ),method="ward.D")$order])
r1<-hclust( dist(b2),method="ward.D") ##do rowclusters based on unadj data
c3<-hclust( as.dist( 1-cor(b1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(b2) ),method="ward.D"),2)
c5<-hclust( as.dist( 1-cor(b2) ),method="ward.D")

sample_anno<-data.frame(adj500=as.character(c1),
  unadj500=as.character(c4),
  PAM50=sampleAnno$pam50.full
  #TNBC=tnbcClass$TNBCtype
  #umap=factor(clusters.umap[samples_use,"class"]),
  #hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-colnames(betaData)
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj500=c("1"="#E41A1C","2"="#377EB8"),
    adj500=c("1"="#E41A1C","2"="#377EB8"),
    #adj5000=c("a"="red","b"="pink","c"="darkgreen","d"="orange","e"="grey"),
    #umap = c("1" = "#5977ff", "2" = "#f74747"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white")
    #TNBC = c("BL1"="#E41A1C","BL2"="#377EB8","IM"="#4DAF4A","LAR"="#984EA3","M"="#FF7F00","MSL"="#FFFF33","NA"="#666666","UNS"="#A65628")
    #,hrd3 = c("[0.0,0.2)" ="#FEE0D2" , "[0.2,0.7)" ="#FC9272" ,"[0.7,1.0]"="#EF3B2C" )
  )

tiff(paste0(HOME,"/fig2_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(b1,cluster_rows = r1,cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",fontsize=18,cutree_cols=2
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/fig2_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta_annotations.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(b1,cluster_rows = r1,cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",fontsize=18,cutree_cols=2
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=TRUE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(b2) ),method="ward.D"),2)
c2<-unique(c1[hclust( as.dist( 1-cor(b2) ),method="ward.D")$order])
r1<-hclust( dist(b2),method="ward.D")
c3<-hclust( as.dist( 1-cor(b2) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(b1) ),method="ward.D"),2)

sample_anno<-data.frame(unadj500=as.character(c1),
  adj500=as.character(c4),
  PAM50=sampleAnno$pam50.full
  #TNBC=tnbcClass$TNBCtype
  #umap=factor(clusters.umap[samples_use,"class"]),
  #hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-colnames(betaData)
sample_anno<-sample_anno[,ncol(sample_anno):1]

tiff(paste0(HOME,"/fig2_random500_heatmap_noAnno_pear_eucl_unadjClust_unadjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(b2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",cutree_cols=2,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

##do composite plot
a2<-image_read(paste0(HOME,"/fig2_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta.tiff"))
a1<-image_read(paste0(HOME,"/fig2_random500_heatmap_noAnno_pear_eucl_unadjClust_unadjBeta.tiff"))

a3<-image_read(paste0(HOME,"/fig2_top100percentBySd_confusionStats_basalVsLuminalSplitIn100randomSets.tiff"))
a3<-image_crop(a3,"6000x3000")

a11<-image_read(paste0(HOME,"/fig2_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta_annotations.tiff"))

a11<-image_crop(a11,"1000x6000+4300")

tiff(paste0(HOME,"/fig2_random500_heatmap_noAnno_white.tiff"),width=.2*500,height=12*500,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a9<-image_read(paste0(HOME,"/fig2_random500_heatmap_noAnno_white.tiff"))

out<-image_append(c(a9,
  a1,
  a9,
  a2,
  a9,
  a11
  ),stack = F)
out<-image_scale(out,"6000x")

out<-image_append(c(out,a3),stack=T)

image_write(out, path = paste0(HOME,"/fig2_random500_heatmap_noAnno_unadj_adj_combined.tiff"), format = "tiff")

rm(a1,a2,a3,a9,a11,out)
rm(c1,c2,c3,c4,c5,r1,b,b1,b2,cl)

gc()

################################################################################
################################################################################
##Create sFig image - With infiniumPurify added

##Panel 1 - Correction top5k, unadj clust 3-panel                             ##

##check stats
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))
temp4<-beta_norm[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))
temp5<-InfiniumPurify(tumor.data=temp3,normal.data=temp4,purity=fracTum)


table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp3,1,function(x) sum(is.na(x))))
#   0
#5000

table(annoObj[rownames(temp1),"chr"])
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4 
 #  557   272   213   255   156   144    98   169   272    69   205   365   130    43    54   248   192 
 # chr5  chr6  chr7  chr8  chr9 
 #  370   392   386   353    57 

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
r1<-hclust( dist(temp1),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp3) ),method="ward.D"),5)
c5<-cutree( hclust( as.dist( 1-cor(temp5) ),method="ward.D"),5)

sample_anno<-data.frame(adj5000=as.character(c1),
  infPur5000=as.character(c5),
  unadj5000=as.character(c4),
  ER=sampleAnno$ER,
  PR=sampleAnno$PR,
  HER2=sampleAnno$HER2,
  TNBC=as.character(as.integer(sampleAnno$TNBC)),
  PAM50=sampleAnno$pam50.full,stringsAsFactors=FALSE
  )
rownames(sample_anno)<-colnames(temp1)
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    infPur5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    ER = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    PR = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    HER2 = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A",
      "NA"="white"
      ),
    TNBC = c("1"="black","0"="lightgrey"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white"),stringsAsFactors=FALSE
  )

tiff(paste0(HOME,"/fig3_top5k_heatmap_pear_eucl_adjClust_infiniumPurifyBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp5,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="top 5000 by sd, \"infiniumPurify\", adj clust , pearCol/euclRow",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_infiniumPurifyBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp5,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="top 5000 by sd, \"infiniumPurify\", adj clust , pearCol/euclRow",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_unadjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp3,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="top 5000 by sd, unadj data, adj clust , pearCol/euclRow",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="top 5000 by sd, adj data, adj clust , pearCol/euclRow",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_InferredNormalBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="top 5000 by sd, \"inferred normal\", adj clust , pearCol/euclRow",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/fig3_splot_top5k_heatmap_whitePadding.tiff"),width=4115/500,height=250,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a6<-image_read(paste0(HOME,"/fig3_splot_top5k_heatmap_whitePadding.tiff"))

##cut away anno legend from 2/3 paste together to form image..
a1<-image_read(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_InferredNormalBeta.tiff"))
a1<-image_crop(a1,"4115x6375+0+125")
a1<-image_append(c(a6,a1),stack = T)
a1<-image_annotate(a1, "Inferred normal", size = 200, gravity = "north", color = "black")

a2<-image_read(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_unadjBeta.tiff"))
a2<-image_crop(a2,"4115x6375+0+125")
a2<-image_append(c(a6,a2),stack = T)
a2<-image_annotate(a2, "Unadjusted beta", size = 200, gravity = "north", color = "black")

a11<-image_read(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_infiniumPurifyBeta.tiff"))
a11<-image_crop(a11,"4115x6375+0+125")
a11<-image_append(c(a6,a11),stack = T)
a11<-image_annotate(a11, "InfiniumPurify", size = 200, gravity = "north", color = "black")

a3<-image_read(paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_adjBeta.tiff"))
a3<-image_crop(a3,"5000x6375+0+125")
a3<-image_append(c(a6,a3),stack = T)
a3<-image_annotate(a3, "Adjusted beta", size = 200, gravity = "north", color = "black")

tiff(paste0(HOME,"/fig3_top5k_heatmap_whitePadding.tiff"),width=.25*500,height=6375/500,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a4<-image_read(paste0(HOME,"/fig3_top5k_heatmap_whitePadding.tiff"))

tiff(paste0(HOME,"/fig3_top5k_heatmap_whitePadding2.tiff"),width=4150,height=.1*500,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a5<-image_read(paste0(HOME,"/fig3_top5k_heatmap_whitePadding2.tiff"))

out<-image_append(c(a4,a1,a4,
  a2,a4,a11,a4,
  a3
  ),stack = F)
out<-image_scale(out,"4150x")

out<-image_append(c(a5,out
  ),stack = T)

image_write(out, path = paste0(HOME,"/fig3_splot_top5k_heatmap_pear_eucl_adjClust_combined.tiff"), format = "tiff")

rm(a1,a2,a3,a4,a5,a11,out)
rm(c1,c2,c3,c4,r1,temp1,temp2,temp3,temp4)

gc()

################################################################################ 
################################################################################
##Create Fig 3 panels and image - Real life correction example

##Panel 1 - Correction top5k, unadj clust 3-panel                             ##

##check stats
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))
temp4<-beta_norm[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))
#temp5<-InfiniumPurify(tumor.data=temp3,normal.data=temp4,purity=fracTum)

table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp3,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp4,1,function(x) sum(is.na(x))))
#    0    1    2    3    4    5    6    7    8   10   11   12   17   40   60 
# 4085  470  200  110   59   31   16    8    7    3    6    2    1    1    1 

table(annoObj[rownames(temp1),"chr"])
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4 
 #  557   272   213   255   156   144    98   169   272    69   205   365   130    43    54   248   192 
 # chr5  chr6  chr7  chr8  chr9 
 #  370   392   386   353    57 

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-hclust( as.dist( 1-cor(temp4) ),method="ward.D")
r1<-hclust( dist(temp1),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp3) ),method="ward.D"),5)

sample_anno<-data.frame(adj5000=as.character(c1),
  unadj5000=as.character(c4),
  ER=sampleAnno$ER,
  PR=sampleAnno$PR,
  HER2=sampleAnno$HER2,
  TNBC=as.character(as.integer(sampleAnno$TNBC)),
  PAM50=sampleAnno$pam50.full,stringsAsFactors=FALSE
  )
rownames(sample_anno)<-colnames(temp1)
sample_anno<-sample_anno[,ncol(sample_anno):1]


my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    ER = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    PR = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    HER2 = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A",
      "NA"="white"
      ),
    TNBC = c("1"="black","0"="lightgrey"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080"),stringsAsFactors=FALSE
  )

#tiff(paste0(HOME,"/fig3_top5k_heatmap_pear_eucl_adjClust_unadjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
tiff(paste0(HOME,"/plos_figure3_panel_1.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
pheatmap(temp3,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="",cutree_cols=5,fontsize = 12,cellwidth = .9
  ,annotation_col=sample_anno,annotation_colors=my_colour#,annotation_row=r_anno
)
dev.off()

#tiff(paste0(HOME,"/fig3_top5k_heatmap_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
tiff(paste0(HOME,"/plos_figure3_panel_2.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="",cutree_cols=5,fontsize = 12,cellwidth = .9
  ,annotation_col=sample_anno,annotation_colors=my_colour#,annotation_row=r_anno
)
dev.off()

#tiff(paste0(HOME,"/fig3_top5k_heatmap_pear_eucl_adjClust_InferredNormalBeta.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
tiff(paste0(HOME,"/plos_figure3_panel_3.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="",cutree_cols=5,fontsize = 12,cellwidth = .9
  ,annotation_col=sample_anno,annotation_colors=my_colour#,annotation_row=r_anno
)
dev.off()

##Fix panel

##ONE
fig<-image_read(paste0(HOME,"/plos_figure3_panel_1.tiff"))

fig<-image_crop(fig,"2470x4500+700")
fig<-image_annotate(fig,"Unadjusted", gravity="north", size = 76, color = "black",font="sans",weight=700)
#fig<-image_scale(fig,"750x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_1_adj.tiff"))

##TWO
fig<-image_read(paste0(HOME,"/plos_figure3_panel_2.tiff"))

fig<-image_crop(fig,"2470x4500+700")
fig<-image_annotate(fig,"Adjusted", gravity="north", size = 76, color = "black",font="sans",weight=700)
#fig<-image_scale(fig,"750x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_2_adj.tiff"))

##THREE
fig<-image_read(paste0(HOME,"/plos_figure3_panel_3.tiff"))

fig<-image_crop(fig,"2470x4500+700")
fig<-image_annotate(fig,"Inferred normal", gravity="north", size = 76, color = "black",font="sans",weight=700)
fig<-image_annotate(fig, "A", size = 76, gravity="northwest", color = "black",font="sans",weight=700)
#fig<-image_scale(fig,"750x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_3_adj.tiff"))

##ANNOTATIONS
fig<-image_read(paste0(HOME,"/plos_figure3_panel_1.tiff"))
fig<-image_crop(fig,"700x4500+3170")
#fig<-image_scale(fig,"750x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_4_adj.tiff"))

##MERGE
fig1<-image_read(paste0(HOME,"/plos_figure3_panel_1_adj.tiff"))
fig2<-image_read(paste0(HOME,"/plos_figure3_panel_2_adj.tiff"))
fig3<-image_read(paste0(HOME,"/plos_figure3_panel_3_adj.tiff"))
fig4<-image_read(paste0(HOME,"/plos_figure3_panel_4_adj.tiff"))

fig<-image_append(c(fig3,fig1,fig2,fig4),stack = F)
fig<-image_scale(fig,"2250x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_A.tiff"))

##NEW PANEL

##Panel 2 - beta density pre/post                                             ##
temp4<-do.call("rbind",lapply(res,function(x) x$y.tum))
#rownames(temp4)<-rownames(testDat2)
temp5<-do.call("rbind",lapply(res,function(x) x$y.orig))
#rownames(temp5)<-rownames(testDat2)

table(apply(temp4,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp5,1,function(x) sum(is.na(x))))
#   0
#5000

quantile(testDat2)
#    0%   25%   50%   75%  100% 
# 0.000 0.152 0.468 0.741 1.000 

quantile(temp4)
#    0%   25%   50%   75%  100% 
# 0.000 0.061 0.529 0.916 1.000 

##Plot histogram of betas before and after correction
tiff(paste0(HOME,"/plos_figure3_panel_b.tiff"),width=7.5*300,height=7.5*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
#tiff(paste0(HOME,"/fig3_top5k_betaDistribution_tumors_beforeAfterCorrection.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")

#par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot(1,xlim=range(round(density(temp4)$x,1)),ylim=range(round(density(temp4)$y,1))+c(0,.5),type="n",las=1,axes=F,
  xlab="beta",ylab="density"
)
lines(density(temp4),col=2,lwd=4)
lines(density(testDat2),col=1,lwd=4)
axis(1,lwd=3,las=1,at=seq(0,1,by=.2))
axis(2,lwd=3,las=1)
legend("topright",legend=c("Unadjusted beta","Adjusted beta"),col=c(1,2),lwd=4,bty="n")
dev.off()

fig<-image_read(paste0(HOME,"/plos_figure3_panel_b.tiff"))
fig<-image_annotate(fig, "B", size = 76, gravity="northwest", color = "black",font="sans",weight=700)
fig<-image_scale(fig,"400x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))

# ##Plot histogram of betas before and after correction
# tiff(paste0(HOME,"/fig3_top5k_tumorFracVsMeanCorrection.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# par(font=2,font.axis=2,font.lab=2,font.sub=2)
# plot(fracTum,
#     apply(abs(temp4-temp5),2,mean),
#   pch=16,cex=.5,xlim=c(0,1),ylim=c(0,.5),
#   main="Average beta correction vs global correlation",
#   xlab="tumor fraction",
#   ylab="mean absolute beta difference pre-post correction",
#   axes=F#,type="n"
# )
# axis(1,lwd=2,las=1)
# axis(2,lwd=2,las=1)
# dev.off()

rm(temp4,temp5)

##Panel 3 - correl inferred vs actual                                         ##
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))

temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))

ff<-intersect( rownames(temp2) , rownames(beta_norm) )

length(ff)
#[1] 5000

#plot( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]) )

table(apply(beta_norm_adj[ff,],1,function(z) any(is.na(z))))
# FALSE  TRUE 
#  4085   915 

table(apply(beta_norm_adj[ff,],1,function(z) sum(is.na(z))))
#    0    1    2    3    4    5    6    7    8   10   11   12   17   40   60 
# 4085  470  200  110   59   31   16    8    7    3    6    2    1    1    1 

plot( rowMeans(temp2[ff,]),rowMeans(beta_norm_adj[ff,],na.rm=TRUE) )
dev.off()

##will change slightly if rerun - not deterministic..
cor( rowMeans(temp2[ff,]),rowMeans(beta_norm_adj[ff,],na.rm=TRUE) )
#[1] 0.9250216
(cor( rowMeans(temp2[ff,]),rowMeans(temp1[ff,]),method="spe" ))
#[1] 0.133569
(cor( rowMeans(temp2[ff,]),rowMeans(beta_norm_adj[ff,],na.rm=TRUE),method="spe" ))
#[1] 0.8848514
(cor( rowMeans(temp2[ff,]),rowMeans(temp1[ff,]),method="pe" ))
#[1] 0.01715253
(sf<-cor( rowMeans(temp2[ff,]),rowMeans(beta_norm_adj[ff,],na.rm=TRUE),method="pe" ))
#[1] 0.9250216
(fs<-cor.test( rowMeans(temp2[ff,]),rowMeans(beta_norm_adj[ff,],na.rm=TRUE),method="pe" )$p.value)
#[1] 0

length(ff)
#[1] 5000

temp4<-beta_norm_adj[ff,]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))

pdf(paste0(HOME,"/fig3_top5kBySd_betaNormals_inferredVsActual.pdf"),width=8,height=8,useDingbats=F)
par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot( rowMeans(temp2[ff,]),rowMeans(temp4[ff,],na.rm=TRUE),pch=16
  ,main=""
  ,xlab="mean inferred normal beta",ylab="mean normal beta 450k GSE67919"
  ,type="n",las=1,axes=F,xlim=c(0,1),ylim=c(0,1)
 )
points(rowMeans(temp2[ff,]),rowMeans(temp4[ff,]),pch=16)
text(.1,.9,paste0("r=",round(sf,2)," | p<2.2e-16 \n",length(ff)," CpGs"))
abline(lm(rowMeans(temp4[ff,])~rowMeans(temp2[ff,])),lwd=2,col=2)
axis(1,lwd=2,las=1,at=seq(0,1,by=.2))
axis(2,lwd=2,las=1,at=seq(0,1,by=.2))
dev.off()

tiff(paste0(HOME,"/fig3_top5kBySd_betaNormals_inferredVsActual.tiff"),,width=8*500,height=8*500,units="px",res=500,compression="lzw")
par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot( rowMeans(temp2[ff,]),rowMeans(temp4[ff,],na.rm=TRUE),pch=16
  ,main=""
  ,xlab="mean inferred normal beta",ylab="mean normal beta 450k GSE67919"
  ,type="n",las=1,axes=F,xlim=c(0,1),ylim=c(0,1)
 )
points(rowMeans(temp2[ff,]),rowMeans(temp4[ff,]),pch=16)
text(.1,.9,paste0("r=",round(sf,2)," | p<2.2e-16 \n",length(ff)," CpGs"))
abline(lm(rowMeans(temp4[ff,])~rowMeans(temp2[ff,])),lwd=2,col=2)
axis(1,lwd=2,las=1,at=seq(0,1,by=.2))
axis(2,lwd=2,las=1,at=seq(0,1,by=.2))
dev.off()

rm(ff,fs)

##combine
a1<-image_read(paste0(HOME,"/fig3_top5k_betaDistribution_tumors_beforeAfterCorrection.tiff"))
a1<-image_scale(a1,"2000x")

a2<-image_read(paste0(HOME,"/fig3_top5kBySd_betaNormals_inferredVsActual.tiff"))
a2<-image_scale(a2,"2000x")

out<-image_append(c(a1,a2
  ),stack = T)
out<-image_scale(out,"x2000")

image_write(out, path = paste0(HOME,"/fig3_top5kBySd_betaDistAndInfVAct_combined.tiff"), format = "tiff")

#fig<-image_annotate(fig, "B", size = 38, location=paste0("+",2250*.3+40,"+",40), color = "black",font="sans",weight=700)


##Panel 4 - boxplot purity 5-group adj/unadj                                  ##
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))

temp2<-do.call("rbind",lapply(res,function(x) x$y.orig))

c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-cutree( hclust( as.dist( 1-cor(temp2) ),method="ward.D"),5)

tiff(paste0(HOME,"/fig3_top5kBySd_clusterVpurity_adj.tiff"),,width=8*500,height=8*500,units="px",res=500,compression="lzw")
par(font=2,font.axis=2,font.lab=2,font.sub=2)#,fig=c(0.25,.75,0,1),mar=c(5.1,4.1,0.1,2.1))
boxplot(fracTum~c1,col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
  ,main=""
  ,xlab="adjusted beta - 5 group",ylab="ABSOLUTE purity"
  ,type="n",las=1,axes=F,xlim=c(0,6),ylim=c(0,1)
  )
axis(1,lwd=0,las=1,at=seq(1,5,by=1),font=2)
axis(2,lwd=2,las=1,at=seq(0,1,by=.2),font=2)
text(2,.05,paste0("p(kw-test)=",signif(kruskal.test(fracTum~c1)$p.value,2)))
dev.off()

tiff(paste0(HOME,"/fig3_top5kBySd_clusterVpurity_unadj.tiff"),,width=8*500,height=8*500,units="px",res=500,compression="lzw")
par(font=2,font.axis=2,font.lab=2,font.sub=2)#,fig=c(0.25,.75,0,1),mar=c(5.1,4.1,0.1,2.1))
boxplot(fracTum~c2,col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
  ,main=""
  ,xlab="unadjusted beta - 5 group",ylab="ABSOLUTE purity"
  ,type="n",las=1,axes=F,xlim=c(0,6),ylim=c(0,1)
  )
axis(1,lwd=0,las=1,at=seq(1,5,by=1),font=2)
axis(2,lwd=2,las=1,at=seq(0,1,by=.2),font=2)
text(2,.05,paste0("p(kw-test)=",signif(kruskal.test(fracTum~c2)$p.value,2)))
dev.off()

a1<-image_read(paste0(HOME,"/fig3_top5kBySd_clusterVpurity_adj.tiff"))
# a1<-image_crop(a1,"3000x2000+1000")
# a1<-image_crop(a1,"2000x2000")
a1<-image_scale(a1,"x2000")

a2<-image_read(paste0(HOME,"/fig3_top5kBySd_clusterVpurity_unadj.tiff"))
# a2<-image_crop(a2,"3000x2000+1000")
# a2<-image_crop(a2,"2000x2000")
a2<-image_scale(a2,"x2000")

out<-image_append(c(a2,a1
  ),stack = T)
out<-image_scale(out,"x2000")

image_write(out, path = paste0(HOME,"/fig3_top5kBySd_clusterVpurity_combined.tiff"), format = "tiff")

##Panel 5 - barplot PAM50

tiff(paste0(HOME,"/fig3_top5kBySd_barplotPam50_adj.tiff"),,width=8*500,height=8*500,units="px",res=500,compression="lzw")
par(font=2,font.axis=2,font.lab=2,font.sub=2)
tt<-table(sampleAnno$pam50.full,factor(c1))
barplot(t( t(tt)/colSums(tt)),col=c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white"),
  ,main="",legend.text=F
  ,xlab="adjusted beta - 5 group",ylab="Fraction of cluster"
  ,las=1,axes=F,xlim=c(0,6),ylim=c(0,1)
  )
#axis(1,lwd=0,las=1,at=seq(1,5,by=1),font=2)
axis(2,lwd=2,las=1,at=seq(0,1,by=.2),font=2)
dev.off()

tiff(paste0(HOME,"/fig3_top5kBySd_barplotPam50_unadj.tiff"),,width=8*500,height=8*500,units="px",res=500,compression="lzw")
par(font=2,font.axis=2,font.lab=2,font.sub=2)
tt<-table(sampleAnno$pam50.full,factor(c2))
barplot(t( t(tt)/colSums(tt)),col=c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white"),
  ,main="",legend.text=T
  ,xlab="unadjusted beta - 5 group",ylab="Fraction of cluster"
  ,las=1,axes=F,xlim=c(0,6),ylim=c(0,1)
  )
#axis(1,lwd=0,las=1,at=seq(1,5,by=1),font=2)
axis(2,lwd=2,las=1,at=seq(0,1,by=.2),font=2)
dev.off()

a1<-image_read(paste0(HOME,"/fig3_top5kBySd_barplotPam50_adj.tiff"))
# a1<-image_crop(a1,"3000x2000+1000")
# a1<-image_crop(a1,"2000x2000")
a1<-image_scale(a1,"x2000")

a2<-image_read(paste0(HOME,"/fig3_top5kBySd_barplotPam50_unadj.tiff"))
# a2<-image_crop(a2,"3000x2000+1000")
# a2<-image_crop(a2,"2000x2000")
a2<-image_scale(a2,"x2000")

out<-image_append(c(a2,a1
  ),stack = T)
out<-image_scale(out,"x2000")

image_write(out, path = paste0(HOME,"/fig3_top5kBySd_barplotPam50_combined.tiff"), format = "tiff")

##combine again
a1<-image_read(paste0(HOME,"/fig3_top5kBySd_clusterVpurity_combined.tiff"))
# a1<-image_crop(a1,"3000x2000+1000")
# a1<-image_crop(a1,"2000x2000")
a1<-image_scale(a1,"x2000")

a2<-image_read(paste0(HOME,"/fig3_top5kBySd_barplotPam50_combined.tiff"))
# a2<-image_crop(a2,"3000x2000+1000")
# a2<-image_crop(a2,"2000x2000")
a2<-image_scale(a2,"x2000")

out<-image_append(c(a1,a2
  ),stack = F)
#out<-image_scale(out,"x2000")

image_write(out, path = paste0(HOME,"/fig3_top5kBySd_clusterStat_combined.tiff"), format = "tiff")


##MERGE
image_write(image_blank(width=44, height=400, color = "white"),paste0(HOME,"/plos_figure3_panel_b_blank.tiff"),format="tiff",density=300)

fig1<-image_read(paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))
fig2<-image_read(paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))
fig3<-image_read(paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))
fig4<-image_read(paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))
fig5<-image_read(paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))
fig6<-image_read(paste0(HOME,"/plos_figure3_panel_b_adj.tiff"))
fig7<-image_read(paste0(HOME,"/plos_figure3_panel_b_blank.tiff"))

fig<-image_append( c( image_append(c(fig1,fig7,fig2,fig7,fig3,fig7,fig7),stack = F),
  image_append(c(fig4,fig7,fig5,fig7,fig6,fig7,fig7),stack = F) ),stack = T )
fig<-image_scale(fig,"2250x")
image_write(fig,paste0(HOME,"/plos_figure3_panel_BG.tiff"))

##MERGE FINAL
fig1<-image_read(paste0(HOME,"/plos_figure3_panel_A.tiff"))
fig2<-image_read(paste0(HOME,"/plos_figure3_panel_BG.tiff"))
fig<-image_append( c(fig1,fig2),stack=T)
fig<-image_scale(fig,"2250x")
image_write(fig,paste0(HOME,"/plos_figure3_final.tiff"))

################################################################################
################################################################################
##Create Fig 4 panels and image

##Panel 1 - Autosomes vs X?

all.equal(names(fracTum),colnames(betaOrig))
#[1] TRUE

all.equal(rownames(annoObj),rownames(betaOrig))
#[1] TRUE

all.equal(rownames(annoObj),rownames(betaNorm))
#[1] TRUE

all(rownames(annoObj) %in% rownames(beta_norm))
#[1] TRUE

table(annoObj$chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
# 41098 21206 25633 21624 10540 13204 13119 18627 24356  5280 21984 30550  9354 
# chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY 
#  3378  7222 22315 17646 21279 31170 25024 18194  8605  9718   242 

table(annoObj$illuminaCpG_CpH_Probe)
#     cg     ch 
# 418858   2510 

betaData<-betaOrig[ !annoObj$chr %in% c("chrY") &
   annoObj$illuminaCpG_CpH_Probe == "cg"
,]

betaData2<-betaAdj[ !annoObj$chr %in% c("chrY") &
   annoObj$illuminaCpG_CpH_Probe == "cg"
,]

betaNorm2<-betaNorm[ !annoObj$chr %in% c("chrY") &
   annoObj$illuminaCpG_CpH_Probe == "cg"
,]

beta_norm2<-beta_norm_adj[rownames(annoObj),]
beta_norm2<-beta_norm2[ !annoObj$chr %in% c("chrY") &
   annoObj$illuminaCpG_CpH_Probe == "cg"
,]

annoObj2<-annoObj[ !annoObj$chr %in% c("chrY") &
   annoObj$illuminaCpG_CpH_Probe == "cg"
,]

all.equal(rownames(annoObj2),rownames(betaData))
#[1] TRUE

all.equal(rownames(annoObj2),rownames(betaData2))
#[1] TRUE

##panel 1 - X-beta in normals
tiff(paste0(HOME,"/plos_figure4_panel_a.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,6.1,4.1,0.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=30,new=F)
#tiff(paste0(HOME,"/fig4_xChromAdjStats_panel1.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
#par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot(density(as.vector(betaNorm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 & annoObj2$saxonovClass == "HCG",]))
  ,main="",type="n"
  ,xlab="X chromosome HCG promoter methylation (N=3623 CpGs)",ylab="Density"
  ,las=1,axes=F,xlim=c(0,1),ylim=c(0,2.5)
  )
abline(v=c(0,.25,.5,.75,1),lty=c(1,2,1,2,1),lwd=6,col="lightgrey")
lines(density((betaNorm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 &  annoObj2$saxonovClass == "HCG",])),lwd=6,col="#E41A1C",lty=1)
lines(density((beta_norm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 & annoObj2$saxonovClass == "HCG",]),na.rm=T),lwd=6,col="#377EB8",lty=1)
axis(1,lwd=6,las=1,at=seq(0,1,by=.25),font=2)
axis(2,lwd=6,las=1,font=2,at=seq(0,2,by=1))
legend("topright",legend=c("inferred normal","GSE67919 normal"),bty="n",col=c("#E41A1C","#377EB8"),lwd=6)

text(.85,.075,labels="Hyper > 0.75")
text(.5,.075,labels="0.25 < Xi < 0.75")
text(.15,.075,labels="Hypo < 0.25")

ai<-table(apply(betaNorm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 & annoObj2$saxonovClass == "HCG",],1,function(x) {
  if( median(x)<.25) { 
    return("hypo") 
  }
  if( median(x)>.75) { 
    return("hyper") 
  } else {
    return("Xi")
  }
}))
text(.40,1.25,paste("inferred normal",paste(paste(names(ai),round(ai/sum(ai),2),sep=": ")," (N=",ai,")",sep="",collapse="\n"),sep="\n"))

bi<-table(apply(beta_norm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 & annoObj2$saxonovClass == "HCG",],1,function(x) {
  if( median(x,na.rm=T)<.25) { 
    return("hypo") 
  }
  if( median(x,na.rm=T)>.75) { 
    return("hyper") 
  } else {
    return("Xi")
  }
}))
text(.40,.75,paste("GSE67919 normal",paste(paste(names(bi),round(bi/sum(bi),2),sep=": ")," (N=",bi,")",sep="",collapse="\n"),sep="\n"))

ai<-apply(betaNorm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 &  annoObj2$saxonovClass == "HCG",],1,function(x) {
  if( median(x)<.25) { 
    return("hypo") 
  }
  if( median(x)>.75) { 
    return("hyper") 
  } else {
    return("Xi")
  }
})
bi<-apply(beta_norm2[ annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1 &  annoObj2$saxonovClass == "HCG",],1,function(x) {
  if( median(x,na.rm=T)<.25) { 
    return("hypo") 
  }
  if( median(x,na.rm=T)>.75) { 
    return("hyper") 
  } else {
    return("Xi")
  }
})
names(ai)<-sub("_.+","",annoObj2$namePromOverlap[annoObj2$chr == "chrX"  & annoObj2$hasPromOverlap==1 & annoObj2$saxonovClass == "HCG"])
names(bi)<-sub("_.+","",annoObj2$namePromOverlap[annoObj2$chr == "chrX" & annoObj2$hasPromOverlap==1  & annoObj2$saxonovClass == "HCG"])
text(.40,.4,paste0("concordance=",round(100*sum(diag(table(ai,bi))/length(ai)),1),"%"))

##compare with published list of Xi-escape genes (N=75)
  ##https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5507-6#Sec22
library(readxl)
library(httr)
url1<-"https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-019-5507-6/MediaObjects/12864_2019_5507_MOESM6_ESM.xlsx"
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
df <- read_xlsx(tf,skip=1)
colnames(df)<-sub("\\W","",colnames(df))
str(df)
df<-df$GeneSymbol

text(.7,2.2,paste0("Katsir & Linial 2019 Xi-escape overlap"))
text(.65,2.1,paste0("Hypo: ",round(100 *length(intersect(sub("ENSG.+\\|","",names(ai)[ai=="hypo"]),df))/sum(ai=="hypo"),1),"%"),pos=4)
text(.65,2,paste0("Xi: ",round(100 *length(intersect(sub("ENSG.+\\|","",names(ai)[ai=="Xi"]),df))/sum(ai=="Xi"),1),"%"),pos=4)
text(.65,1.9,paste0("Hyper: ",round(100 *length(intersect(sub("ENSG.+\\|","",names(ai)[ai=="hyper"]),df))/sum(ai=="hyper"),1),"%"),pos=4)
dev.off()

##panel 2 - Beta dist in unadj vs adj
tiff(paste0(HOME,"/plos_figure4_panel_b.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,6.1,4.1,0.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=30,new=F)
# tiff(paste0(HOME,"/fig4_xChromAdjStats_panel2.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot(density((betaData[ annoObj2$chr == "chrX" ,]))
  ,main="",type="n"
  ,xlab="X chromosome HCG CpG methylation (N=6736 CpGs)",ylab="Density"
  ,las=1,axes=F,xlim=c(0,1),ylim=c(0,3.5)
  )
abline(v=c(0,.25,.5,.75,1),lty=c(1,2,1,2,1),lwd=6,col="lightgrey")
lines(density((betaData[ annoObj2$chr == "chrX" & annoObj2$saxonovClass == "HCG",])),lwd=6,col="#377EB8",lty=1)
lines(density((betaData2[ annoObj2$chr == "chrX" & annoObj2$saxonovClass == "HCG",])),lwd=6,col="#E41A1C",lty=1)
axis(1,lwd=6,las=1,at=seq(0,1,by=.25),font=2)
axis(2,lwd=6,las=1,font=2,at=seq(0,3.5,by=1))
legend("topright",legend=c("adjusted beta","unadjusted beta"),bty="n",col=c("#E41A1C","#377EB8"),lwd=6)
dev.off()

a1<-image_read(paste0(HOME,"/plos_figure4_panel_a.tiff"))
a1<-image_annotate(a1, "A", size = 152, location=paste0("+",80,"+",80), color = "black",font="sans",weight=700)
a1<-image_scale(a1,"1125x")

a2<-image_read(paste0(HOME,"/plos_figure4_panel_b.tiff"))
a2<-image_annotate(a2, "B", size = 152, location=paste0("+",80,"+",80), color = "black",font="sans",weight=700)
a2<-image_scale(a2,"1125x")

out<-image_append(c(a1,a2
  ),stack = T)
#out<-image_scale(out,"2075x")

image_write(out, path = paste0(HOME,"/plos_figure4_final.tiff"), format = "tiff")

################################################################################
################################################################################
##Create Fig 5 panels and image

##check stats
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))
temp4<-beta_norm[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))
#temp5<-InfiniumPurify(tumor.data=temp3,normal.data=temp4,purity=fracTum)

table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp3,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp4,1,function(x) sum(is.na(x))))
#    0    1    2    3    4    5    6    7    8   10   11   12   17   40   60 
# 4085  470  200  110   59   31   16    8    7    3    6    2    1    1    1 

table(annoObj[rownames(temp1),"chr"])
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4 
 #  557   272   213   255   156   144    98   169   272    69   205   365   130    43    54   248   192 
 # chr5  chr6  chr7  chr8  chr9 
 #  370   392   386   353    57 

##load TFmatrix
load(file=paste0(DATA,"/annotateFeatures/","object_noChromFilter_421368x340_encodeTfbsAnnotations.RData"))

##Which are most enriched
Bb<-(colSums(tfMat[rownames(temp1),])/5000) / (colSums(tfMat)/nrow(tfMat))

length(Bb)
#[1] 340
Bb[order(Bb)]
#     CREBBP      NR0B1      PLRG1      RBM17      SAFB2      SMAD2      SRSF9     ZNF507       ZZZ3       KAT8     ZBTB7B     ZNF574      ASH1L       HSF1      NR2C2 
# 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.06713245 0.13705050 0.18539133 0.19237169 0.19415643 0.19872149 0.20081080 0.21710674 
#     ZNF579      ZBTB1     SREBF1     TRIP13       ELK1       NFIB    SNRNP70       SIX5      KDM5A      GATA4      DEAF1      SIN3B     SREBF2     ZBTB11      NR2C1 
# 0.22766175 0.23672360 0.23876616 0.23952705 0.24018308 0.24151580 0.24498140 0.24709852 0.24828427 0.25047444 0.25397523 0.25969872 0.26174316 0.26271190 0.26572221 
#      U2AF2       ZHX1    NEUROD1       ZHX2      CLOCK      THAP1       TBX3     ZNF444      DNMT1   HNRNPUL1       MBD2        MYB      PYGO2       NFYA      MYBL2 
# 0.27880310 0.28204803 0.28830772 0.28850858 0.28921240 0.28949187 0.29517898 0.29550188 0.30533913 0.30559464 0.30862844 0.30958979 0.31270353 0.31423388 0.32330026 
#        FUS       IRF2      HNF1A      BRCA1      CCAR2      HCFC1      U2AF1      KDM4B      LCORL        PML     ZNF639   C11orf30       E2F8       IRF1       ETV4 
# 0.32596047 0.32611519 0.32848801 0.32898536 0.32931626 0.33040578 0.33291103 0.33520237 0.33662313 0.33798046 0.33844819 0.33921015 0.34382158 0.35437338 0.35511585 
#       E2F4      PTBP1    HNRNPH1      TAF9B     FIP1L1      RBM15      PHF20      KAT2B     ZNF407        SKI       IRF3      ZMIZ1      WHSC1     GABPB1       MTA1 
# 0.35780500 0.35875745 0.36422950 0.36442032 0.36534116 0.36587670 0.36655521 0.36780622 0.36941421 0.37071390 0.37568567 0.38083615 0.38132851 0.38179542 0.38290865 
#       ELF4    CREB3L1       KLF5      NCOA1     ZNF207       E2F7      CDC5L       BRD4     SETDB1       E2F1      ZFP91     NUFIP1      CREB1      SMAD1        RLF 
# 0.38797300 0.38909986 0.38990940 0.39028978 0.39051027 0.39457011 0.39613587 0.40117644 0.41338221 0.41387791 0.41486673 0.41502731 0.41669251 0.41963330 0.42127003 
#      NFRKB     ZNF217      SRSF4       RFX5   ARHGAP35        JUN      RBM25       E4F1      TFAP4       SAFB    ZKSCAN1      ZBTB2      RBM22       IRF5      SMAD5 
# 0.42145720 0.42219646 0.42396963 0.42488144 0.42501095 0.42594840 0.42834849 0.43060161 0.43493865 0.43503536 0.43573429 0.43670259 0.43908004 0.43926818 0.44368921 
#      EP400      ESRRA      GMEB1      FOXK2      STAT1       RFX1      NCOA4      KDM5B    SMARCB1     PRDM10      RBM14      NR2F1       AFF1      FOXP1       SOX6 
# 0.44640317 0.44813140 0.45032698 0.45173002 0.45359828 0.45420430 0.46413470 0.46627403 0.47124101 0.47341415 0.47557089 0.47939871 0.48088473 0.48179285 0.48276162 
#        SRF      PRPF4     ZNF282      STAT2       LEF1      NCOA6    ZNF512B       BRD9      RUNX1    TBL1XR1     ZBTB40       CBX1    SUPT20H     HNRNPL       AGO1 
# 0.48330666 0.48455683 0.48463567 0.48536416 0.48544309 0.48630790 0.48696004 0.48710303 0.48849378 0.48957568 0.49035112 0.49085250 0.49282807 0.49485700 0.49597099 
#        ZFX      SAP30      GATA1      HNF4G     PHF21A      GATA3       ETS1      TAF15       CBFB      RBM34      MEIS2     WRNIP1      NR3C1     NFE2L2       AGO2 
# 0.49602075 0.49626368 0.49807393 0.50588481 0.50711685 0.50897917 0.50907716 0.51353008 0.51901004 0.51935251 0.52027538 0.52229062 0.52271228 0.52627927 0.52650264 
#      HNF4A     ZNF184      HDAC3      MIER1      ZNF24       ESR1      HDAC1     ARID3A     TCF7L2      FOXM1     NFATC3        MNT      NCOR1      PCBP1     TRIM24 
# 0.52837288 0.52914905 0.53227587 0.53423514 0.53451230 0.53895024 0.53907693 0.54374937 0.54559886 0.54587365 0.54643747 0.54784986 0.54796694 0.54961043 0.54980076 
#       PBX3      XRCC5       ZEB2       RXRA     ZNF830       CBX3      NCOA2      FOXA2        RB1       CHD4       BCL3     POU2F2     CC2D1A       NFYB       BCOR 
# 0.55056470 0.55213586 0.55257825 0.55305616 0.55376918 0.55412778 0.55443158 0.55526783 0.55587254 0.56170818 0.56302834 0.56450568 0.56606952 0.56633212 0.56896550 
#       MTA3       MAFF      CTBP1      DACH1       PAX8      MEF2B      FOSL2       MYNN       UBTF       ATF4      FOSL1       SKIL    GATAD2A      RBM39     ZNF318 
# 0.57077197 0.57407874 0.57567810 0.57824516 0.57927373 0.57942877 0.57963301 0.58130559 0.58176059 0.58280498 0.58327341 0.58500024 0.59041019 0.59321741 0.59454052 
#      FOXA1    SMARCA5       THRA      NR2F2     POLR2G      SOX13    ZSCAN29      KLF16      NR2F6       TAL1     TRIM22       ELF1       RELB       ATF3       CREM 
# 0.59535401 0.59628146 0.59662726 0.59688715 0.60128020 0.60140879 0.60375298 0.60435861 0.60490377 0.60830791 0.61006726 0.61078765 0.61100223 0.61543114 0.61642872 
#      RAD51    ZNF280A     CHAMP1       ARNT       ATF7       TCF7      ARID2    SMARCE1      RCOR1        NBN     RBFOX2       MITF     HNRNPK        MYC      ZBTB5 
# 0.61938051 0.61965882 0.62006865 0.62136246 0.62159317 0.62610401 0.62825108 0.62838409 0.62886857 0.62934876 0.63092253 0.63243057 0.63372915 0.63601136 0.63602717 
#      TEAD4      GABPA     ZNF687       TAF7       NRF1     TRIM28      NCOA3      PCBP2     ZNF274     ZBTB7A       CUX1       USF2    HNRNPLL     ZNF384     ZNF592 
# 0.63729075 0.63808594 0.64854978 0.64950214 0.65374443 0.65754849 0.65808507 0.65818973 0.66226798 0.66265162 0.66730932 0.66924667 0.67016367 0.67052040 0.67054444 
#    GATAD2B       ATF2       DPF2    SMARCC2       MAFK      ZBED1       NFE2     STAT5A      MEF2C       SIX4      TCF12      SIRT6       ZEB1       EGR1      BACH1 
# 0.67273774 0.67323942 0.67335307 0.67446693 0.67637780 0.67851234 0.67891726 0.67911634 0.68462940 0.69020147 0.69215766 0.69238873 0.69687519 0.69893231 0.70068286 
#       USF1      ZMYM3       E2F6      TBX21       PAX5       MCM5     PKNOX1       IRF4      KDM1A       PHF8     ZNF143        MGA       NFIC     ZBTB8A       JUND 
# 0.70085438 0.70109671 0.70252585 0.70501763 0.70696968 0.71418305 0.71477906 0.71701527 0.71963309 0.72392683 0.72447910 0.72474927 0.72495307 0.72545423 0.72593732 
#       MCM7     ARID1B       REST      MLLT1     ZNF263      EHMT2       CHD1        MAX       MTA2       PHB2     ZNF316      IKZF1       HES1       ETV6        FOS 
# 0.73582496 0.73812925 0.74075681 0.74622611 0.74763123 0.74907472 0.75047775 0.75202803 0.75498291 0.76122183 0.76494957 0.76507236 0.77575421 0.77593472 0.77933333 
#    BHLHE40      HDAC2      IKZF2       MXI1     HMBOX1        SP1    L3MBTL2    SMARCA4      MEF2A       CBX5      STAT3     ZBTB33    CBFA2T2      CEBPB    CBFA2T3 
# 0.78240793 0.78322041 0.78362806 0.78387893 0.79634869 0.80841203 0.81263043 0.81560307 0.82205416 0.82269597 0.82366361 0.83177828 0.83293674 0.83415219 0.83461098 
#        EED        TBP      RUNX3      SIN3A        YY1     NFATC1       SMC3       SPI1       MCM2       TAF1      NANOG      EP300       BMI1      RAD21      NFXL1 
# 0.83467804 0.83843876 0.84030245 0.84081304 0.84161188 0.84581719 0.86460297 0.87965257 0.89414960 0.90166748 0.90318802 0.93102734 0.93600410 0.97573832 0.98196454 
#     BCL11A       JUNB      XRCC3       CTCF       BATF       EBF1     POLR2A      ASH2L      HDAC6      RBBP5      EWSR1      GATA2      SRSF7       RNF2      KDM4A 
# 0.98951389 1.00004460 1.01128320 1.01647820 1.01785634 1.03158479 1.04118682 1.04247583 1.05195895 1.07309225 1.10402096 1.10696576 1.17046667 1.18861042 1.23716004 
#        ATM       ZNF8       MCM3      COPS2       CHD7      KAT2A      SUZ12       EZH2       CBX8       CBX2 
# 1.28990204 1.30512212 1.32193882 1.45299310 1.47222276 2.02580769 2.12930064 2.21950214 2.41570516 2.61349032 
rm(Bb)
##Seems very driven by CGIs

##Panel 1 - As heatmap

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-hclust( as.dist( 1-cor(temp3) ),method="ward.D")
r1<-hclust( dist(temp1),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp3) ),method="ward.D"),5)

sample_anno<-data.frame(adj5000=as.character(c1),
  #unadj5000=as.character(c4),
  #ER=sampleAnno$ER,
  #PR=sampleAnno$PR,
  #HER2=sampleAnno$HER2,
  TNBC=as.character(as.integer(sampleAnno$TNBC)),
  PAM50=sampleAnno$pam50.full,stringsAsFactors=FALSE
  )
rownames(sample_anno)<-colnames(temp1)
sample_anno<-sample_anno[,ncol(sample_anno):1]

r_anno<-data.frame(feature=as.character(sub(" .+","",annoObj[rownames(temp1),"featureClass"])),
  island=as.character(annoObj[rownames(temp1),"cgiClass"]),
  Atac=as.character(annoObj[rownames(temp1),"hasAtacOverlap"]),
  ESR1=as.character(tfMat[rownames(temp1),"ESR1"]),
  FOXA1=as.character(tfMat[rownames(temp1),"FOXA1"]),
  GATA3=as.character(tfMat[rownames(temp1),"GATA3"]),
  EZH2=as.character(tfMat[rownames(temp1),"EZH2"]),
  SUZ12=as.character(tfMat[rownames(temp1),"SUZ12"]),stringsAsFactors=FALSE
  )
rownames(r_anno)<-rownames(temp1)
r_anno<-r_anno[,ncol(r_anno):1]

my_colour = list(#unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    # ER = c("[Not Available]"="#FFFF33",
    #   "[Not Evaluated]"="#FF7F00",
    #   "Equivocal"="#377EB8",
    #   "Indeterminate"="#984EA3",
    #   "Negative"="#E41A1C",
    #   "Positive"="#4DAF4A"
    #   ),
    # PR = c("[Not Available]"="#FFFF33",
    #   "[Not Evaluated]"="#FF7F00",
    #   "Equivocal"="#377EB8",
    #   "Indeterminate"="#984EA3",
    #   "Negative"="#E41A1C",
    #   "Positive"="#4DAF4A"
    #   ),
    # HER2 = c("[Not Available]"="#FFFF33",
    #   "[Not Evaluated]"="#FF7F00",
    #   "Equivocal"="#377EB8",
    #   "Indeterminate"="#984EA3",
    #   "Negative"="#E41A1C",
    #   "Positive"="#4DAF4A"
    #   ),
    TNBC = c("1"="black","0"="lightgrey"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080"),
    feature=c("distal"=brewer.pal(9,"Blues")[c(5)],"promoter"=brewer.pal(9,"Reds")[c(5)],"proximal"="white"),
    island=c("ocean"=brewer.pal(9,"Blues")[c(5)],"shore"="white","cgi"=brewer.pal(9,"Reds")[c(5)]),
    Atac=c("0"="white","1"=brewer.pal(9,"Greens")[c(5)]),
    ESR1=c("0"="white","1"="orange"),
    FOXA1=c("0"="white","1"="orange"),
    GATA3=c("0"="white","1"="orange"),
    EZH2=c("0"="white","1"="purple"),
    SUZ12=c("0"="white","1"="purple")
    ,stringsAsFactors=FALSE
  )

tiff(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_adjBeta_withRowAnno.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
#pdf(paste0(HOME,"/final_fig1.pdf"),width=12,height=12,useDingbats=FALSE)
#tiff(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_adjBeta_withRowAnno.tiff"),width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="",cutree_cols=5,fontsize = 12,cellwidth = .8,# cellheight =.1 
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_row=r_anno
)
dev.off()

tiff(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_unadjBeta_withRowAnno.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
pheatmap(temp3,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F,treeheight_row =0,treeheight_col =0
  ,main="",cutree_cols=5,fontsize = 12,cellwidth = .8
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_row=r_anno
)
dev.off()

##first panel
a1<-image_read(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_unadjBeta_withRowAnno.tiff"))
a1<-image_crop(a1,"500x4500+690")
a1<-image_scale(a1,"x2000")
#image_write(a2, path = paste0(HOME,"/plos_figure5_temp2.tiff"), format = "tiff")

a2<-image_read(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_adjBeta_withRowAnno.tiff"))
a2<-image_crop(a2,"500x4500+3360")
a2<-image_scale(a2,"x2000")
#image_write(a2, path = paste0(HOME,"/plos_figure5_temp.tiff"), format = "tiff")

a3<-image_read(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_adjBeta_withRowAnno.tiff"))
a3<-image_crop(a3,"2180x4500+1190")
a3<-image_annotate(a3, "Adjusted", size = 76, gravity="North", color = "black",font="sans",weight=700)
a3<-image_scale(a3,"x2000")
#image_write(a3, path = paste0(HOME,"/plos_figure5_temp.tiff"), format = "tiff")

a4<-image_read(paste0(HOME,"/plos_figure5_top5k_heatmap_pear_eucl_adjClust_unadjBeta_withRowAnno.tiff"))
a4<-image_crop(a4,"2180x4500+1190")
a4<-image_annotate(a4, "Unadjusted", size = 76, gravity="North", color = "black",font="sans",weight=700)
a4<-image_scale(a4,"x2000")
#image_write(a4, path = paste0(HOME,"/plos_figure5_temp.tiff"), format = "tiff")

out<-image_append(c(a1,a3,image_blank(width=25, height=2000, color = "white"),
  a4,a2
  ),stack = F)
out<-image_scale(out,"2250x")

image_write(out, path = paste0(HOME,"/plos_figure5_panel_a.tiff"), format = "tiff",depth=8)

###Calculate beta difference pre/post correction for top5k by 5-cluster split and PAM50
  #1. islands with ezh2+suz12
  #2. ATAC with FOXA1+GATA3
  #3. distal non-ATAC/TF

##1.
sel<-annoObj[rownames(temp1),"cgiClass"]=="cgi" & tfMat[rownames(temp1),"EZH2"]==1 & tfMat[rownames(temp1),"SUZ12"]==1

table(sel)
# sel
# FALSE  TRUE 
#  4119   881 

##
tiff(paste0(HOME,"/plos_figure5_panel_b.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0.5,1),mar=c(6.1,5.1,3.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=.8,ps=30,new=F)
# tiff(paste0(HOME,"/plos_figure5_panel_b.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# par(mar=c(6.1, 4.1, 3.1 ,2.1),fig=c(0,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2)
plot(1,cex.main=1.
  ,type="n",xlab=""
  ,main="CGI with ENCODE EZH2+SUZ12 (N=881)",ylab="Beta"
  ,las=1,axes=F,xlim=c(0,35),ylim=c(0,1)
  )
abline(h=c(0,.25,.5,.75,1),lty=c(1,2,1,2,1),lwd=3,col="lightgrey")
abline(v=c(16,30),lty=c(2),lwd=3,col=1)
##loop meth
for(j in 1:5) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      pal<-c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j]
      tmp<-temp1
      kk<-paste0("adj ",j)
    } else {
      pal<-paste0(c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j],"50")
      tmp<-temp3
      kk<-paste0("unadj ",j)
    }
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,c1==j])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,c1==j]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,c1==j]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,c1==j]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,c1==j]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-1)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)

  }
}
##loop pam50
for(j in 19:16) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)]
      tmp<-temp1
      kk<-paste0("adj ",i)
    } else {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-paste0(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",i)
    }
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,sampleAnno$pam50.full==i])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-15)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)
  }
}
##loop tnbc
for(j in 30:31) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-abs(30-j)
      pal<-c("0"="#000000","1"="#E41A1C")[abs(29-j)]
      tmp<-temp1
      kk<-paste0("adj ",c("Lum","TNBC")[i+1])
    } else {
      i<-abs(30-j)
      pal<-paste0(c("0"="#000000","1"="#E41A1C")[abs(29-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",c("Lum","TNBC")[i+1])
    }
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-29)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)
  }
}
axis(2,lwd=2,las=1,at=seq(0,1,by=.5),font=2)
#legend("topleft",legend=c("adjusted beta","unadjusted beta"),bty="n",col=c("#E41A1C","#377EB8"),lwd=3)

##plus densplot

##loop meth
for(j in 1:5) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      pal<-c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j]
      tmp<-temp1
      kk<-paste0("adj ",j)
    } else {
      pal<-paste0(c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j],"50")
      tmp<-temp3
      kk<-paste0("unadj ",j)
    }
  if(j==1) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(0,.2,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==2) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.2,.4,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==3) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(0,.2,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==4) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.2,.4,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==5) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(0,.2,.03,.18),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,c1==j])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=j,col=pal,bty="n")
  }
}

##loop pam50
for(j in 19:16) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)]
      tmp<-temp1
      kk<-paste0("adj ",i)
    } else {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-paste0(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",i)
    }
  if(i=="Basal") par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.4,.6,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="LumB") par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.6,.8,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="LumA") par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.4,.6,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="Her2") par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.6,.8,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,sampleAnno$pam50.full==i])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=i,col=pal,bty="n")
  }
}
##loop tnbc
for(j in 30:31) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-abs(30-j)
      pal<-c("0"="#000000","1"="#E41A1C")[abs(29-j)]
      tmp<-temp1
      kk<-paste0("adj ",c("Lum","TNBC")[i+1])
    } else {
      i<-abs(30-j)
      pal<-paste0(c("0"="#000000","1"="#E41A1C")[abs(29-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",c("Lum","TNBC")[i+1])
    }
  if(i==0) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.8,1,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i==1) par(mar=c(0.5, 2.5, 0.5 ,0.5),fig=c(.8,1,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=c("Lum","TNBC")[i+1],col=pal,bty="n")
  }
}

dev.off()

##2.
sel<-tfMat[rownames(temp1),"FOXA1"]==1 & tfMat[rownames(temp1),"GATA3"]==1 & annoObj[rownames(temp1),"hasAtacOverlap"]==1

table(sel)
# sel
# FALSE  TRUE 
#  4948    52 

##
tiff(paste0(HOME,"/plos_figure5_panel_c.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0.5,1),mar=c(6.1,5.1,3.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=.8,ps=30,new=F)
# tiff(paste0(HOME,"/fig5_panel_foxaGata.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# par(mar=c(6.1, 4.1, 3.1 ,2.1),fig=c(0,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2)
plot(1,cex.main=1.
  ,type="n",xlab=""
  ,main="TCGA ATAC with ENCODE FOXA1+GATA3 (N=52)",ylab="Beta"
  ,las=1,axes=F,xlim=c(0,35),ylim=c(0,1)
  )
abline(h=c(0,.25,.5,.75,1),lty=c(1,2,1,2,1),lwd=3,col="lightgrey")
abline(v=c(16,30),lty=c(2),lwd=3,col=1)
##loop meth
for(j in 1:5) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      pal<-c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j]
      tmp<-temp1
      kk<-paste0("adj ",j)
    } else {
      pal<-paste0(c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j],"50")
      tmp<-temp3
      kk<-paste0("unadj ",j)
    }
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,c1==j])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,c1==j]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,c1==j]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,c1==j]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,c1==j]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-1)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)

  }
}
##loop pam50
for(j in 19:16) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)]
      tmp<-temp1
      kk<-paste0("adj ",i)
    } else {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-paste0(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",i)
    }
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,sampleAnno$pam50.full==i])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-15)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)
  }
}
##loop tnbc
for(j in 30:31) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-abs(30-j)
      pal<-c("0"="#000000","1"="#E41A1C")[abs(29-j)]
      tmp<-temp1
      kk<-paste0("adj ",c("Lum","TNBC")[i+1])
    } else {
      i<-abs(30-j)
      pal<-paste0(c("0"="#000000","1"="#E41A1C")[abs(29-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",c("Lum","TNBC")[i+1])
    }
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-29)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)
  }
}
axis(2,lwd=2,las=1,at=seq(0,1,by=.5),font=2)
#legend("topleft",legend=c("adjusted beta","unadjusted beta"),bty="n",col=c("#E41A1C","#377EB8"),lwd=3)

##plus densplot

##loop meth
for(j in 1:5) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      pal<-c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j]
      tmp<-temp1
      kk<-paste0("adj ",j)
    } else {
      pal<-paste0(c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j],"50")
      tmp<-temp3
      kk<-paste0("unadj ",j)
    }
  if(j==1) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(0,.2,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==2) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.2,.4,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==3) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(0,.2,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==4) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.2,.4,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==5) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(0,.2,.03,.18),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,c1==j])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=j,col=pal,bty="n")
  }
}

##loop pam50
for(j in 19:16) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)]
      tmp<-temp1
      kk<-paste0("adj ",i)
    } else {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-paste0(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",i)
    }
  if(i=="Basal") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.4,.6,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="LumB") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.6,.8,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="LumA") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.4,.6,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="Her2") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.6,.8,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,sampleAnno$pam50.full==i])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=i,col=pal,bty="n")
  }
}
##loop tnbc
for(j in 30:31) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-abs(30-j)
      pal<-c("0"="#000000","1"="#E41A1C")[abs(29-j)]
      tmp<-temp1
      kk<-paste0("adj ",c("Lum","TNBC")[i+1])
    } else {
      i<-abs(30-j)
      pal<-paste0(c("0"="#000000","1"="#E41A1C")[abs(29-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",c("Lum","TNBC")[i+1])
    }
  if(i==0) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.8,1,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i==1) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.8,1,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=c("Lum","TNBC")[i+1],col=pal,bty="n")
  }
}

dev.off()

##3. 
sel<-sub(" .+","",annoObj[rownames(temp1),"featureClass"])=="distal" & rowSums(tfMat[rownames(temp1),])==0

table(sel)
# sel
# FALSE  TRUE 
#  4524   476 

##
tiff(paste0(HOME,"/plos_sFigure3.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0.5,1),mar=c(6.1,5.1,3.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=.8,ps=30,new=F)
#tiff(paste0(HOME,"/fig5_panel_distalNoTf.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
#par(mar=c(6.1, 4.1, 3.1 ,2.1),fig=c(0,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2)
plot(1,cex.main=1.
  ,type="n",xlab=""
  ,main="Distal without ENCODE TFBS (N=476)",ylab="Beta"
  ,las=1,axes=F,xlim=c(0,35),ylim=c(0,1)
  )
abline(h=c(0,.25,.5,.75,1),lty=c(1,2,1,2,1),lwd=3,col="lightgrey")
abline(v=c(16,30),lty=c(2),lwd=3,col=1)
##loop meth
for(j in 1:5) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      pal<-c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j]
      tmp<-temp1
      kk<-paste0("adj ",j)
    } else {
      pal<-paste0(c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j],"50")
      tmp<-temp3
      kk<-paste0("unadj ",j)
    }
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,c1==j])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,c1==j]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,c1==j]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,c1==j]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-1)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,c1==j]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-1)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)

  }
}
##loop pam50
for(j in 19:16) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)]
      tmp<-temp1
      kk<-paste0("adj ",i)
    } else {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-paste0(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",i)
    }
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,sampleAnno$pam50.full==i])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-15)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,sampleAnno$pam50.full==i]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-15)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)
  }
}
##loop tnbc
for(j in 30:31) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-abs(30-j)
      pal<-c("0"="#000000","1"="#E41A1C")[abs(29-j)]
      tmp<-temp1
      kk<-paste0("adj ",c("Lum","TNBC")[i+1])
    } else {
      i<-abs(30-j)
      pal<-paste0(c("0"="#000000","1"="#E41A1C")[abs(29-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",c("Lum","TNBC")[i+1])
    }
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(median(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i])),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),.25),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,.25),rep(quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),.75),2),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(.25,.25),y=quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),c(.25,.75)),lwd=4,col=pal)
    lines(x=j+(2*(j-29)+1*(k-1))+c(-.25,-.25),y=quantile(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i]),c(.25,.75)),lwd=4,col=pal)
    axis(1,at=j+(2*(j-29)+1*(k-1)),labels=kk,las=2,lwd=0,font=2,line=-1)
  }
}
axis(2,lwd=2,las=1,at=seq(0,1,by=.5),font=2)
#legend("topleft",legend=c("adjusted beta","unadjusted beta"),bty="n",col=c("#E41A1C","#377EB8"),lwd=3)

##plus densplot

##loop meth
for(j in 1:5) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      pal<-c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j]
      tmp<-temp1
      kk<-paste0("adj ",j)
    } else {
      pal<-paste0(c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[j],"50")
      tmp<-temp3
      kk<-paste0("unadj ",j)
    }
  if(j==1) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(0,.2,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==2) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.2,.4,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==3) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(0,.2,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==4) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.2,.4,.18,.32),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(j==5) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(0,.2,.03,.18),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,c1==j])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=j,col=pal,bty="n")
  }
}

##loop pam50
for(j in 19:16) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)]
      tmp<-temp1
      kk<-paste0("adj ",i)
    } else {
      i<-names(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3"))[abs(15-j)]
      pal<-paste0(c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3")[abs(15-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",i)
    }
  if(i=="Basal") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.4,.6,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="LumB") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.6,.8,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="LumA") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.4,.6,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i=="Her2") par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.6,.8,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,sampleAnno$pam50.full==i])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=i,col=pal,bty="n")
  }
}
##loop tnbc
for(j in 30:31) {
  for(k in 1:2) {
    #orig
    if( k %%2==0) {
      i<-abs(30-j)
      pal<-c("0"="#000000","1"="#E41A1C")[abs(29-j)]
      tmp<-temp1
      kk<-paste0("adj ",c("Lum","TNBC")[i+1])
    } else {
      i<-abs(30-j)
      pal<-paste0(c("0"="#000000","1"="#E41A1C")[abs(29-j)],"50")
      tmp<-temp3
      kk<-paste0("unadj ",c("Lum","TNBC")[i+1])
    }
  if(i==0) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.8,1,.33,.48),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  if(i==1) par(mar=c(0.5, 2.5, 0.5 ,.5),fig=c(.8,1,.18,.33),font=2,font.axis=2,font.lab=2,font.sub=2,new=TRUE)
  plot(1,
    ,type="n",xlab="Beta"
    ,main="",ylab="Density"
    ,las=1,axes=F,xlim=c(-.1,1.1),ylim=c(0,7.5)
    )
  abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
  axis(1,lwd=0,las=1,at=seq(0,1,by=.5),font=2,cex.axis=.75,line=-.8)
  axis(2,at=seq(0,5,by=2.5),lwd=0,las=1,font=2,cex.axis=.75,line=-.8)
  lines(x=density(as.vector(tmp[sel,as.integer(sampleAnno$TNBC)==i])),lwd=4,col=pal)
  text(x=.5,y=4.5,labels=c("Lum","TNBC")[i+1],col=pal,bty="n")
  }
}

dev.off()

##first panel
a1<-image_read(paste0(HOME,"/plos_figure5_panel_a.tiff"))

a2<-image_read(paste0(HOME,"/plos_figure5_panel_b.tiff"))
a2<-image_annotate(a2, "B", size = 152, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)
a2<-image_scale(a2,"1125x")
image_write(a2,paste0(HOME,"/plos_figure5_panel_b_final.tiff"))

a3<-image_read(paste0(HOME,"/plos_figure5_panel_c.tiff"))
a3<-image_annotate(a3, "C", size = 152, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)
a3<-image_scale(a3,"1125x")
image_write(a3,paste0(HOME,"/plos_figure5_panel_c_final.tiff"))

a4<-image_read(paste0(HOME,"/plos_sFigure3.tiff"))
a4<-image_scale(a4,"1125x")
image_write(a4,paste0(HOME,"/plos_sFigure3.tiff"))

a2<-image_read(paste0(HOME,"/plos_figure5_panel_b_final.tiff"))
a3<-image_read(paste0(HOME,"/plos_figure5_panel_c_final.tiff"))

out<-image_append(c(a2,image_blank(width=50, height=1125, color = "white"),a3
  ),stack = F)
out<-image_scale(out,"2250x")

out<-image_append(c(a1,out
  ),stack = T)
out<-image_scale(out,"x2625")

# image_write(image_blank(width=131, height=2625, color = "none"),density=300,path = paste0(HOME,"/plos_figure5_tmp.tiff"),format = "tiff")
# a5<-image_read(paste0(HOME,"/plos_figure5_tmp.tiff"),density=300)

# out<-image_append(c(a5,
#   out,
#   a5
#   ),stack = F)
# out<-image_scale(out,"2250x")

image_write(out, path = paste0(HOME,"/plos_figure5_final.tiff"), format = "tiff",depth=8)

################################################################################
################################################################################

##Same source as original probe mapping info but now migrated by creator [wanding.zhou@pennmedicine.upenn.edu]
#http://zwdzwd.github.io/InfiniumAnnotation

tmp<-tempfile()
download.file("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM450/HM450.hg19.manifest.tsv.gz",tmp)

p.info<-read.table(tmp,sep="\t",header=T,as.is=T)
rownames(p.info)<-p.info$probeID

all(rownames(annoObj) %in% p.info$probeID)
#[1] TRUE

length(intersect(p.info$probeID , rownames(annoObj)))
#[1] 421368

p.info<-p.info[intersect(rownames(annoObj) , p.info$probeID),]

dim(p.info)
#[1] 421368     57

all(rownames(annoObj) == p.info$probeID)
#[1] TRUE

unlink(tmp)

################################################################################
################################################################################
##Create CpG O/E based analysis

# ##O/E by type
# tiff(paste0(HOME,"/fig6_panel_a1.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# par(mar=c(5.1, 5.6, 4.1 ,.6),fig=c(0,1,0,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=FALSE)
# plot(1,type="n",main="",sub="",xlab="",ylab="",
#   axes=F,xlim=c(0,1.6),ylim=c(0,1)
#   )
# abline(h=c(0),lty=c(1),lwd=3,col="lightgrey")
# axis(1,lwd=2,las=1,at=seq(0,1.5,by=.5),font=2,cex.axis=2,line=-0.25)
# axis(2,at=seq(0,1,by=.5),lwd=2,las=1,font=2,cex.axis=2,line=0.25)
# legend("topright",legend=c("I","II"),col=brewer.pal(9,"Oranges")[c(7,5)],bty="n",pch=15,cex=2)
# mtext("CpG O/E",side=1,at=.75,line=2.5,font=2,cex=2)
# mtext("Scaled density",side=2,line=4,font=2,cex=2)

# dd<-density(annoObj[,"saxonovOE"][p.info$designType=="I"])
# dd$y<-dd$y / max(dd$y)
# lines(dd,col=brewer.pal(9,"Oranges")[c(7)],lwd=3,lty=1)

# dd<-density(annoObj[,"saxonovOE"][p.info$designType=="II"])
# dd$y<-dd$y / max(dd$y)
# lines(dd,col=brewer.pal(9,"Oranges")[c(5)],lwd=3,lty=1)
# dev.off()

# ##type by bin
# tiff(paste0(HOME,"/fig6_panel_a2.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# par(mar=c(5.1, 5.6, 4.1 ,.6),fig=c(0,1,0,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=FALSE)
# cgBins<-cut(annoObj[,"saxonovOE"],breaks=quantile(annoObj[,"saxonovOE"],seq(0,1,length.out=101)),labels=F,include.lowest=TRUE)
# barplot( t( table(cgBins,p.info$designType) / rowSums(table(cgBins,p.info$designType)) ) ,col=brewer.pal(9,"Oranges")[c(7,5)],axes=F,axisnames=FALSE,xlim=c(0,140))
# legend("topright",legend=c("I","II"),col=brewer.pal(9,"Oranges")[c(7,5)],bty="n",pch=15,cex=2)
# axis(1,lwd=2,las=1,at=seq(0,120,by=60),labels=c("1","50","100"),font=2,cex.axis=2,line=-0.25)
# axis(2,at=seq(0,1,by=.5),lwd=2,las=1,font=2,cex.axis=2,line=0.25)
# mtext("100 CpG O/E bins",side=1,at=60,line=2.5,font=2,cex=2)
# mtext("Fraction",side=2,line=4,font=2,cex=2)
# dev.off()

# ##meth by bin+type
# tiff(paste0(HOME,"/fig6_panel_a3.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
# cgBins<-cut(annoObj[,"saxonovOE"],breaks=quantile(annoObj[,"saxonovOE"],seq(0,1,length.out=101)),labels=F,include.lowest=TRUE)
# par(mar=c(5.1, 6.1, 4.1 ,6.1),fig=c(0,1,0,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=FALSE)
# plot(1,type="n",main="",sub="",xlab="",ylab="",
#   axes=F,xlim=c(0,100),ylim=c(0,1)
#   )
# ##I/II density
# dd<-density(cgBins[p.info$designType=="I"])
# dd$y<-dd$y / max(dd$y)
# lines(dd,col="grey",lwd=3,lty=1)
# dd<-density(cgBins[p.info$designType=="II"])
# dd$y<-dd$y / max(dd$y)
# lines(dd,col="grey",lwd=3,lty=2)
# ##legends
# axis(1,lwd=2,las=1,at=seq(0,100,by=50),font=2,cex.axis=2,line=-0.25)
# axis(2,at=seq(0,1,by=.5),lwd=2,las=1,font=2,cex.axis=2,line=0.25)
# axis(4,at=seq(0,1,by=.5),lwd=2,las=1,font=2,cex.axis=2,line=0.25,col="grey",col.axis="grey")
# mtext("100 CpG O/E bins",side=1,at=50,line=2.5,font=2,cex=2)
# mtext("Median beta",side=2,line=4,font=2,cex=2)
# mtext("Scaled density",side=4,line=4,font=2,cex=2,col="grey")
# legend("topright",legend=paste0(rep(c("adj","unadj"),each=2),c(" I"," II")),col=rep(brewer.pal(9,"Oranges")[c(7,5)],2),bty="n",lty=c(1,1,2,2),lwd=3,cex=2)
# legend(x=65,y=.45,legend=c(" I"," II"),col="grey",bty="n",lty=c(1,2),lwd=3,cex=2,text.col="grey")
# #adj
# lines(unlist(lapply(split(as.vector(betaAdj[p.info$designType=="I",]),rep(cgBins[p.info$designType=="I"],ncol(betaAdj))),median)),col=brewer.pal(9,"Oranges")[7],lty=1,lwd=3)
# lines(unlist(lapply(split(as.vector(betaAdj[p.info$designType=="II",]),rep(cgBins[p.info$designType=="II"],ncol(betaAdj))),median)),col=brewer.pal(9,"Oranges")[5],lty=1,lwd=3)
# #orig
# lines(unlist(lapply(split(as.vector(betaOrig[p.info$designType=="I",]),rep(cgBins[p.info$designType=="I"],ncol(betaAdj))),median)),col=brewer.pal(9,"Oranges")[7],lty=2,lwd=3)
# lines(unlist(lapply(split(as.vector(betaOrig[p.info$designType=="II",]),rep(cgBins[p.info$designType=="II"],ncol(betaAdj))),median)),col=brewer.pal(9,"Oranges")[5],lty=2,lwd=3)

# dev.off()

# ##first panel
# a1<-image_read(paste0(HOME,"/fig6_panel_a1.tiff"))

# a2<-image_read(paste0(HOME,"/fig6_panel_a2.tiff"))

# a3<-image_read(paste0(HOME,"/fig6_panel_a3.tiff"))

# out<-image_append(c(a1,a2,a3
#   ),stack = F)
# out<-image_scale(out,"4150x")

# image_write(out, path = paste0(HOME,"/fig6_panel_a.tiff"), format = "tiff")

################################################################################
################################################################################
##Do analysis of stability

##do multicore
no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

clusterEvalQ(cl, {
  library("flexmix")
})

##check avg delta by perturbation
unlist(lapply(seq(0.01,.19,by=.02),function(x) {
  mean(unlist(lapply(1:1000,function(y) mean(abs(fracTum-(fracTum-rnorm(length(fracTum),mean=0,sd=x)))) )) )
}))
 # [1] 0.00797237 0.02390659 0.03991191 0.05573767 0.07181002 0.08782423
 # [7] 0.10383002 0.11980487 0.13541895 0.15180635
##check proportion outside 0-1
unlist(lapply(seq(0.01,.19,by=.02),function(x) {
  mean(unlist(lapply(1:1000,function(y) { z<-fracTum-rnorm(length(fracTum),mean=0,sd=x) ; sum(z>1 | z<0)/length(z)   })) ) 
}))
 # [1] 0.002423810 0.003103175 0.005709524 0.010680952 0.017671429 0.026779365
 # [7] 0.037568254 0.049919048 0.064147619 0.079257143

##number of iterations
N<-10

###Change purity but run same set
str(testDat2)
 # num [1:5000, 1:630] 0.341 0.057 0.012 0.022 0 0.004 0.021 0.892 0.022 0.163 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg09248054" "cg25340711" "cg06443533" "cg16601494" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

#tiff(paste0(HOME,"/fig7_panel_a.tiff"),width=8*500,height=8*500,units="px",res=500,compression="lzw")
#par(mar=c(5.1, 6.1, 4.1 ,6.1),fig=c(0,1,0,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=FALSE)
tiff(paste0(HOME,"/plos_figure6_panel_a.tiff"),width=7.5*300,height=7.5*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)

plot(1,type="n",main="",sub="",xlab="SD of noise",ylab="Similarity",
  axes=F,xlim=c(0,.2),ylim=c(0,1)
  )
abline(v=seq(0.02,.18,by=.02),col="lightgray",lwd=2,lty=2)
legend("bottomleft",legend=c("FlexMix group call concordance","Adjusted vs perturbed adjustment correlation","Normal vs inferred normal correlation"),
  bty="n",col=c("black","red","orange"),pch=16)
axis(2,las=1)
axis(1,las=3,at=seq(0.01,.19,by=.02))

sink(file=paste0(HOME,"/perturbationTestInTop5000_groupConcordanceStats.txt"))
###RUN
set.seed(20211208)
for (i in 1:N) {
  cat("Run",i," of ",N,"\n")
  ##betaData filtered for chrX/Y
  betaRun<-cbind(seed=1:nrow(testDat2),testDat2)
  betaNames<-colnames(testDat2)
  ##add fudgefactor
  fracTum2<-fracTum + rnorm(length(fracTum),mean=0,sd=seq(0.01,.19,by=.02)[i])
  fracTum2[fracTum2>1]<-1
  fracTum2[fracTum2<0]<-0
  b<-parRapply(cl = cl, betaRun, adjustBeta,purity=fracTum2,snames=betaNames,seed=TRUE)
  ##gather adjusted matrix
  temp1<-do.call("rbind",lapply(b,function(x) x$y.tum))
  temp2<-do.call("rbind",lapply(b,function(x) x$y.orig))
  temp3<-do.call("rbind",lapply(b,function(x) x$y.norm))
  temp4<-beta_norm[rownames(testDat2),]
  temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))
  temp5<-betaAdj[rownames(testDat2),]

  grFilt<-unlist(lapply(b,function(x) max(x$groups))) == unlist(lapply(res,function(x) max(x$groups)))
  cat("Fraction of class calls with same N",sum(grFilt)/length(grFilt),"\n\n")
  temp6<-do.call("rbind",lapply(b[grFilt],function(x) x$groups))
  temp7<-do.call("rbind",lapply(res[grFilt],function(x) x$groups))
  
  xx<-unlist(lapply(1:sum(grFilt),function(x) sum(apply(table(temp6[x,],temp7[x,]),2,max))/ncol(temp6) ))
  points(seq(0.01,.19,by=.02)[i]-.005,median(xx),col=1,pch=16)
  lines(rep(seq(0.01,.19,by=.02)[i],2)-.005,quantile(xx,c(.025,.975)),col=1,lwd=2)

  points(seq(0.01,.19,by=.02)[i],median(diag(cor(temp1,temp5))),pch=16,col=2)
  lines(rep(seq(0.01,.19,by=.02)[i],2),quantile(diag(cor(temp1,temp5)),c(.025,.975)),lwd=2,col=2)

  points(seq(0.01,.19,by=.02)[i]+.005,median(cor(temp3,rowMeans(temp4))),pch=16,col="orange")
  lines(rep(seq(0.01,.19,by=.02)[i],2)+.005,quantile(cor(temp3,rowMeans(temp4)),c(.025,.975)),lwd=2,col="orange")

  #points(seq(0.01,.19,by=.02)[i],median(diag(cor(t(temp1),t(temp5)))),pch=16,col=2)
  #lines(rep(seq(0.01,.19,by=.02)[i],2),quantile(diag(cor(t(temp1),t(temp5))),c(.025,.975)),lwd=2,col=2)

  save(b,file=paste0(HOME,"/perturbationTestInTop5000_iteration",i,"_runResults.RData"))
  cat("#######","\n\n")
}
sink()
dev.off()
rm(b,i,N,betaRun,betaNames)

fig<-image_scale(image_read(paste0(HOME,"/plos_figure6_panel_a.tiff")),"1125x")
#fig<-image_annotate(fig, "A", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)
image_write(fig,paste0(HOME,"/plos_figure6_panel_a.tiff"))
image_destroy(fig)
rm(fig)


################################################################################
##Show correlation individual tumor/normal to mean normal

temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.orig))
temp3<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp4<-beta_norm_adj[rownames(temp1),]
temp4<-t(apply(temp4,1,function(x) { if(any(is.na(x))) { x[is.na(x)]<-median(x[!is.na(x)]) } ; return(x) } ))
temp5<-betaAdj[rownames(temp1),]

#pdf(paste0(HOME,"/sFig1_final.pdf"),width=8,height=8,useDingbats=FALSE)
#par(mar=c(5.1, 6.1, 4.1 ,6.1),fig=c(0,1,0,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=FALSE)
tiff(paste0(HOME,"/plos_sFig1.tiff"),width=7.5*300,height=7.5*300,units="px",res=300,compression="lzw")
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)

plot(1,type="n",main="",sub="",xlab="Sample correlation to normal breast methylation (N=630)",ylab="Scaled density",
  axes=F,xlim=c(-1,1),ylim=c(0,1)
  )
axis(2,las=1)
axis(1,las=1,at=c(-1,0,1))

rescorrs<-NULL
for(i in 1:ncol(temp1)) {
  rescorrs<-c(rescorrs,cor(temp2[,i],rowMeans(temp4)))
}
xx<-density(rescorrs)
xx$y<-xx$y/max(xx$y)
lines(xx,col="black",lwd=3)

rescorrs2<-NULL
for(i in 1:ncol(temp1)) {
  rescorrs2<-c(rescorrs2,cor(temp1[,i],rowMeans(temp4)))
}
xx<-density(rescorrs2)
xx$y<-xx$y/max(xx$y)
lines(xx,col="red",lwd=3)

rescorrs3<-NULL
for(i in 1:ncol(temp1)) {
   rescorrs3<-c(rescorrs3,cor(temp3[,i],rowMeans(temp4)))
}
xx<-density(rescorrs3)
xx$y<-xx$y/max(xx$y)
lines(xx,col="orange",lwd=3)

legend("topleft",legend=paste(c("Unadjusted tumor (median r=","Adjusted tumor (median r=","Inferred normal (median r="),
  c(signif(median(rescorrs),2),
    signif(median(rescorrs2),2),
    signif(median(rescorrs3),2)),")",sep=""),
  bty="n",col=c("black","red","orange"),pch=16,cex=.75)

abline(h=0,col="lightgray",lwd=2,lty=1)
dev.off()

fig<-image_scale(image_read(paste0(HOME,"/plos_sFig1.tiff")),"1125x")
#fig<-image_annotate(fig, "A", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)
image_write(fig,paste0(HOME,"/plos_sFig1.tiff"))
image_destroy(fig)
rm(fig)

################################################################################
################################################################################
##Use MethylResolver to estimate pre/post purity

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
##Show additional examples of beta correciton for fig1

################################################################################
##Show wether CN-status influences presence of 3 groups

################################################################################
##Show how IM-group dissolves in corrected data


save.image(paste0(HOME,"/workspace_tcgaBrca_full.RData"))










###END