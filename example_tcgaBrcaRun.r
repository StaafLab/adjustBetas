#####======================================================================#####
### Correct TCGA BRCA beta values for infiltration
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

##Git for data download and preprocessing of TCGA BRCA data (Staaf-lab pipeline)
#https://github.com/StaafLab/tcgaBrca
#resluting data object can be obtained on demand from MA/Staaf-lab to avoid extensive download/run times

##Git for beta adjustment
#https://github.com/StaafLab/adjustBetas

################################################################################
##Set directory paths

##set/create own home directory below:
HOME<-"~/hdd1/tcgaBrca"
MANIFEST<-"~/Documents/tcgaBrca/manifest"

GIT<-"~/Documents/adjustBetas"

##tumor type
TUMOR_TYPE<-"brca"

################################################################################
##load required packages

if(!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel") }

library(doParallel)

if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel") }

library(parallel)

##source - flexmix loaded on source
source(paste0(GIT,"/function_correctBetas.r"))

################################################################################
##Get data object and purity vector

load(file=paste0(HOME,"/","finalWorkspace_atacCnGexMeWes_withAnnotations.RData"))

ls()
#  [1] "adjustBeta" "annoObj"    "betaOrig"   "DATA"       "dataAtac"   "dataCn"     "dataMut"    "dataSeg"   
#  [9] "gexCounts"  "gexFpkm"    "gexUq"      "GIT"        "HOME"       "MANIFEST"   "sampleAnno" "TUMOR_TYPE"

fracTum<-read.table(file=paste0(DATA,"/me/sampleTumorDnaFractionVector.txt"),header=TRUE,as.is=T,sep="\t")

str(fracTum)
# 'data.frame':   630 obs. of  2 variables:
#  $ sample            : chr  "TCGA-3C-AAAU-01" "TCGA-3C-AALI-01" "TCGA-3C-AALJ-01" "TCGA-3C-AALK-01" ...
#  $ tumor_dna_fraction: num  0.85 0.73 0.69 0.62 0.64 0.53 0.89 0.47 0.56 0.65 ...

rownames(fracTum)<-fracTum$sample

identical(rownames(fracTum),sub(".$","",colnames(betaOrig)))
#[1] TRUE

rownames(fracTum)<-colnames(betaOrig)

################################################################################
###Beta correction on top 5000 most varying probes

##filter X/Y and get top5k
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

##create data set and run function
testDat2<-betaData[varF,]

str(testDat2)
 # num [1:5000, 1:630] 0.341 0.057 0.012 0.022 0 0.004 0.021 0.892 0.022 0.163 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg09248054" "cg25340711" "cg06443533" "cg16601494" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

################################################################################
###correct betas using multicore 

#need to feed seed to function so fully reproducible

no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

##estimated runtime is ~500 sec on 7 cores

clusterEvalQ(cl, {
  library("flexmix")
})

##add rng seed to each row and pass to each paralell instance
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

table( unlist(lapply(1:length(res),function(x) { all( unlist(res[[x]],use.names=FALSE) == unlist(res[[x]],use.names=FALSE) ) }) ) )
# TRUE 
# 5000 
table( unlist(lapply(1:length(res),function(x) { all( unlist(res[[x]],use.names=FALSE) == unlist(res2[[x]],use.names=FALSE) ) }) ) )
# TRUE 
# 5000 

##same results across runs.
rm(betaRun,betaNames,testDat,testDat2,res2,cl,no_cores,varF)

################################################################################
## object "res" is a list containing collected parameters of fitted lines as well as estimates of TME and pure tumor methylation

##One CpG output list entry contains following:
str(res[1])
# List of 1
#  $ cg09248054:List of 10
#   ..$ y.norm          : Named num [1:630] 0.251 0.104 0.077 0.512 0.236 0.433 0.358 0.774 0.326 0.293 ...
#   .. ..- attr(*, "names")= chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#   ..$ y.tum           : Named num [1:630] 0.357 0.045 0.019 0.618 0.919 0.539 1 0.88 1 0.976 ...
#   .. ..- attr(*, "names")= chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#   ..$ y.orig          : Named num [1:630] 0.341 0.061 0.037 0.578 0.673 0.489 0.966 0.824 0.709 0.737 ...
#   .. ..- attr(*, "names")= chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#   ..$ groups          : int [1:630] 3 1 1 3 2 3 2 3 2 2 ...
#   ..$ n.groups        : int 3
#   ..$ med.norm        : num 0.318
#   ..$ glob.cor        : num 0.246
#   ..$ avg.betaDiff    : num -0.14
#   ..$ model.intercepts: num [1:3] 0.065 1.011 0.542
#   ..$ model.slopes    : num [1:3] 0.058 -0.683 -0.106

##Simple way to extract matrices from object:

#adjusted tumor betas:
betaAdj<-do.call("rbind",lapply(res,function(x) x$y.tum))

str(betaAdj)
 # num [1:5000, 1:630] 0.357 0 0.009 0.018 0 0 0.018 1 0.022 0.036 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg09248054" "cg25340711" "cg06443533" "cg16601494" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

#inferred normal
betaNormal<-do.call("rbind",lapply(res,function(x) x$y.norm))

#original beta values
betaOld<-do.call("rbind",lapply(res,function(x) x$y.orig))

q("no")
###END