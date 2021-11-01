#####======================================================================#####
### Short paper on beta calibration based on tumor content
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncology & Pathology

################################################################################
##create work directories in default location

##work
HOME<-"~/hdd1"
DATA<-"~/hdd1/luTnbc"
GIT<-"~/Documents/adjustBetas"

##home
HOME<-"I:/data/luTnbc"
DATA<-"I:/data/luTnbc"
GIT<-"F:/gitProjects/adjustBetas"

##Trim sample set to only include 235 with data on all levels
TrimTo235<-TRUE

################################################################################
##load required packages

if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer") }

library("RColorBrewer")

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges") }

library("GenomicRanges")

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap") }

library("pheatmap")

if(!requireNamespace("flexmix", quietly = TRUE)) {
  install.packages("flexmix") }

library("flexmix")

################################################################################
##load data

##main data
load(paste0(HOME,"/preprocess/platformOverlaps/20190319_workspace_850kWithGexAndGeneModels.RData"))

#object - "betaAdj" - purity adjusted betas
load(file=paste0(HOME,"/preprocess/adjustBetas/200918_data850k_760405x236_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaValues.RData"))

#object - "betaNorm" - inferred normal methylation
load(file=paste0(HOME,"/preprocess/adjustBetas/200918_data850k_760405x236_minfiNormalized_ringnerAdjusted_purityAdjusted_normalBetaValues.RData"))

##feature annotations - repetitive element overlaps not yet incorporated
#load(file=paste0(HOME,"/preprocess/annotateFeatures/20201021_object_noChromFilter_760405x45_expandedAnnotations.RData"))
load(file=paste0(HOME,"/preprocess/annotateFeatures/20210124_object_noChromFilter_760405x50_annoObjectWithEncodeCre.RData"))

##get extra sample annotations from Staaf et al. NatMed
load(file=paste0(HOME,"/preprocess/annotateSamples/20201020_object_tnbcClass.RData"))
rownames(tnbcClass)<-tnbcClass$TumorAssay

##TFBS-CpG overlaps matrix - separate as size is big - load if needed
#load(file=paste0(HOME,"/preprocess/annotateFeatures/20210124_object_noChromFilter_760405x340_encodeTfbsAnnotations.RData"))

ls()
#  [1] "annoObj"        "betaAdj"        "betaNew"        "betaNorm"      
#  [5] "clinAnno"       "geneCoords"     "getProm"        "gexAnno"       
#  [9] "gexCoords"      "gexFpkm"        "gexTpm"         "HOME"          
# [13] "makePromObject" "probeAnno"      "sampleAnno"     "sampleSets"    
# [17] "tnbcClass"      "TrimTo235"     

##objects are row/column-matched
# "annoObj" - main probe annotation object
# "betaAdj" - purity adjusted betas
# "betaNew" - basic normalized + Ringner adjusted data
# "betaNorm" - inferred normal methylation
# "clinAnno" - clinical annotations from Staaf et al. Nat Med
# "geneCoords" - GRanges with gene coodinates
# "getProm" - function for getting promoter-object
# "gexAnno" - gene annotations for GEX matrix
# "gexCoords" - GRanges with gene coodinates
# "gexFpkm" - FPKM GEX data from JS | Raw Fpkm rounded to 3 decimals
# "gexTpm" - FPKM -> TPM converted data | is round( log2(TPM+1) ,3). REF: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# "HOME" - Path to "luTnbc"-folder on local machine
# "makePromObject" - used with getProm 
# "probeAnno"  - GRanges with probe coodinates
# "sampleAnno" - sample matrix for methylation data
# "sampleSets" - sample sets on all available data levels including overlaps
# "tnbcClass" - TNBC-type and PAM50 for tumors

##GEX and Methylation matrices are not column matched. Both have 236 samples but only 235 of theese overlap on both levels.

str(sampleSets)
# List of 7
#  $ samplesAll    : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ samplesMethWgs: chr [1:236] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ samplesGexWgs : chr [1:236] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...
#  $ samplesMeth   : chr [1:236] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ samplesGex    : chr [1:236] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...
#  $ samplesWgs    : chr [1:252] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ samplesClin   : chr [1:237] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...

dim(clinAnno)
#[1] 237 194

str(tnbcClass)
# 'data.frame':   237 obs. of  3 variables:
#  $ TumorAssay: chr  "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...
#  $ TNBCtype  : chr  "UNS" "LAR" "BL1" "UNS" ...
#  $ PAM50_AIMS: chr  "Basal" "Her2" "Basal" "Basal" ...

all.equal(rownames(tnbcClass),rownames(clinAnno)) 
#[1] TRUE

################################################################################
##trim to only include select samples with data on all levels 

if(TrimTo235) { ##open TrimTo235

##Apply these steps when you want to look at intersecting samples between all data levels
  
  ##If you want to look only on one level you can skip this step and analyze all data on one level
  ##Methylation consensus clustering has been performed on all 236 Tumors with 850K data

sel<-sampleSets$samplesAll

##Meth
all.equal(colnames(betaNorm),colnames(betaAdj)) 
#[1] TRUE
all.equal(colnames(betaNorm),colnames(betaNew)) 
#[1] TRUE
all.equal(colnames(betaNorm),rownames(sampleAnno)) 
#[1] TRUE

##Gex
all.equal(colnames(gexFpkm),colnames(gexTpm)) 
#[1] TRUE

##Gex-Meth
all.equal(colnames(gexFpkm),colnames(betaNorm)) 
#[1] "235 string mismatches"
sum(colnames(gexFpkm) %in% colnames(betaNorm)) 
#[1] 235

##Clinical
all(colnames(gexFpkm) %in% rownames(clinAnno)) 
#[1] TRUE
all(colnames(betaNorm) %in% rownames(clinAnno)) 
#[1] TRUE
all.equal(rownames(tnbcClass),rownames(clinAnno)) 
#[1] TRUE

length(sel)
#[1] 235

head(sel)
#[1] "PD31028a" "PD31029a" "PD31030a" "PD31031a" "PD31032a" "PD31033a"

##trim meth
betaNorm<-betaNorm[,sel]
betaAdj<-betaAdj[,sel]
betaNew<-betaNew[,sel]
sampleAnno<-sampleAnno[sel,]
##trim gex
gexFpkm<-gexFpkm[,sel]
gexTpm<-gexTpm[,sel]
##trim clinical 
clinAnno<-clinAnno[sel,]
tnbcClass<-tnbcClass[sel,]

##Meth
all.equal(colnames(betaNorm),colnames(betaAdj)) 
#[1] TRUE
all.equal(colnames(betaNorm),colnames(betaNew)) 
#[1] TRUE
all.equal(colnames(betaNorm),rownames(sampleAnno)) 
#[1] TRUE

##Gex
all.equal(colnames(gexFpkm),colnames(gexTpm)) 
#[1] TRUE

##Gex-Meth
all.equal(colnames(gexFpkm),colnames(betaNorm)) 
#[1] [1] TRUE

##Clinical
all.equal(colnames(gexFpkm),rownames(clinAnno)) 
#[1] TRUE
all.equal(colnames(betaNorm),rownames(clinAnno)) 
#[1] TRUE
all.equal(rownames(tnbcClass),rownames(clinAnno)) 
#[1] TRUE

dim(betaNorm)
#[1] 760405    235

dim(gexTpm)
#[1] 18776   235

rm(sel)

} ##close TrimTo235

################################################################################
##remove non CpG-probes from analyses as well as remove X/Y (or customize probe set)

##use annotations to focus analyses on different CpG-sets/contexts
str(annoObj) 
# 'data.frame':   760405 obs. of  50 variables:
#  $ illuminaID                  : chr  "cg21870274" "cg09499020" "cg16535257" "cg06325811" ...
#  $ chr                         : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start                       : num  69591 817115 817316 860948 864703 ...
#  $ end                         : num  69592 817116 817317 860949 864704 ...
#  $ hasUCSCknownGeneOverlap     : num  1 0 0 0 0 1 1 0 0 0 ...
#  $ nameUCSCknownGeneOverlap    : chr  "OR4F5" "" "" "" ...
#  $ numberUCSCknownGeneOverlap  : num  1 0 0 0 0 1 1 0 0 0 ...
#  $ hasUCSCknownGeneUp5kbOverlap: num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneDn5kbOverlap: num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxPromOverlap   : num  0 1 1 0 0 1 1 0 0 0 ...
#  $ hasUCSCknownTxUp5kbOverlap  : num  0 0 0 0 0 1 0 0 0 0 ...
#  $ hasUCSCknownTxDn5kbOverlap  : num  1 1 0 1 1 0 1 0 0 0 ...
#  $ hasGeneOverlap              : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ nameGeneOverlap             : chr  "79501|OR4F5" "" "" "" ...
#  $ numberGeneOverlap           : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasUp5kbOverlap             : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ nameUp5kbOverlap            : chr  "" "" "" "" ...
#  $ numberUp5kbOverlap          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasDn5kbOverlap             : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ nameDn5kbOverlap            : chr  "79501|OR4F5" "" "" "" ...
#  $ numberDn5kbOverlap          : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasPromOverlap              : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ namePromOverlap             : chr  "" "" "" "" ...
#  $ numberPromOverlap           : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ isPromMostVariable          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ namePromMostVariable        : chr  "" "" "" "" ...
#  $ numberPromMostVariable      : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasAtacOverlap              : num  0 0 0 0 0 1 0 0 0 0 ...
#  $ nameAtacOverlap             : chr  "" "" "" "" ...
#  $ numberAtacOverlap           : num  0 0 0 0 0 1 0 0 0 0 ...
#  $ isAtacMostVariable          : num  0 0 0 0 0 1 0 0 0 0 ...
#  $ isDistal                    : num  0 1 1 1 1 1 1 1 1 1 ...
#  $ isGenic                     : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ isDistalNonGenic            : num  0 1 1 1 1 1 1 1 1 1 ...
#  $ isGeneBody                  : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ weberOE                     : num  0.248 0.361 0.352 0.313 0.25 ...
#  $ weberClass                  : chr  "LCP" "ICP" "ICP" "ICP" ...
#  $ saxonovOE                   : num  0.238 0.302 0.342 0.277 0.268 ...
#  $ saxonovClass                : chr  "LCG" "LCG" "LCG" "LCG" ...
#  $ isCgi                       : num  0 0 0 0 0 1 0 0 0 0 ...
#  $ isShore                     : num  0 0 0 0 1 0 0 0 0 0 ...
#  $ isOcean                     : num  1 1 1 1 0 0 1 1 1 1 ...
#  $ cgiClass                    : chr  "ocean" "ocean" "ocean" "ocean" ...
#  $ featureClass                : chr  "proximal dn" "distal" "distal" "distal" ...
#  $ illuminaCpG_CpH_Probe       : chr  "cg" "cg" "cg" "cg" ...
#  $ encodeCre                   : num  0 1 1 1 0 1 0 0 0 0 ...
#  $ encodeCreCategory           : chr  "" "PLS" "PLS" "CTCF-only" ...
#  $ encodeCreCtcf               : num  0 1 1 1 0 1 0 0 0 0 ...
#  $ encodeChipMeanNcellPeaks    : num  0 1.5 1.29 0 0 ...
#  $ encodeChipUniqueTfPeaks     : num  0 8 14 0 0 30 0 9 7 7 ...

##keep autosomes, finter nonCpG-probes
sel<-which( annoObj$chr %in% paste("chr",c(1:22),sep="") &
  sub("(..).+","\\1",annoObj$illuminaID) == "cg"
  )

##if all combos match then all rows are matched..
all.equal(rownames(betaNorm),intersect(rownames(betaAdj),rownames(betaNew)) )
#[1] TRUE
all.equal(rownames(betaNorm),intersect(rownames(annoObj),names(probeAnno)) )
#[1] TRUE
all.equal(colnames(betaNorm),intersect(colnames(betaAdj),colnames(betaNew)) )
#[1] TRUE

length(sel)
#[1] 741145

betaNorm<-betaNorm[sel,]
betaAdj<-betaAdj[sel,]
betaNew<-betaNew[sel,]
annoObj<-annoObj[sel,]
probeAnno<-probeAnno[sel]

all.equal(rownames(betaNorm),intersect(rownames(betaAdj),rownames(betaNew)) )
#[1] TRUE
all.equal(rownames(betaNorm),intersect(rownames(annoObj),names(probeAnno)) )
#[1] TRUE
all.equal(colnames(betaNorm),intersect(colnames(betaAdj),colnames(betaNew)) )
#[1] TRUE

dim(betaNorm)
#[1] 741145

##check chromosomes+probe type post filter 
table(annoObj$chr,sub("(..).+","\\1",annoObj$illuminaID))
  #          cg
  # chr1  72114
  # chr10 36826
  # chr11 43686
  # chr12 39628
  # chr13 18228
  # chr14 26033
  # chr15 24885
  # chr16 32568
  # chr17 38908
  # chr18 13344
  # chr19 33141
  # chr2  57247
  # chr20 20642
  # chr21  8481
  # chr22 15777
  # chr3  43730
  # chr4  31990
  # chr5  39632
  # chr6  47131
  # chr7  40175
  # chr8  33818
  # chr9  23161

rm(sel)

################################################################################
### Ready to use - Switch to adjustBetas folder

ls()
#  [1] "annoObj"        "betaAdj"        "betaNew"        "betaNorm"      
#  [5] "clinAnno"       "geneCoords"     "getProm"        "gexAnno"       
#  [9] "gexCoords"      "gexFpkm"        "gexTpm"         "HOME"          
# [13] "makePromObject" "probeAnno"      "sampleAnno"     "sampleSets"    
# [17] "tnbcClass"      "TrimTo235"    

##use annoObj to focus on CpG subsets e.g. "atac distal"
#sel<- (annoObj$featureClass=="distal" | annoObj$featureClass=="distal body" ) & annoObj$hasAtacOverlap==1
##Then filter for e.g. variance to get top5000 or similar

rm(TrimTo235)

##home
HOME<-"I:/data"
DATA<-"I:/data/tcgaBrca"
GIT<-"F:/gitProjects/adjustBetas"

##create data directories
if( !file.exists( paste0(HOME,"/","adjustBetas2") )   ) {
  dir.create(paste0(HOME,"/","adjustBetas2/me"),recursive=TRUE)
}

HOME<-"~/hdd1/adjustBetas2"
HOME<-"I:/data/adjustBetas2"

################################################################################
##Make fig1

str(getProm("BRCA1"))
# List of 6
#  $ chrom     : chr "chr17"
#  $ gene      : chr "672|BRCA1"
#  $ tssStart  : int 43125483
#  $ tssDist   : int [1:30] 440 286 225 135 118 110 107 105 73 71 ...
#  $ probes    : chr [1:30] "cg19531713" "cg19088651" "cg08386886" "cg24806953" ...
#  $ probeStart: int [1:30] 43125042 43125196 43125257 43125347 43125364 43125372 43125375 43125377 43125409 43125411 ...

plot(clinAnno$ASCAT_TUM_FRAC,clinAnno$BATTENBERG_TUMOUR_FRAC,
  pch=16,xlab="ASCAT tumor fraction [WGS]",xlab="BATTENBERG tumor fraction [WGS]",
  )

tiff(paste0(HOME,"/fig1_conceptPlot.tiff"),width=12*500,height=12*500,units="px",res=500,compression="lzw")
par(fig=c(0,.3,.7,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=F)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro meth%",
  main="Lund TNBC cohort",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
legend("topleft",legend=c("Pyro neg","Pyro pos"),bty="n",col=1:2,pch=16,cex=1.25)
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

##Use flexmix to get pops
  #do modeling
  y2<-clinAnno$BRCA1_PromMetPc
  set.seed(12345)
  #y2<-y2+rnorm(length(y2),mean=1,sd=.05)
  x<-clinAnno$BATTENBERG_TUMOUR_FRAC 
  set.seed(12345)
  model <- stepFlexmix(y2 ~ x,k = 1:3, nrep = 10,verbose = FALSE)
  model <- getModel(model, "BIC")
  #get clusters
  cl<-clusters(model)

par(fig=c(.3,.6,.7,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro",
  main="Flexmix K=1,2,3 nrep=10, BIC",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.25)
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
points(y=clinAnno$BRCA1_PromMetPc[cl==3],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==3],pch=16,col=3)

##Use flexmix to get pops
  #do modeling
  y2<-clinAnno$BRCA1_PromMetPc
  set.seed(12345)
  y2<-y2+rnorm(length(y2),mean=0,sd=.25)
  x<-clinAnno$BATTENBERG_TUMOUR_FRAC 
  set.seed(12345)
  model <- stepFlexmix(y2 ~ x,k = 1:3, nrep = 10,verbose = FALSE)
  model <- getModel(model, "BIC")
  #get clusters
  cl<-clusters(model)

par(fig=c(.6,.9,.7,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,y2,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro + N(0,0.25)",
  main="Flexmix K=1,2,3 nrep=10, BIC",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.25)
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
points(y=clinAnno$BRCA1_PromMetPc[cl==3],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==3],pch=16,col=3)

# dev.off()

set.seed(12345)
# quantile(rnorm(length(y2),mean=0,sd=.25))
#          0%         25%         50%         75%        100% 
# -0.59508952 -0.14090821  0.03467275  0.20604706  0.66394707 

set.seed(12345)
mean(rnorm(length(y2),mean=0,sd=.25))
#[1] 0.03116207

set.seed(12345)
sd(rnorm(length(y2),mean=0,sd=.25))
#[1] 0.2649375

################################################################################
##Show line fits and intercepts

##panel D
par(fig=c(.1,.45,.35,.7),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro",
  main="Model \"normal breast\" methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

l1<-lm(clinAnno$BRCA1_PromMetPc[cl==1]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2<-lm(clinAnno$BRCA1_PromMetPc[cl==2]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
abline( l1,col=1,lwd=3)

points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
abline( l2,col=2,lwd=3)

legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.25)
text(x=.75,y=25,labels=paste0("L1 y=",round(coefficients(l1)[2],2),"x+",round(coefficients(l1)[1],2)),font=2,cex=1)
text(x=.75,y=15,labels=paste0("L2 y=",round(coefficients(l2)[2],2),"x+",round(coefficients(l2)[1],2)),font=2,cex=1)

##panel E
par(fig=c(.45,.8,.35,.7),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro",
  main="Inferred \"normal breast\" methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

l1<-lm(clinAnno$BRCA1_PromMetPc[cl==1]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2<-lm(clinAnno$BRCA1_PromMetPc[cl==2]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

y3<-clinAnno$BRCA1_PromMetPc
y3[cl==1]<-round(coefficients(l1)[1],2)+residuals(l1)
y3[cl==2]<-round(coefficients(l2)[1],2)+residuals(l2)
y3[y3>100]<-100
y3[y3<0]<-0

l1n<-lm(y3[cl==1]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2n<-lm(y3[cl==2]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

abline( l1,lwd=3,col="grey")
abline( l2,lwd=3,col="grey")

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col="grey")
points(y=y3[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col="grey")
points(y=y3[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)

abline( l1n,lwd=3,col=1)
abline( l2n,lwd=3,col=2)

legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.25)
text(x=.20,y=65,labels=paste0("L1 y=",round(coefficients(l1n)[2],2),"x+",round(coefficients(l1n)[1],2)),font=2,cex=1)
text(x=.20,y=55,labels=paste0("L2 y=",round(coefficients(l2n)[2],2),"x+",round(coefficients(l2n)[1],2)),font=2,cex=1)

##panel F
par(fig=c(.1,.45,0,.35),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)
plot(1-clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="1- WGS tumor fraction",ylab="BRCA1 pyro",
  main="Model pure tumor methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

l1<-lm(clinAnno$BRCA1_PromMetPc[cl==1]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2<-lm(clinAnno$BRCA1_PromMetPc[cl==2]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
abline( l1,col=1,lwd=3)

points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
abline( l2,col=2,lwd=3)

legend("topright",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.25)
text(x=.25,y=25,labels=paste0("L1 y=",round(coefficients(l1)[2],2),"x+",round(coefficients(l1)[1],2)),font=2,cex=1)
text(x=.25,y=10,labels=paste0("L2 y=",round(coefficients(l2)[2],2),"x+",round(coefficients(l2)[1],2)),font=2,cex=1)

##panel G
par(fig=c(.45,.8,0,.35),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)
plot(1-clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="1- WGS tumor fraction",ylab="BRCA1 pyro",
  main="Corrected tumor methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.2)

l1<-lm(clinAnno$BRCA1_PromMetPc[cl==1]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2<-lm(clinAnno$BRCA1_PromMetPc[cl==2]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

y3<-clinAnno$BRCA1_PromMetPc
y3[cl==1]<-round(coefficients(l1)[1],2)+residuals(l1)
y3[cl==2]<-round(coefficients(l2)[1],2)+residuals(l2)
y3[y3>100]<-100
y3[y3<0]<-0

l1n<-lm(y3[cl==1]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2n<-lm(y3[cl==2]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col="grey")
points(y=y3[cl==1],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
abline( l1,lwd=3,col="grey")
abline( l1n,lwd=3,col=1)

points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col="grey")
points(y=y3[cl==2],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
abline( l2,lwd=3,col="grey")
abline( l2n,lwd=3,col=2)

legend("right",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.25)
text(x=.20,y=25,labels=paste0("L1 y=",round(coefficients(l1n)[2],2),"x+",round(coefficients(l1n)[1],2)),font=2,cex=1)
text(x=.20,y=10,labels=paste0("L2 y=",round(coefficients(l2n)[2],2),"x+",round(coefficients(l2n)[1],2)),font=2,cex=1)

dev.off()


###END