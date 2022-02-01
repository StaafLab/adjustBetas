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
if( !file.exists( paste0(HOME,"/","adjustBetas_paper") )   ) {
  dir.create(paste0(HOME,"/","adjustBetas_paper/me"),recursive=TRUE)
}

HOME<-"~/hdd1/adjustBetas_paper"
HOME<-"I:/data/adjustBetas_paper"

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

# plot(clinAnno$ASCAT_TUM_FRAC,clinAnno$BATTENBERG_TUMOUR_FRAC,
#   pch=16,xlab="ASCAT tumor fraction [WGS]",ylab="BATTENBERG tumor fraction [WGS]",
#   )

tiff(paste0(HOME,"/plos_figure1.tiff"),width=15*300,height=15*300,units="px",res=300,compression="lzw")
#pdf(paste0(HOME,"/final_fig1.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,.3,.7,1),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=F)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro meth%",
  main="SCAN-B TNBC cohort",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
legend("topleft",legend=c("Pyro neg","Pyro pos"),bty="n",col=1:2,pch=16,cex=1.0)
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

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

par(fig=c(.3,.6,.7,1),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro",
  main="Flexmix K=1,2,3 nrep=10, BIC",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.0)
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

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

par(fig=c(.6,.9,.7,1),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,y2,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro + N(0,0.25)",
  main="Flexmix K=1,2,3 nrep=10, BIC",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.0)
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
points(y=clinAnno$BRCA1_PromMetPc[cl==3],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==3],pch=16,col=3)

#dev.off()

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
par(fig=c(.1,.45,.35,.7),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro",
  main="Model \"normal breast\" methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

l1<-lm(clinAnno$BRCA1_PromMetPc[cl==1]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2<-lm(clinAnno$BRCA1_PromMetPc[cl==2]~c(clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
abline( l1,col=1,lwd=3)

points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
abline( l2,col=2,lwd=3)

legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.0)
text(x=.4,y=20,labels=paste0("L1 y=",round(coefficients(l1)[2],2),"x+",round(coefficients(l1)[1],2)),font=2,cex=1,pos=4)
text(x=.4,y=10,labels=paste0("L2 y=",round(coefficients(l2)[2],2),"x+",round(coefficients(l2)[1],2)),font=2,cex=1,pos=4)

##panel E
par(fig=c(.45,.8,.35,.7),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=T)
plot(clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="WGS tumor fraction",ylab="BRCA1 pyro",
  main="Inferred \"normal breast\" methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

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

legend("topleft",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.0)
text(x=0.4,y=20,labels=paste0("L1 y=",round(coefficients(l1n)[2],2),"x+",round(coefficients(l1n)[1],2)),font=2,cex=1,pos=4)
text(x=0.4,y=10,labels=paste0("L2 y=",round(coefficients(l2n)[2],2),"x+",round(coefficients(l2n)[1],2)),font=2,cex=1,pos=4)

##panel F
par(fig=c(.1,.45,0,.35),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=T)
plot(1-clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="1- WGS tumor fraction",ylab="BRCA1 pyro",
  main="Model pure tumor methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

l1<-lm(clinAnno$BRCA1_PromMetPc[cl==1]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1]))
l2<-lm(clinAnno$BRCA1_PromMetPc[cl==2]~c(1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2]))

points(y=clinAnno$BRCA1_PromMetPc[cl==1],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==1],pch=16,col=1)
abline( l1,col=1,lwd=3)

points(y=clinAnno$BRCA1_PromMetPc[cl==2],x=1-clinAnno$BATTENBERG_TUMOUR_FRAC[cl==2],pch=16,col=2)
abline( l2,col=2,lwd=3)

legend("topright",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.0)
text(x=0.0,y=20,labels=paste0("L1 y=",round(coefficients(l1)[2],2),"x+",round(coefficients(l1)[1],2)),font=2,cex=1,pos=4)
text(x=0.0,y=10,labels=paste0("L2 y=",round(coefficients(l2)[2],2),"x+",round(coefficients(l2)[1],2)),font=2,cex=1,pos=4)

##panel G
par(fig=c(.45,.8,0,.35),mar=c(5.1,5.1,4.1,1.1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=T)
plot(1-clinAnno$BATTENBERG_TUMOUR_FRAC,clinAnno$BRCA1_PromMetPc,ylim=c(0,100),xlim=c(0,1),
  pch=16,xlab="1- WGS tumor fraction",ylab="BRCA1 pyro",
  main="Corrected tumor methylation",type="n",
  col=1+clinAnno$BRCA1_PromMetPc_Class,axes=F
  )
axis(1,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.0)
axis(2,at=seq(0,100,length.out=5),lwd=2,las=1,cex=1.0)

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

legend("right",legend=paste0("flexmix ",1:max(cl)),bty="n",col=1:max(cl),pch=16,cex=1.0)
text(x=0.0,y=20,labels=paste0("L1 y=",round(coefficients(l1n)[2],2),"x+",round(coefficients(l1n)[1],2)),font=2,cex=1,pos=4)
text(x=0.0,y=10,labels=paste0("L2 y=",round(coefficients(l2n)[2],2),"x+",round(coefficients(l2n)[1],2)),font=2,cex=1,pos=4)

dev.off()

fig<-image_scale(image_read(paste0(HOME,"/plos_figure1.tiff")),"2250x")

fig<-image_annotate(fig, "A", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)
fig<-image_annotate(fig, "B", size = 38, location=paste0("+",2250*.3+40,"+",40), color = "black",font="sans",weight=700)
fig<-image_annotate(fig, "C", size = 38, location=paste0("+",2250*.6+40,"+",40), color = "black",font="sans",weight=700)

fig<-image_annotate(fig, "D", size = 38, location=paste0("+",2250*.1+40,"+",2250*.3+40), color = "black",font="sans",weight=700)
fig<-image_annotate(fig, "E", size = 38, location=paste0("+",2250*.45+40,"+",2250*.3+40), color = "black",font="sans",weight=700)
fig<-image_annotate(fig, "F", size = 38, location=paste0("+",2250*.1+40,"+",2250*.65+40), color = "black",font="sans",weight=700)
fig<-image_annotate(fig, "G", size = 38, location=paste0("+",2250*.45+40,"+",2250*.65+40), color = "black",font="sans",weight=700)

image_write(fig,paste0(HOME,"/plos_figure1_final.tiff"))

################################################################################
##Perturb purity and estimate concordance with BRCA1-pyro

str(getProm())
# List of 6
#  $ chrom     : chr "chr17"
#  $ gene      : chr "672|BRCA1"
#  $ tssStart  : int 43125483
#  $ tssDist   : int [1:30] 440 286 225 135 118 110 107 105 73 71 ...
#  $ probes    : chr [1:30] "cg19531713" "cg19088651" "cg08386886" "cg24806953" ...
#  $ probeStart: int [1:30] 43125042 43125196 43125257 43125347 43125364 43125372 43125375 43125377 43125409 43125411 ...
 
which.max(apply(betaNew[getProm()$probes,],1,sd))
# cg09441966 
#         13 

##source - flexmix loaded on source
source(paste0(GIT,"/function_correctBetas.r"))

##get ascat tumor
fracA<-clinAnno[,"ASCAT_TUM_FRAC"]

##get battenberg tumor
fracB<-clinAnno[,"BATTENBERG_TUMOUR_FRAC"]

quantile(fracA)
#  0%  25%  50%  75% 100%
#0.10 0.29 0.44 0.60 0.92
quantile(fracB)
#       0%       25%       50%       75%      100%
#0.1134026 0.3628369 0.5144925 0.6699751 0.9270524

cor(fracA,fracB)
#[1] 0.9086183
cor(fracA,fracB,method="spe")
#[1] 0.9151116

fracTum<-round((fracA+fracB)/2,3)
rm(fracA,fracB)

##magnitude of change
unlist(lapply(seq(0.0,.5,by=.02),function(x) {
  mean(unlist(lapply(1:1000,function(y) mean(abs(fracTum-(fracTum-rnorm(length(fracTum),mean=0,sd=x)))) )) )
}))
#  [1] 0.00000000 0.01597357 0.03188370 0.04788066 0.06376187 0.07975185
#  [7] 0.09578673 0.11196906 0.12753086 0.14345450 0.15954574 0.17559817
# [13] 0.19130312 0.20749602 0.22317164 0.23915439 0.25507570 0.27099644
# [19] 0.28675604 0.30291128 0.31867986 0.33451630 0.35195350 0.36630076
# [25] 0.38345657 0.40048723

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

resAcc<-matrix(nrow=50,ncol=26)
resSens<-matrix(nrow=50,ncol=26)
resSpec<-matrix(nrow=50,ncol=26)
resBeta<-matrix(nrow=50,ncol=26)

for(i in 1:ncol(resAcc)) {
  cat(i,"\n")
  for(j in 1:nrow(resAcc)) {
  cat(".")
  fracTum2<-fracTum+rnorm(length(fracTum),mean=0,sd=seq(0.0,.5,by=.02)[i])
  fracTum2[fracTum2>1]<-1
  fracTum2[fracTum2<0]<-0
  b<-adjustBeta(methylation=betaNew["cg09441966",],purity=fracTum2,snames=colnames(betaNew),seed=FALSE)
    #function(methylation=NULL,purity=NULL,snames=NULL,nmax=3,nrep=3,seed=TRUE)

  gr.n<-b$groups
  gr.n<-1+ as.integer(gr.n == which.max(unlist(lapply(1:3,function(x) sum(gr.n[clinAnno$BRCA1_PromMetPc_Class==1]==x)))))
  #table(gr.n,clinAnno$BRCA1_PromMetPc_Class+1)
  #plot(fracTum2,betaNew["cg09441966",],col=1+clinAnno$BRCA1_PromMetPc_Class,pch=16)
  #plot(fracTum2,betaNew["cg09441966",],col=b$groups,pch=16)

  resAcc[j,i]<-getAcc(gr.n,clinAnno$BRCA1_PromMetPc_Class+1)

  resSens[j,i]<-getSens(gr.n,clinAnno$BRCA1_PromMetPc_Class+1)

  resSpec[j,i]<-getSpec(gr.n,clinAnno$BRCA1_PromMetPc_Class+1)

  resBeta[j,i]<-b$avg.betaDiff
  }
  cat("\n")
}

#pdf(paste0(HOME,"/final_fig6_panel_b.pdf"),width=8,height=8,useDingbats=FALSE)
tiff(paste0(HOME,"/plos_figure6_panel_b.tiff"),width=7.5*300,height=7.5*300,units="px",res=300,compression="lzw")

par(mar=c(5.1, 6.1, 4.1 ,0.1),fig=c(0,1,0,1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=FALSE)
plot(fracTum,betaNew["cg09441966",],col=clinAnno$BRCA1_PromMetPc_Class+1,pch=16,
  xlab="WGS purity estimate",ylab="BRCA1 cg09441966",xlim=c(0,1),ylim=c(0,1),
  axes=F
)
legend("topleft",legend=c("Pyro unMe","Pyro Me"),pch=16,col=1:2,bty="n")
axis(2,las=1)
axis(1,las=1)

set.seed(12345)
accc<-NULL
for( i in 1:10) {
  b<-adjustBeta(methylation=betaNew["cg09441966",],purity=fracTum,snames=colnames(betaNew),seed=FALSE)
  gr.n<-b$groups
  gr.n<-1+ as.integer(gr.n == which.max(unlist(lapply(1:3,function(x) sum(gr.n[clinAnno$BRCA1_PromMetPc_Class==1]==x)))))
  accc<-c(accc,getAcc(gr.n,clinAnno$BRCA1_PromMetPc_Class+1))
}
text(.7,.25,labels=paste0("Mean of 10 repetitions\n accuracy: ",signif(mean(accc),3),
    "\n range(",paste(signif(range(accc),3),collapse="-"),")"),font=2)

dev.off()

image_write(image_scale(image_read(paste0(HOME,"/plos_figure6_panel_b.tiff")),"1125x"),paste0(HOME,"/plos_figure6_panel_b.tiff"))
image_annotate(fig, "A", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)


#pdf(paste0(HOME,"/final_fig6_panel_c.pdf"),width=8,height=8,useDingbats=FALSE)
tiff(paste0(HOME,"/plos_figure6_panel_c.tiff"),width=7.5*300,height=7.5*300,units="px",res=300,compression="lzw")

par(mar=c(5.1, 6.1, 4.1 ,0.1),fig=c(0,1,0,1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=FALSE)
plot(1,pch=16,type="n",
  xlab="SD of noise",ylab="Accuracy vs Pyro",xlim=c(0,.5),ylim=c(0.9,1),
  axes=F
)
points(seq(0.0,.5,by=.02),colMeans(resAcc),pch=16)
for(i in 1:ncol(resAcc) ) {
  lines(rep(seq(0.0,.5,by=.02)[i],2),quantile(resAcc[,i],c(.025,.975)))
}
axis(2,las=1)
axis(1,las=1,at=seq(0.0,.5,by=.05))
text(.2,.92,labels="Mean + CI95 50 repetitions per bin",font=2)

dev.off()

image_write(image_scale(image_read(paste0(HOME,"/plos_figure6_panel_c.tiff")),"1125x"),paste0(HOME,"/plos_figure6_panel_c.tiff"))

#pdf(paste0(HOME,"/final_fig6_panel_d.pdf"),,width=8,height=8,useDingbats=FALSE)
tiff(paste0(HOME,"/plos_figure6_panel_d.tiff"),width=7.5*300,height=7.5*300,units="px",res=300,compression="lzw")

par(mar=c(5.1, 6.1, 4.1 ,0.1),fig=c(0,1,0,1),family="sans",font=2,font.axis=2,font.lab=2,font.sub=2,cex=1,cex.main=1,cex.sub=1,cex.lab=1,cex.axis=1,ps=18,new=FALSE)
plot(1,type="n",
  xlab="SD of noise",ylab="Mean beta shift post-correction",xlim=c(0,.5),ylim=c(0,-.15),
  axes=F
)
points(seq(0.0,.5,by=.02),colMeans(resBeta),pch=16)
for(i in 1:ncol(resAcc) ) {
  lines(rep(seq(0.0,.5,by=.02)[i],2),quantile(resBeta[,i],c(.025,.975)))
}
axis(2,las=1)
axis(1,las=1,at=seq(0.0,.5,by=.05))
text(.3,-.12,labels="Mean + CI95 50 repetitions per bin",font=2)

dev.off()

image_write(image_scale(image_read(paste0(HOME,"/plos_figure6_panel_d.tiff")),"1125x"),paste0(HOME,"/plos_figure6_panel_d.tiff"))

# ##first panel
a1<-image_read(paste0(HOME,"/plos_figure6_panel_a.tiff"))
a1<-image_annotate(a1, "A", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)

a2<-image_read(paste0(HOME,"/plos_figure6_panel_b.tiff"))
a2<-image_annotate(a2, "B", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)

a3<-image_read(paste0(HOME,"/plos_figure6_panel_c.tiff"))
a3<-image_annotate(a3, "C", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)

a4<-image_read(paste0(HOME,"/plos_figure6_panel_d.tiff"))
a4<-image_annotate(a4, "D", size = 38, location=paste0("+",40,"+",40), color = "black",font="sans",weight=700)

out1<-image_append(c(a1,a2),stack = F)
out2<-image_append(c(a3,a4),stack = F)
out<-image_append(c(out1,out2),stack = T)
out<-image_scale(out,"2250x")
image_write(out, path = paste0(HOME,"/plos_figure6_final.tiff"), format = "tiff")

save.image(file=paste0(HOME,"/figures_workspace_luTnbc.RData"))

q("no")
###END