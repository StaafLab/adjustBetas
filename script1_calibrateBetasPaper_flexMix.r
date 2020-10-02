#####======================================================================#####
### Short paper on beta calibration based on tumor content
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncology & Pathology

################################################################################
##create work directories in default location

##work
dir.create("~/hdd1/adjustBetas")
HOME<-"~/hdd1/adjustBetas"
GIT<-"~/Documents/adjustBetas"

##home
 dir.create("I:/data/adjustBetas")
 HOME<-"I:/data/adjustBetas"
 GIT<-"F:/gitProjects/adjustBetas"

##home
#setwd("I:/Lund University/Staaf_lab - Documents/tnbc_project/mattias/")

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
##load data

##Load GSE67919 - 96 Normal breast samples
load(paste0(HOME,"/data/","GSE67919_Annotations.RData"))
load(paste0(HOME,"/data/","GSE67919_Beta.RData"))

ls()
#[1] "adjustBeta"  "annotations" "beta"        "GIT"         "HOME" 

annotations_norm<-annotations
beta_norm<-beta
rm(annotations,beta)

load(paste0(HOME,"/data/","20191021_workspace_tcgaAtacBrca.RData"))

ls()
#  [1] "adjustBeta"       "annotations_norm" "atacManifest"     "atacObjects"     
#  [5] "beta_norm"        "betaNew"          "clinAnno"         "clusters.umap"   
#  [9] "geneCoords"       "getProm"          "gexAnno"          "gexCoords"       
# [13] "gexFpkm"          "gexTpm"           "GIT"              "HOME"            
# [17] "makePromObject"   "nonCodingObjects" "probeAnno"        "sampleAnno"      
# [21] "samples_use"      "sampleSets"       "stanfordIdToTcga" "tnbcClass"     

################################################################################
### Get tumor purity

str(sampleSets)
#List of 7
# $ samplesAll    : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
# $ samplesMethWgs: chr [1:236] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
# $ samplesGexWgs : chr [1:236] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...
# $ samplesMeth   : chr [1:236] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
# $ samplesGex    : chr [1:236] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...
# $ samplesWgs    : chr [1:252] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
# $ samplesClin   : chr [1:237] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...

all.equal(sampleSets$samplesMeth,colnames(betaNew))
#[1] TRUE

all.equal(sampleSets$samplesAll,samples_use)
#[1] TRUE

##get ascat tumor
fracA<-clinAnno[sampleSets$samplesAll,"ASCAT_TUM_FRAC"]

##get battenberg tumor
fracB<-clinAnno[sampleSets$samplesAll,"BATTENBERG_TUMOUR_FRAC"]

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

################################################################################
###get test data

##filter X and get top5k
betaData<-betaNew[,samples_use]

all.equal(sampleSets$samplesAll,colnames(betaData))
#[1] TRUE

table(as.character(seqnames(probeAnno[(rownames(betaData))])))
# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20
#72332 36967 43803 39746 18321 26109 24980 32643 38993 13371 33194 57496 20711
#chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY
# 8514 15817 43851 32123 39743 47249 40311 33945 23242 16640   304

betaData<-betaData[as.character(seqnames(probeAnno[(rownames(betaData))]))!="chrX",]
betaData<-betaData[as.character(seqnames(probeAnno[(rownames(betaData))]))!="chrY",]

table(as.character(seqnames(probeAnno[(rownames(betaData))])))
# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20
#72332 36967 43803 39746 18321 26109 24980 32643 38993 13371 33194 57496 20711
#chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9
# 8514 15817 43851 32123 39743 47249 40311 33945 23242

varF<-apply(betaData,1,sd)

varF<-varF >= quantile(varF,1-(5000/nrow(betaData)))

table(varF)
#varF
# FALSE   TRUE
#738461   5000

varF<-rownames(betaData)[varF]

length(intersect(varF,rownames(beta_norm)))
#[1] 3694

length(intersect(rownames(betaData),rownames(beta_norm)))
#[1] 383671

length(intersect(rownames(betaNew),rownames(beta_norm)))
#[1] 392729

##create data set and run function
set.seed(20191203)
testDat2<-betaData[varF,]

str(testDat2)
 # num [1:5000, 1:235] 0.008 0.044 0.163 0.705 0.65 0.067 0.057 0.704 0.06 0 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg06712559" "cg17928920" "cg09248054" "cg27541454" ...
 #  ..$ : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...

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

clusterSetRNGStream(cl, 20200918) ##will not make exactly replicable..
system.time(
	res<-parRapply(cl = cl, testDat2, adjustBeta,purity=fracTum,snames=colnames(testDat2))
)

table(unlist(lapply(res,function(x) x$n.groups)))
   # 1    2    3 
   # 1  790 4209 

#save(res,file=paste0(HOME,"/2020918_object_adjustedBetas.RData"))

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

temp4<-temp4[!apply(temp4,1,function(x) any(is.na(x))),]

quantile(testDat2)
#   0%   25%   50%   75%  100%
#0.000 0.111 0.409 0.686 1.000

quantile(temp4)
#    0%   25%   50%   75%  100% 
# 0.000 0.035 0.496 0.964 1.000 

##Plot histogram of betas before and after correction
pdf(paste0(HOME,"/20191203_top5k_betaDistribution_tumors_beforeAfterCorrection.pdf"),width=8,height=8,useDingbats=F)

par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot(1,xlim=range(round(density(temp4)$x,1)),ylim=range(round(density(temp4)$y,1))+c(0,.5),type="n",las=1,axes=F,
  xlab="beta",ylab="density"
)
lines(density(temp4),col=2,lwd=2)
lines(density(testDat2),col=1,lwd=2)
axis(1,lwd=2,las=1,at=seq(0,1,by=.2))
axis(2,lwd=2,las=1)
legend("topright",legend=c("unadjusted beta","adjusted beta"),col=c(1,2),lwd=2,bty="n")
dev.off()

##Plot histogram of betas before and after correction
pdf(paste0(HOME,"/20191203_top5k_correlationTumorFrac_beforeCorrection.pdf"),width=8,height=8,useDingbats=F)
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

################################################################################      ##HÄR
###plot top5k clusters - unadjuster order

# testDat2<-betaNew[varF,samples_use]

# str(testDat2)
# # num [1:5000, 1:235] 0.008 0.044 0.163 0.705 0.65 0.067 0.057 0.704 0.06 0 ...
# # - attr(*, "dimnames")=List of 2
# #  ..$ : chr [1:5000] "cg06712559" "cg17928920" "cg09248054" "cg27541454" ...
# #  ..$ : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...

##check stats
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
#rownames(temp1)<-rownames(testDat2)
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
#rownames(temp2)<-rownames(testDat2)

table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000

##remove NA-rows
testDat2<-testDat2[!apply(temp1,1,function(x) any(is.na(x))),]
temp2<-temp2[!apply(temp1,1,function(x) any(is.na(x))),]
temp1<-temp1[!apply(temp1,1,function(x) any(is.na(x))),]

table(seqnames(probeAnno[(rownames(temp1))]))
# chr1 chr22  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9 chr10 chr11 chr12
#  565    64   415   254   160   355   500   341   366    99   304   201   222
#chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21  chrX  chrY
#  108   126    94   132   209    95   234   108    48     0     0

testDat<-testDat2[,samples_use]
temp1<-temp1[,samples_use]
temp2<-temp2[,samples_use]

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(testDat) ),method="ward.D"),5)
c2<-unique(c1[hclust( as.dist( 1-cor(testDat) ),method="ward.D")$order])
r1<-hclust( dist(testDat),method="ward.D")
c3<-hclust( as.dist( 1-cor(testDat) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)

sample_anno<-data.frame(unadj5000=c1,
  adj5000=c4,
  AIMS=tnbcClass$PAM50_AIMS,
  TNBC=tnbcClass$TNBCtype,
  #umap=factor(clusters.umap[samples_use,"class"]),
  hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-samples_use
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    #umap = c("1" = "#5977ff", "2" = "#f74747"),
    AIMS = c("Basal" = "red" , "Her2" = "pink" , "LumA" = "darkgreen" , "LumB" = "orange" , "Normal" = "grey"),
    TNBC = c("BL1"="#E41A1C","BL2"="#377EB8","IM"="#4DAF4A","LAR"="#984EA3","M"="#FF7F00","MSL"="#FFFF33","NA"="#666666","UNS"="#A65628")
    ,hrd3 = c("[0.0,0.2)" ="#FEE0D2" , "[0.2,0.7)" ="#FC9272" ,"[0.7,1.0]"="#EF3B2C" )
  )

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_unadjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(testDat,clustering_distance_rows = "euclidean",#"correlation"
  clustering_distance_cols = "correlation", clustering_method = "ward.D",show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_adjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_normalBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

##infiltration estimeates by group
tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_unadjBeta_forInfiltrationEstimate.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(testDat,clustering_distance_rows = "euclidean",#"correlation"
  clustering_distance_cols = "correlation", clustering_method = "ward.D",show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_adjBeta_forInfiltrationEstimate.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3,
  show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_normBeta_forInfiltrationEstimate.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3,
  show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", unadj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_unadj_ascat_battenberg_forInfiltrationEstimate.tiff"),width=6*500,height=8*500,units="px",res=500,compression="lzw")
par(mfrow=c(2,1),font.lab=2,font=2,lwd=2,font.axis=2)
boxplot(clinAnno[names(c1),"BATTENBERG_TUMOUR_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Battenberg tumor%"
  )
boxplot(clinAnno[names(c1),"ASCAT_TUM_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Ascat tumor%"
  )
dev.off()

pdf(paste0(HOME,"/20191203_top5k_heatmap_unadj_ascat_battenberg_forInfiltrationEstimate.pdf"),width=12,height=12,useDingbats=F)
par(mfrow=c(2,1),font.lab=2,font=2,lwd=2,font.axis=2)
boxplot(clinAnno[names(c1),"BATTENBERG_TUMOUR_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Battenberg tumor%",las=1
  )
boxplot(clinAnno[names(c1),"ASCAT_TUM_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Ascat tumor%",las=1
  )
dev.off()

library(magick)
a1<-image_read(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_unadjBeta_forInfiltrationEstimate.tiff"))
a2<-image_read(paste0(HOME,"/20191203_top5k_heatmap_unadj_ascat_battenberg_forInfiltrationEstimate.tiff"))
a3<-image_read(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_adjBeta_forInfiltrationEstimate.tiff"))
a4<-image_read(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_unadjClust_normBeta_forInfiltrationEstimate.tiff"))

image_write(image_scale(image_append(c(a1,a2)),5000), path = paste0(HOME,"/20191203_top5k_unadjClust_unadjBeta_combinedIfiltration.tiff"), format = "tiff")
image_write(image_scale(image_append(c(a3,a2)),5000), path = paste0(HOME,"/20191203_top5k_unadjClust_adjBeta_combinedIfiltration.tiff"), format = "tiff")
image_write(image_scale(image_append(c(a4,a2)),5000), path = paste0(HOME,"/20191203_top5k_unadjClust_normBeta_combinedIfiltration.tiff"), format = "tiff")

rm(a1,a2,a3,a4)

rm(sample_anno,my_colour,c1,c2,c3,c4,r1)

################################################################################
###plot top5000 clusters - adjusted order

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-unique(c1[hclust( as.dist( 1-cor(temp1) ),method="ward.D")$order])
r1<-hclust( dist(temp1),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(testDat) ),method="ward.D"),5)

sample_anno<-data.frame(adj5000=c1,
  unadj5000=c4,
  AIMS=tnbcClass$PAM50_AIMS,
  TNBC=tnbcClass$TNBCtype,
  #umap=factor(clusters.umap[samples_use,"class"]),
  hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-samples_use
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    #umap = c("1" = "#5977ff", "2" = "#f74747"),
    AIMS = c("Basal" = "red" , "Her2" = "pink" , "LumA" = "darkgreen" , "LumB" = "orange" , "Normal" = "grey"),
    TNBC = c("BL1"="#E41A1C","BL2"="#377EB8","IM"="#4DAF4A","LAR"="#984EA3","M"="#FF7F00","MSL"="#FFFF33","NA"="#666666","UNS"="#A65628")
    ,hrd3 = c("[0.0,0.2)" ="#FEE0D2" , "[0.2,0.7)" ="#FC9272" ,"[0.7,1.0]"="#EF3B2C" )
  )

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp1,clustering_distance_rows = "euclidean",#"correlation"
  clustering_distance_cols = "correlation", clustering_method = "ward.D",show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_unadjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(testDat,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, adj clust ,pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_normalBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

##infiltration estimeates by group
tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_adjBeta_forInfiltrationEstimate.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp1,clustering_distance_rows = "euclidean",#"correlation"
  clustering_distance_cols = "correlation", clustering_method = "ward.D",show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_unadjBeta_forInfiltrationEstimate.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(testDat,cluster_rows = r1, cluster_cols = c3,
  show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_normBeta_forInfiltrationEstimate.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3,
  show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, \"inferred normal\", adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_ascat_battenberg_forInfiltrationEstimate.tiff"),width=6*500,height=8*500,units="px",res=500,compression="lzw")
par(mfrow=c(2,1),font.lab=2,font=2,lwd=2,font.axis=2)
boxplot(clinAnno[names(c1),"BATTENBERG_TUMOUR_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Battenberg tumor%"
  )
boxplot(clinAnno[names(c1),"ASCAT_TUM_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Ascat tumor%"
  )
dev.off()

pdf(paste0(HOME,"/20191203_top5k_heatmap_ascat_battenberg_forInfiltrationEstimate.pdf"),width=12,height=12,useDingbats=F)
par(mfrow=c(2,1),font.lab=2,font=2,lwd=2,font.axis=2)
boxplot(clinAnno[names(c1),"BATTENBERG_TUMOUR_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Battenberg tumor%",las=1
  )
boxplot(clinAnno[names(c1),"ASCAT_TUM_FRAC"]~factor(c1,levels=c2),varwidth=T,ylim=c(0,1),
  col=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00")[c2],ylab="Ascat tumor%",las=1
  )
dev.off()

library(magick)
a1<-image_read(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_adjBeta_forInfiltrationEstimate.tiff"))
a2<-image_read(paste0(HOME,"/20191203_top5k_heatmap_ascat_battenberg_forInfiltrationEstimate.tiff"))
a3<-image_read(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_unadjBeta_forInfiltrationEstimate.tiff"))
a4<-image_read(paste0(HOME,"/20191203_top5k_heatmap_pear_eucl_adjClust_normBeta_forInfiltrationEstimate.tiff"))

image_write(image_scale(image_append(c(a1,a2)),5000), path = paste0(HOME,"/20191203_top5k_adjClust_adjBeta_combinedIfiltration.tiff"), format = "tiff")
image_write(image_scale(image_append(c(a3,a2)),5000), path = paste0(HOME,"/20191203_top5k_adjClust_unadjBeta_combinedIfiltration.tiff"), format = "tiff")
image_write(image_scale(image_append(c(a4,a2)),5000), path = paste0(HOME,"/20191203_top5k_adjClust_normBeta_combinedIfiltration.tiff"), format = "tiff")

rm(a1,a2,a3,a4)

rm(sample_anno,my_colour,c1,c2,c3,c4,r1)

################################################################################
###Check correlation "inferred normal" to true normal

ff<-intersect( rownames(temp2) , rownames(beta_norm) )

length(ff)
#[1] 3694

#plot( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]) )

table(apply(beta_norm[ff,],1,function(z) any(is.na(z))))
#FALSE  TRUE 
# 2902   792 

table(apply(beta_norm[ff,],1,function(z) sum(is.na(z))))
#   0    1    2    3    4    5    6    7    8    9   10   11   12   14   16   17   21   24 
#2902  344  168  125   54   36   25   12   11    4    4    2    2    1    1    1    1    1 

ff2<-!apply(beta_norm[ff,],1,function(z) any(is.na(z)))

ff<-ff[ff2]

plot( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]) )

##will change slightly if rerun - not deterministic..
cor( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]) )
#[1] 0.763079
(cor( rowMeans(temp2[ff,]),rowMeans(temp1[ff,]),method="spe" ))
#[1] 0.100357
(cor( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]),method="spe" ))
#[1] 0.6649534
(cor( rowMeans(temp2[ff,]),rowMeans(temp1[ff,]),method="pe" ))
#[1] -0.03620098
(sf<-cor( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]),method="pe" ))
#[1] 0.763079
(fs<-cor.test( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]),method="pe" )$p.value)
#[1] 0

length(ff)
#[1] 2902

pdf(paste0(HOME,"/20191203_top5kBySd_betaNormals_inferredVsActual.pdf"),width=8,height=8,useDingbats=F)
par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot( rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]),pch=16
  ,main="Correlation 450K normal - 850K inferred normal"
  ,xlab="mean inferred normal beta 850k",ylab="mean normal beta 450k GSE67919"
  ,type="n",las=1,axes=F,xlim=c(0,1),ylim=c(0,1)
 )
points(rowMeans(temp2[ff,]),rowMeans(beta_norm[ff,]),pch=16)
text(.1,.9,paste0("r=",round(sf,2)," | p<2.2e-16 \n",length(ff)," CpGs"))
abline(lm(rowMeans(beta_norm[ff,])~rowMeans(temp2[ff,])),lwd=2,col=2)
axis(1,lwd=2,las=1,at=seq(0,1,by=.2))
axis(2,lwd=2,las=1,at=seq(0,1,by=.2))
dev.off()

rm(ff,fs,ff2)

################################################################################
###Save objects

##Save correction object
correctionTop5000<-res

length(correctionTop5000)
#[1] 5000

save(correctionTop5000,file=paste0(HOME,"/20191203_top5k_correctionObject_5000cpgs_235tumors.RData"))

dataTop5000<-testDat
save(dataTop5000,file=paste0(HOME,"/20191203_top5k_testDataBetaMatrix_235withGex.RData"))

dataAdjTop5000<-temp1
save(dataAdjTop5000,file=paste0(HOME,"/20191203_top5k_testDataAdjBetaMatrix_235withGex.RData"))

rm(correctionTop5000,dataTop5000,dataAdjTop5000)

################################################################################
### Do comparative plots

##check stats
testDat<-do.call("rbind",lapply(res,function(x) x$y.orig))

##alluvial of results
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
unique(c1[hclust( as.dist( 1-cor(temp1) ),method="ward.D")$order])
#[1] 4 3 2 5 1

c1<-cutree( hclust( as.dist( 1-cor(testDat) ),method="ward.D"),5)
unique(c1[hclust( as.dist( 1-cor(testDat) ),method="ward.D")$order])
#[1] 3 4 5 2 1

cl<-as.data.frame(cbind(unadjusted=cutree( hclust( as.dist( 1-cor(testDat) ),method="ward.D"),5),
  adjusted=letters[1:5][cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)]
  ),stringsAsFactors=F)

table(cl$unadjusted,cl$adjusted)
  #    a  b  c  d  e
  # 1 23 51  0  0  6
  # 2 46  6 10  0 14
  # 3  0  0  1 13  0
  # 4  0  0 14  8  0
  # 5 16  6  2  0 19

cl<-as.data.frame(table(cl$unadjusted,cl$adjusted))
cl$Var1<-factor(cl$Var1,levels=c("3","4","5","2","1"))
cl$Var2<-factor(cl$Var2,levels=c("d","c","b","e","a"))

(pal<-brewer.pal(5,"Set1")[c(4 ,3 ,2,5 ,1)])
#[1] "#4DAF4A" "#FF7F00" "#377EB8" "#984EA3" "#E41A1C"
#names(pal)<-c("a","b","c","d","e")
(pal2<-brewer.pal(5,"Set1")[c(3,4,5,2,1)])
#[1] "#4DAF4A" "#984EA3" "#FF7F00" "#E41A1C" "#377EB8"

pal<-rev(pal)
pal2<-rev(pal2)

q<-ggplot(cl,
       aes(y = Freq,
           axis1 = Var1, axis2 = Var2,)) +
  geom_alluvium(aes(fill = Var1 ),alpha=.85,       #fill = Var1
                width = 0, knot.pos = 1/5, reverse = T,inherit.aes = T
                ) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/20, reverse = TRUE,colour="black",fill=c(pal2,pal)) +     #colour=c(pal2,pal),fill=c(pal2,pal)
  geom_text(stat = "stratum", infer.label = F, reverse = TRUE,size=10,fontface="bold",label=c(rev(c("c","d","e","b","a")),(c("d","c","b","e","a")))
    ) +
  scale_x_continuous(breaks = 1:2, labels = c("unadjusted", "adjusted")) +
  scale_fill_manual(values=rev(pal2)) +
  #scale_fill_manual(values=rep("lightgrey",5)) +
  #coord_flip() +
  #restitle("") +
  theme_bw() +
  theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    )

pdf(paste0(HOME,"/20191203_top5k_alluvial_pear_eucl_hClust_unadjBeta_to_adjBeta.pdf"),width=12,height=12,useDingbats=F)
q
dev.off()

rm(q,cl)

################################################################################
###Adjusted clustering

##check stats
testDat<-do.call("rbind",lapply(res,function(x) x$y.orig))

##leave out dendrograms

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-unique(c1[hclust( as.dist( 1-cor(temp1) ),method="ward.D")$order])
r1<-hclust( dist(temp1),method="ward.D")
c3<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(testDat) ),method="ward.D"),5)

c1<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c1)))))
c4<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c4)))))

sample_anno<-data.frame(adj5000=c1,
  unadj5000=c4,
  AIMS=tnbcClass$PAM50_AIMS
  #TNBC=tnbcClass$TNBCtype
  #umap=factor(clusters.umap[samples_use,"class"]),
  #hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-samples_use
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("a"="#E41A1C","b"="#377EB8","c"="#4DAF4A","d"="#984EA3","e"="#FF7F00"),
    adj5000=c("a"="#E41A1C","b"="#377EB8","c"="#4DAF4A","d"="#984EA3","e"="#FF7F00"),
    #umap = c("1" = "#5977ff", "2" = "#f74747"),
    AIMS = c("Basal" = "red" , "Her2" = "pink" , "LumA" = "darkgreen" , "LumB" = "orange" , "Normal" = "grey")
    #TNBC = c("BL1"="#E41A1C","BL2"="#377EB8","IM"="#4DAF4A","LAR"="#984EA3","M"="#FF7F00","MSL"="#FFFF33","NA"="#666666","UNS"="#A65628")
    #,hrd3 = c("[0.0,0.2)" ="#FEE0D2" , "[0.2,0.7)" ="#FC9272" ,"[0.7,1.0]"="#EF3B2C" )
  )

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp1,clustering_distance_rows = "euclidean",#"correlation"
  clustering_distance_cols = "correlation", clustering_method = "ward.D",show_rownames=F,show_colnames=F
  ,main=" ",cutree_cols=5,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_adjClust_unadjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(testDat,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main=" ",cutree_cols=5,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_adjClust_normalBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main=" ",cutree_cols=5,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_legend_pear_eucl_adjClust_normalBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",cutree_cols=5,fontsize=18,fontsize_col=16
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=T,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=T
)
dev.off()

rm(sample_anno,my_colour,c1,c2,c3,c4,r1)

dev.off()

################################################################################
###Unadjusted clustering

##do clustering
c1<-cutree( hclust( as.dist( 1-cor(testDat) ),method="ward.D"),5)
c2<-unique(c1[hclust( as.dist( 1-cor(testDat) ),method="ward.D")$order])
r1<-hclust( dist(testDat),method="ward.D")
c3<-hclust( as.dist( 1-cor(testDat) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)

c1<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c1)))))
c4<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c4)))))

sample_anno<-data.frame(unadj5000=c1,
  adj5000=c4,
  AIMS=tnbcClass$PAM50_AIMS
  #TNBC=tnbcClass$TNBCtype
  #umap=factor(clusters.umap[samples_use,"class"]),
  #hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-samples_use
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("a"="#E41A1C","b"="#377EB8","c"="#4DAF4A","d"="#984EA3","e"="#FF7F00"),
    adj5000=c("a"="#E41A1C","b"="#377EB8","c"="#4DAF4A","d"="#984EA3","e"="#FF7F00"),
    #adj5000=c("a"="red","b"="pink","c"="darkgreen","d"="orange","e"="grey"),
    #umap = c("1" = "#5977ff", "2" = "#f74747"),
    AIMS = c("Basal" = "red" , "Her2" = "pink" , "LumA" = "darkgreen" , "LumB" = "orange" , "Normal" = "grey")
    #TNBC = c("BL1"="#E41A1C","BL2"="#377EB8","IM"="#4DAF4A","LAR"="#984EA3","M"="#FF7F00","MSL"="#FFFF33","NA"="#666666","UNS"="#A65628")
    #,hrd3 = c("[0.0,0.2)" ="#FEE0D2" , "[0.2,0.7)" ="#FC9272" ,"[0.7,1.0]"="#EF3B2C" )
  )

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_unadjClust_unadjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(testDat,clustering_distance_rows = "euclidean",#"correlation"
  clustering_distance_cols = "correlation", clustering_method = "ward.D",show_rownames=F,show_colnames=F
  ,main="",cutree_cols=5,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_unadjClust_adjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",cutree_cols=5,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_unadjClust_normalBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",cutree_cols=5,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_legend_pear_eucl_unadjClust_normalBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",cutree_cols=5,fontsize=18,fontsize_col=16
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=T,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=T
)
dev.off()

rm(sample_anno,my_colour,c1,c2,c3,c4,r1)

################################################################################
###Combine all

library(magick)
a1<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_adjClust_unadjBeta.tiff"))
a2<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_adjClust_normalBeta.tiff"))
a3<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_adjClust_adjBeta.tiff"))

a4<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_unadjClust_unadjBeta.tiff"))
a5<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_unadjClust_normalBeta.tiff"))
a6<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_pear_eucl_unadjClust_adjBeta.tiff"))

a7<-image_read(paste0(HOME,"/20191203_top5k_heatmap_ascat_battenberg_forInfiltrationEstimate.tiff"))
a8<-image_read(paste0(HOME,"/20191203_top5k_heatmap_unadj_ascat_battenberg_forInfiltrationEstimate.tiff"))

a11<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_legend_pear_eucl_adjClust_normalBeta.tiff"))
a12<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_legend_pear_eucl_unadjClust_normalBeta.tiff"))

a11<-image_crop(a11,"1000x6000+3850")
a12<-image_crop(a12,"1000x6000+3850")

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_white.tiff"),width=.2*500,height=12*500,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a9<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_white.tiff"))

tiff(paste0(HOME,"/20191203_top5k_heatmap_noAnno_white2.tiff"),width=18*500,height=.4*500,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a10<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_white2.tiff"))

out<-image_append(c(a9,
  a1,
  a9,
  a2,
  a9,
  a3,
  a11,
  a9
  ),stack = F)
out<-image_scale(out,"9000x")

out2<-image_append(c(a9,
  a4,
  a9,
  a5,
  a9,
  a6,
  a12,
  a9
  ),stack = F)
out2<-image_scale(out2,"9000x")

out3<-image_append(c(out2,a10,out),stack = T)

image_write(out3, path = paste0(HOME,"/20191203_top5k_heatmap_noAnno_unadj_adj_combined.tiff"), format = "tiff")

rm(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,out,out2,out3)

gc()

################################################################################
###Do plot with pairwise correlations "inferred normal" to normal

ff<-intersect( rownames(temp2) , rownames(beta_norm) )

ff2<-!apply(beta_norm[ff,],1,function(z) any(is.na(z)))

ff<-ff[ff2]

pdf(paste0(HOME,"/20191203_top5kBySd_betaNormals_allInferredVsActual.pdf"),width=8,height=8,useDingbats=F)
par(font=2,font.axis=2,font.lab=2,font.sub=2)
plot( density( apply(temp2[ff,],2,function(x) cor(x,rowMeans(beta_norm[ff,]))) ),pch=16,cex=2,cex.lab=1.6,cex.main=2
  ,main="Sample correlations to average 450K normal"
  ,xlab="Pearson correlation",ylab="Density"
  ,type="n",las=1,axes=F,xlim=c(-.4,1)
 )
lines(density( apply(temp1[ff,],2,function(x) cor(x,rowMeans(beta_norm[ff,]))) ),lwd=3,col="grey")
lines(density( apply(temp2[ff,],2,function(x) cor(x,rowMeans(beta_norm[ff,]))) ),lwd=3,col=1)

abline(h=0,lwd=3,col=1)
abline(v=median(apply(temp1[ff,],2,function(x) cor(x,rowMeans(beta_norm[ff,])))),lwd=3,col="grey",lty=2)
abline(v=median(apply(temp2[ff,],2,function(x) cor(x,rowMeans(beta_norm[ff,])))),lwd=3,col=1,lty=2)

(sf<-median(apply(temp2[ff,],2,function(x) cor(x,rowMeans(beta_norm[ff,])))))
#[1] 0.6661728
text(-.2,3.5,paste0("median r=",round(sf,2)),cex=1.6)
legend("topleft",legend=c("tumor","inferred normal"),col=c("grey",1),lwd=3,bty="n",cex=1.6)

axis(1,lwd=2,las=1,at=seq(-.4,1,by=.2),cex.axis=1.6)
axis(2,lwd=2,las=1,cex.axis=1.6)
dev.off()

rm(sf)

################################################################################
###Do basal vs luminal split in 100 * 500 random CpG sets
  ##extract -logP(fisher) hclust 2-split vs PAM50 basal-luminal
  ##plot distribution of p-values for corrected vs uncorrected comparison

##number of lists
N<-100
resMat<-matrix(ncol=3,nrow=N,dimnames=list(1:N,c("p.raw","p.adj","p.dic")))

# ##generate 500 random CpG sets
# set.seed(20191214)
# varF<-apply(betaData[,samples_use],1,sd)
# varF<-varF > quantile(varF,.5)
# p_list<-lapply(1:N,function(x) sample(rownames(betaData)[varF],500) )

##generate 500 random CpG sets
set.seed(20200904)
varF<-apply(betaData[,samples_use],1,sd)
varF<-varF > 0
p_list<-lapply(1:N,function(x) sample(rownames(betaData)[varF],500) )

fracTum2<-fracTum[match(samples_use,colnames(betaData))] 

##do multicore
no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

clusterEvalQ(cl, {
  library("flexmix")
})

clusterSetRNGStream(cl, 20200918) ##will not make exactly replicable..

for (i in 1:N) {
 	cat(i," of ",N,"\n")
 	##betaData filtered for chrX/Y
 	b<-parRapply(cl = cl, betaData[p_list[[i]],samples_use], adjustBeta,purity=fracTum2,snames=samples_use)
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

  	resMat[i,1]<- -log10(fisher.test(table(b2,tnbcClass$PAM50_AIMS != "Basal"))$p.value)
  	resMat[i,2]<- -log10(fisher.test(table(b1,tnbcClass$PAM50_AIMS != "Basal"))$p.value)
  	resMat[i,3]<- -log10(fisher.test(table(b3,tnbcClass$PAM50_AIMS != "Basal"))$p.value)
}
rm(b,b1,b2,b3,i,N,cl)

save(p_list,file=paste0(HOME,"/20191214_top100percentBySd_basalVsLuminalSplitIn100randomSets_UsedProbeSets.RData"))
rm(p_list)

pdf(paste0(HOME,"/20191214_top100percentBySd_basalVsLuminalSplitIn100randomSets.pdf"),width=8,height=8,useDingbats=F)
par(font=2,font.axis=2,font.lab=2,font.sub=2)

plot(density(resMat[,2]),xlim=c(-5,max(resMat)+c(5)),pch=16,cex=2,cex.lab=1.6,cex.main=2
  ,main="Discrimination of PAM50 Basal vs Rest split"
  ,xlab="-log10(p,Fisher test), 100 iterations",ylab="Density"
  ,type="n",las=1,axes=F,
  ylim=c(0,max( c(density(resMat[,1])$y,density(resMat[,2])$y,density(resMat[,3])$y) ) )
)
lines(density(resMat[,1]),col="grey",lwd=3)
lines(density(resMat[,2]),col=1,lwd=3)
lines(density(resMat[,3]),col=2,lwd=3)

abline(v=median(resMat[,1]),lwd=3,col="grey",lty=2)
abline(v=median(resMat[,2]),lwd=3,col=1,lty=2)
abline(v=median(resMat[,3]),lwd=3,col=2,lty=2)
abline(h=0,lwd=3,col=1)

legend("topright",legend=c("unadjusted beta","adjusted beta","dichotomized beta"),col=c("grey",1,2),lwd=3,bty="n",cex=1.6)

axis(1,lwd=2,las=1,cex.axis=1.6)
axis(2,lwd=2,las=1,cex.axis=1.6)
dev.off()

save(resMat,file=paste0(HOME,"/20191214_top100percentBySd_basalVsLuminalSplitIn100randomSets_FisherPVals.RData"))
rm(resMat)

rm(varF)

################################################################################
###Do plot with one of the random 500 iterations

##4 figure panels
  ##1. heatmap uncorrected
  ##2. heatmap corrected
  ##3. alluvial 4-g
  ##.4 basal/luminal v 2-group

##choose one of the random iters
  ##best change pre-post?!
load(file=paste0(HOME,"/20191214_top100percentBySd_basalVsLuminalSplitIn100randomSets_FisherPVals.RData"))

load(file=paste0(HOME,"/20191214_top100percentBySd_basalVsLuminalSplitIn100randomSets_UsedProbeSets.RData"))

iii<-which.max( resMat[,2]-resMat[,1] )
iii<-p_list[[iii]]

##do multicore
no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

clusterEvalQ(cl, {
  library("flexmix")
})

clusterSetRNGStream(cl, 20200918) ##will not make exactly replicable..

b<-parRapply(cl = cl, betaData[iii,samples_use], adjustBeta,purity=fracTum2,snames=samples_use) 

rm(resMat,iii,p_list)

##adjusted
b1<-do.call("rbind",lapply(b,function(x) x$y.tum))

##unadjusted
b2<-do.call("rbind",lapply(b,function(x) x$y.norm))


##do clustering
c1<-cutree( hclust( as.dist( 1-cor(b1) ),method="ward.D"),2)
c2<-unique(c1[hclust( as.dist( 1-cor(b1) ),method="ward.D")$order])
r1<-hclust( dist(b1),method="ward.D") ##do rowclusters based on unadj data
c3<-hclust( as.dist( 1-cor(b1) ),method="ward.D")
c4<-cutree( hclust( as.dist( 1-cor(b2) ),method="ward.D"),2)
c5<-hclust( as.dist( 1-cor(b2) ),method="ward.D")

c1<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c1)))))
c4<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c4)))))

sample_anno<-data.frame(adj500=c1,
  unadj500=c4,
  AIMS=tnbcClass$PAM50_AIMS
  #TNBC=tnbcClass$TNBCtype
  #umap=factor(clusters.umap[samples_use,"class"]),
  #hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-samples_use
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj500=c("a"="#E41A1C","b"="#377EB8"),
    adj500=c("a"="#E41A1C","b"="#377EB8"),
    #adj5000=c("a"="red","b"="pink","c"="darkgreen","d"="orange","e"="grey"),
    #umap = c("1" = "#5977ff", "2" = "#f74747"),
    AIMS = c("Basal" = "red" , "Her2" = "pink" , "LumA" = "darkgreen" , "LumB" = "orange" , "Normal" = "grey")
    #TNBC = c("BL1"="#E41A1C","BL2"="#377EB8","IM"="#4DAF4A","LAR"="#984EA3","M"="#FF7F00","MSL"="#FFFF33","NA"="#666666","UNS"="#A65628")
    #,hrd3 = c("[0.0,0.2)" ="#FEE0D2" , "[0.2,0.7)" ="#FC9272" ,"[0.7,1.0]"="#EF3B2C" )
  )

tiff(paste0(HOME,"/20191215_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(b1,cluster_rows = r1,cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",fontsize=18,cutree_cols=2
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

tiff(paste0(HOME,"/20191215_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta_annotations.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
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

c1<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c1)))))
c4<-sub("5","e",sub("4","d",sub("3","c",sub("2","b",sub("1","a",c4)))))

sample_anno<-data.frame(unadj500=c1,
  adj500=c4,
  AIMS=tnbcClass$PAM50_AIMS
  #TNBC=tnbcClass$TNBCtype
  #umap=factor(clusters.umap[samples_use,"class"]),
  #hrd3=factor(clinAnno[samples_use,"HRD.3"])
  )
rownames(sample_anno)<-samples_use
sample_anno<-sample_anno[,ncol(sample_anno):1]

tiff(paste0(HOME,"/20191215_random500_heatmap_noAnno_pear_eucl_unadjClust_unadjBeta.tiff"),width=10*500,height=12*500,units="px",res=500,compression="lzw")
pheatmap(b2,cluster_rows = r1, cluster_cols = c3
  ,show_rownames=F,show_colnames=F
  ,main="",cutree_cols=2,fontsize=18
  ,annotation_col=sample_anno,annotation_colors=my_colour,annotation_legend=FALSE,annotation_names_col=F
  ,treeheight_row=0,treeheight_col=0,legend=F
)
dev.off()

##alluvial of results
c1<-cutree( hclust( as.dist( 1-cor(b1) ),method="ward.D"),2)
unique(c1[hclust( as.dist( 1-cor(b1) ),method="ward.D")$order])
#[1] 2 1 
c1<-cutree( hclust( as.dist( 1-cor(b2) ),method="ward.D"),2)
unique(c1[hclust( as.dist( 1-cor(b2) ),method="ward.D")$order])
#[1] 2 1

cl<-as.data.frame(cbind(unadjusted=cutree( hclust( as.dist( 1-cor(b2) ),method="ward.D"),2),
  adjusted=letters[1:2][cutree( hclust( as.dist( 1-cor(b1) ),method="ward.D"),2)]
  ),stringsAsFactors=F)

table(cl$unadjusted,cl$adjusted)
  #     a   b
  # 1 184   7
  # 2   8  36

cl<-as.data.frame(table(cl$unadjusted,cl$adjusted))
cl$Var2<-factor(cl$Var2,levels=c("b","a"))
cl$Var1<-factor(cl$Var1,levels=c("2","1"))

(pal2<-brewer.pal(5,"Set1")[c(2,1)])
(pal<-brewer.pal(5,"Set1")[c(2,1)])

pal<-rev(pal)
pal2<-rev(pal2)

q<-ggplot(cl,
       aes(y = Freq,
           axis1 = Var1, axis2 = Var2,)) +
  geom_alluvium(aes(fill = Var1 ),alpha=.75,       #fill = Var1
                width = 0, knot.pos = 1/5, reverse = TRUE
                ) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/20, reverse = TRUE,colour=c(pal,pal2),fill=c(pal,pal2)) +
  geom_text(stat = "stratum", infer.label = F, reverse = TRUE,size=10,fontface="bold",label=c(rev(c("b","a")),rev(c("b","a")))
    ) +
  scale_x_continuous(breaks = 1:2, labels = c("unadjusted", "adjusted")) +
  scale_fill_manual(values=rev(pal)) +
  #scale_fill_manual(values=rep("lightgrey",5)) +
  #coord_flip() +
  #restitle("") +
  theme_bw() +
  theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    )

pdf(paste0(HOME,"/20191203_random500_alluvial_pear_eucl_hClust_unadjBeta_to_adjBeta.pdf"),width=12,height=12,useDingbats=F)
q
dev.off()
rm(q)

a2<-image_read(paste0(HOME,"/20191215_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta.tiff"))
a1<-image_read(paste0(HOME,"/20191215_random500_heatmap_noAnno_pear_eucl_unadjClust_unadjBeta.tiff"))

a11<-image_read(paste0(HOME,"/20191215_random500_heatmap_noAnno_pear_eucl_adjClust_adjBeta_annotations.tiff"))

a11<-image_crop(a11,"1000x6000+4300")

tiff(paste0(HOME,"/20191215_random500_heatmap_noAnno_white.tiff"),width=.2*500,height=12*500,units="px",res=500,compression="lzw")
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
dev.off()
a9<-image_read(paste0(HOME,"/20191215_random500_heatmap_noAnno_white.tiff"))

out<-image_append(c(a9,
  a1,
  a9,
  a2,
  a9,
  a11
  ),stack = F)
out<-image_scale(out,"6000x")

image_write(out, path = paste0(HOME,"/20191215_random500_heatmap_noAnno_unadj_adj_combined.tiff"), format = "tiff")

rm(a1,a2,a9,a11,out)
rm(c1,c2,c3,c4,c5,r1,b,b1,b2)

gc()

rm(cl)

################################################################################
###Crop previous double image for new plot

a1<-image_read(paste0(HOME,"/20191203_top5k_heatmap_noAnno_unadj_adj_combined.tiff"))

a1
## A tibble: 1 x 7
#  format width height colorspace matte filesize density
#  <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
#1 TIFF    9000   6786 sRGB       FALSE 33047108 500x500

ww<- c(image_info(a1)$width , ceiling( image_info(a1)$height / 2 ))

a1<-image_crop(a1,paste0(ww,collapse="x"))

image_write(a1, path = paste0(HOME,"/20191215_random500_heatmap_noAnno_unadj_adj_combinedTopCropped.tiff"), format = "tiff")

rm(a1,ww)

gc()

################################################################################
###Do plot with all ATAC Promoter CpGs

##choose one promoter for doing figure with pre-post correction results

iii<-intersect(atacObjects$atacMethProbes,names(res))
ii<-match(iii,names(res) )
ij<-match(iii,atacObjects$atacMethProbes )

all(iii==names(res)[ii])
#[1] TRUE
all(iii==atacObjects$atacMethProbes[ij])
ij<-names(atacObjects$atacMethProbes)[ij]

#source("../function_doLmTests_modified_2.r")
#system.time( res2<-apply(testDat2,1,doLmTests) )

#system.time( res2<-apply(testDat[iii,],1,doLmTests) )
    
pdf(paste0(HOME,"/20191211_top5kBySd_allCpGsInAtacPromoters.pdf"),width=8,height=8,useDingbats=F)
par(mfrow=c(2,2),font=2,font.axis=2,font.lab=2,font.sub=2)

for(i in 1:length(iii)) {

  # plot(fracTum,betaData[iii[i],],col=res[[iii[i]]]$groups,pch=16,xlab="mean tumor fraction",ylab="raw beta",axes=F,
  # main=names(res)[ii][i],sub=ij[i]
  # )
  plot(fracTum,res[[iii[i]]]$y.orig,col=res[[iii[i]]]$groups,pch=16,xlab="mean tumor fraction",ylab="raw beta",axes=F,
  main=names(res)[ii][i],sub=ij[i]
  )
  axis(1,lwd=2,las=1,at=seq(-.4,1,by=.2),cex.axis=1.6)
  axis(2,lwd=2,las=1,cex.axis=1.6)
  plot(fracTum,res[[iii[i]]]$y.tum,col=res[[iii[i]]]$groups,pch=16,xlab="mean tumor fraction",ylab="adj beta",axes=F,
    main=iii[i]
  )
  axis(1,lwd=2,las=1,at=seq(-.4,1,by=.2),cex.axis=1.6)
  axis(2,lwd=2,las=1,cex.axis=1.6)

}
dev.off()
rm(i,ii,iii,ij)

################################################################################
###Pick one promoter as example..

str(tnbcClass)
#'data.frame':   235 obs. of  3 variables:
# $ TumorAssay: chr  "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
# $ TNBCtype  : chr  "BL2" "BL1" "M" "IM" ...
# $ PAM50_AIMS: chr  "Basal" "Basal" "Basal" "Basal" ...

table(res[["cg21237687"]]$groups,tnbcClass$TNBCtype)   
  #   BL1 BL2 IM LAR  M MSL NA UNS
  # 1  24   9 30  25 23  12  7  15
  # 2  22  13 16   5 18   2  3  11
   
table(res[["cg21237687"]]$groups,tnbcClass$PAM50_AIMS)   
  #  Basal Her2 LumB Normal
  # 1   104   23    1     17
  # 2    78    8    0      4

table(res[["cg21237687"]]$groups,clinAnno[samples_use,"HRD.3"])   
  #  [0.0,0.2) [0.2,0.7) [0.7,1.0]
  # 1        61         8        76
  # 2        22         5        63

probeN<-"cg21237687"
geneN<-"239|ALOX12"

pdf(paste0(HOME,"/20191216_atacPromoters_ALOX12_betaVsTumFrac_Concept.pdf"),width=16,height=16,useDingbats=F)
par(fig=c(0,.5,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=F)
##1
plot(1-fracTum,res[[probeN]]$y.orig,col=1,xlim=0:1,ylim=0:1,
  pch=16,cex=1.2,cex.lab=1.6,cex.main=2,
  main="ALOX12 promoter methylation",
  xlab="1-tumor fraction",
  ylab="unadjusted beta",
  axes=F
)
points(1-fracTum,res[[probeN]]$y.orig,col=res[[probeN]]$groups,pch=16,cex=1.1)
abline(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==1]~(1-fracTum)[res[[probeN]]$groups==1]),col=1,lwd=3)
abline(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==2]~(1-fracTum)[res[[probeN]]$groups==2]),col=2,lwd=3)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)
axis(2,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)
legend("topright",legend=c("pop1","pop2"),col=1:2,pch=16,bty="n",cex=1.6)
text(.2,.25,
  paste0("pop 1 intercept=",
  round(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==1]~(1-fracTum)[res[[probeN]]$groups==1])$coeff[1],2)
  ),cex=1.6)
text(.2,.65,
  paste0("pop 2 intercept=",
  round(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==2]~(1-fracTum)[res[[probeN]]$groups==2])$coeff[1],2)
  ),cex=1.6)

##2
par(fig=c(.5,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(1-fracTum,res[[probeN]]$y.tum,col=1,xlim=0:1,ylim=0:1,
  pch=16,cex=1.2,cex.lab=1.6,cex.main=2,
  main="Adjusted ALOX12 promoter methylation",
  xlab="1-tumor fraction",
  ylab="adjusted beta",
  axes=F
)
points(1-fracTum,res[[probeN]]$y.tum,col=res[[probeN]]$groups,pch=16,cex=1.1)
abline(lm(res[[probeN]]$y.tum[res[[probeN]]$groups==1]~(1-fracTum)[res[[probeN]]$groups==1]),col=1,lwd=3)
abline(lm(res[[probeN]]$y.tum[res[[probeN]]$groups==2]~(1-fracTum)[res[[probeN]]$groups==2]),col=2,lwd=3)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)
axis(2,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)

##3
par(fig=c(0,.5,0,.5),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(fracTum,res[[probeN]]$y.orig,col=1,xlim=0:1,ylim=0:1,
  pch=16,cex=1.2,cex.lab=1.6,cex.main=2,
  main="ALOX12 promoter methylation",
  xlab="tumor fraction",
  ylab="unadjusted beta",
  axes=F
)
points(fracTum,res[[probeN]]$y.orig,col=res[[probeN]]$groups,pch=16,cex=1.1)
abline(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==1]~(fracTum)[res[[probeN]]$groups==1]),col=1,lwd=3)
abline(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==2]~(fracTum)[res[[probeN]]$groups==2]),col=2,lwd=3)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)
axis(2,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)
legend("topright",legend=c("pop1","pop2"),col=1:2,pch=16,bty="n",cex=1.6)
text(.8,.3,
  paste0("pop 1 intercept=",
  round(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==1]~(fracTum)[res[[probeN]]$groups==1])$coeff[1],2)
  ),cex=1.6)
text(.8,.65,
  paste0("pop 2 intercept=",
  round(lm(res[[probeN]]$y.orig[res[[probeN]]$groups==2]~(fracTum)[res[[probeN]]$groups==2])$coeff[1],2)
  ),cex=1.6)

##4
par(fig=c(.5,1,0,.5),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(fracTum,res[[probeN]]$y.norm,col=1,xlim=0:1,ylim=0:1,
  pch=16,cex=1.2,cex.lab=1.6,cex.main=2,
  main="ALOX12 inferred normal methylation",
  xlab="tumor fraction",
  ylab="inferred normal beta",
  axes=F
)
points(fracTum,res[[probeN]]$y.norm,col=res[[probeN]]$groups,pch=16,cex=1.1)
abline(lm(res[[probeN]]$y.norm[res[[probeN]]$groups==1]~(fracTum)[res[[probeN]]$groups==1]),col=1,lwd=3)
abline(lm(res[[probeN]]$y.norm[res[[probeN]]$groups==2]~(fracTum)[res[[probeN]]$groups==2]),col=2,lwd=3)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)
axis(2,lwd=2,las=1,at=seq(0,1,.2),cex.axis=1.6)

dev.off()

rm(geneN,probeN)

################################################################################
###Do brca1 promoter

table(res[["cg09441966"]]$groups,clinAnno[samples_use,"BRCA1_PromMetPc_Class"])
  #    0   1
  # 1   5  57
  # 2 173   0

#      0   1
#  1 177   2
#  2   1  55

table(res[["cg09441966"]]$y.orig>.3,clinAnno[samples_use,"BRCA1_PromMetPc_Class"])
  #        0   1
  # FALSE 178   7
  # TRUE    0  50

i<-"cg09441966"

pdf(paste0(HOME,"/20191216_atacPromoters_BRCA1_betaVsTumFrac_adjNonAdj.pdf"),width=12,height=12,useDingbats=F)
par(fig=c(0,.5,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=F)
##1
bclass<-as.integer(factor(paste0(res[[i]]$groups,clinAnno[samples_use,"BRCA1_PromMetPc_Class"])))
plot(fracTum,res[[i]]$y.orig,col=1,xlim=0:1,ylim=0:1,
  pch=c(17,16,16)[bclass],cex=1.2,
  main="BRCA1 methylation vs tumor fraction",
  xlab="tumor fraction",
  ylab="unadjusted beta",
  axes=F
)
points(fracTum,res[[i]]$y.orig,col=res[[i]]$groups,pch=c(17,16,16)[bclass],cex=1.1)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex=1.5)
axis(2,lwd=2,las=1,at=seq(0,1,.2),cex=1.5)
legend("topleft",legend=c("pyro concordant","pyro discordant"),col=c(1),pch=c(16,17),bty="n",cex=1.6)
abline(lm(res[[i]]$y.orig[res[[i]]$groups==1]~fracTum[res[[i]]$groups==1]),col=1,lwd=3)
abline(lm(res[[i]]$y.orig[res[[i]]$groups==2]~fracTum[res[[i]]$groups==2]),col=2,lwd=3)
#abline(lm(res[[i]]$y.orig[res[[i]]$groups==3]~fracTum[res[[i]]$groups==3]),col=3,lwd=3)

##2
par(fig=c(.5,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(fracTum,res[[i]]$y.tum,col=1,xlim=0:1,ylim=0:1,
  pch=16,cex=1.2,
  main="BRCA1 methylation vs tumor fraction",
  xlab="tumor fraction",
  ylab="adjusted beta",
  axes=F
)
points(fracTum,res[[i]]$y.tum,col=res[[i]]$groups,pch=16,cex=1.1)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex=1.6)
axis(2,lwd=2,las=1,at=seq(0,1,.2),cex=1.6)
abline(lm(res[[i]]$y.tum[res[[i]]$groups==1]~fracTum[res[[i]]$groups==1]),col=1,lwd=3)
abline(lm(res[[i]]$y.tum[res[[i]]$groups==2]~fracTum[res[[i]]$groups==2]),col=2,lwd=3)
#abline(lm(res[[i]]$y.tum[res[[i]]$groups==3]~fracTum[res[[i]]$groups==3]),col=3,lwd=3)

##3
par(fig=c(0,.5,0,.5),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
bclass<-as.integer(factor(paste0(res[[i]]$groups,clinAnno[samples_use,"BRCA1_PromMetPc_Class"])))
length(res[[i]]$groups)
#[1] 235

table(res[[i]]$groups,clinAnno[samples_use,"BRCA1_PromMetPc_Class"])
  #   0   1
  # 1   5  57
  # 2 173   0

table(paste0(res[[i]]$groups,clinAnno[samples_use,"BRCA1_PromMetPc_Class"]))
# 10  11  20 
#   5  57 173 

plot(fracTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=1,xlim=0:1,ylim=c(0,100),
  pch=c(17,16,16)[bclass],cex=1.2,
  main="Pyro BRCA1 methylation vs tumor fraction",
  xlab="tumor fraction",
  ylab="Pyro BRCA1 methylation percent",
  axes=F
)
#abline(lm(clinAnno[samples_use,"BRCA1_PromMetPc"][res[[i]]$groups==1]~fracTum[res[[i]]$groups==1]),col=1,lwd=3)
#abline(lm(clinAnno[samples_use,"BRCA1_PromMetPc"][res[[i]]$groups==2]~fracTum[res[[i]]$groups==2]),col=2,lwd=3)
#abline(lm(clinAnno[samples_use,"BRCA1_PromMetPc"][res[[i]]$groups==3]~fracTum[res[[i]]$groups==3]),col=3,lwd=3)
points(fracTum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=res[[i]]$groups,pch=c(17,16,16)[bclass],cex=1.1)
points(fracTum[bclass==2],clinAnno[samples_use,"BRCA1_PromMetPc"][bclass==2],col=res[[i]]$groups[bclass==2],pch=c(16),cex=1.1)
points(fracTum[bclass==3],clinAnno[samples_use,"BRCA1_PromMetPc"][bclass==3],col=res[[i]]$groups[bclass==3],pch=c(16),cex=1.1)
points(fracTum[bclass==1],clinAnno[samples_use,"BRCA1_PromMetPc"][bclass==1],col=res[[i]]$groups[bclass==1],pch=c(17),cex=1.1)
#points(fracTum[bclass==4],clinAnno[samples_use,"BRCA1_PromMetPc"][bclass==4],col=3,pch=c(16),cex=1.1)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex=1.6)
axis(2,lwd=2,las=1,at=seq(0,100,20),cex=1.6)
legend("topleft",legend=c("pyro meth","pyro unmeth","discordant"),col=c(1,2,1),pch=c(16,16,17),bty="n",cex=1.5)
abline(h=7,lwd=3,lty=2,col="lightgrey")

##4
par(fig=c(.5,1,0,.5),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(res[[i]]$y.tum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=1,xlim=0:1,ylim=c(0,100),
  pch=c(17,16,16)[bclass],cex=1.2,
  main="Pyro BRCA1 methylation vs adjusted beta",
  xlab="adjusted beta",
  ylab="Pyro BRCA1 methylation percent",
  axes=F
)
points(res[[i]]$y.tum,clinAnno[samples_use,"BRCA1_PromMetPc"],col=res[[i]]$groups,pch=c(17,16,16)[bclass],cex=1.1)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex=1.6)
axis(2,lwd=2,las=1,at=seq(0,100,20),cex=1.6)
legend("topleft",legend=c("concordant (N=230)","discordant (N=5)"),col=1,pch=c(16,17),bty="n",cex=1.5)
rm(bclass)

dev.off()
rm(i)

################################################################################
###Correlation to GEX pre/post correction

str(gexAnno)
#'data.frame':   18776 obs. of  6 variables:
# $ Gene.ID      : chr  "ENSG00000000003.14" "ENSG00000000005.5" "ENSG00000000419.12" "ENSG00000000457.13" ...
# $ Gene.Name    : chr  "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
# $ HGNC         : chr  "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
# $ EntrezGene   : chr  "7105" "64102" "8813" "57147" ...
# $ RefSeq       : chr  "NM_003270.3;NM_001278740.1;NM_001278742.1" "NM_022144.2" "NM_001317035.1;NM_001317034.1;NM_001317036.1;NA" "NA;NM_020423.6;NM_181093.3" ...
# $ transcript_id: chr  "ENST00000373020.8;ENST00000612152.4;ENST00000614008.4" "ENST00000373031.4" "ENST00000371582.8;ENST00000371584.8;ENST00000371588.9;ENST00000413082.1" "ENST00000367770.5;ENST00000367771.10;ENST00000367772.8;ENST00000423670.1" ...

str(gexTpm)
# num [1:18776, 1:236] 3.658 0.343 6.045 1.985 3.137 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:18776] "7105|TSPAN6" "64102|TNMD" "8813|DPM1" "57147|SCYL3" ...
#  ..$ : chr [1:236] "PD35926a" "PD35927a" "PD35928a" "PD35929a" ...

#geneCoords-object has coordinates

length(res)
#[1] 5000

p5000<-probeAnno[names(res)]
 
all(names(p5000)==names(res))
#[1] TRUE

pCoord<-promoters(geneCoords,upstream=1500,downstream=500)

ol5000<-findOverlaps(pCoord,p5000)

str(ol5000)
#Formal class 'SortedByQueryHits' [package "S4Vectors"] with 6 slots
#  ..@ from           : int [1:1455] 30 42 42 49 101 101 101 125 131 202 ...
#  ..@ to             : int [1:1455] 2530 2611 2612 1274 4121 4122 4123 2425 2411 11 ...
#  ..@ nLnode         : int 18776
#  ..@ nRnode         : int 5000
#  ..@ elementMetadata: NULL
#  ..@ metadata       : list()

table(unlist(lapply(split(subjectHits(ol5000),queryHits(ol5000)),length)))
#  1   2   3   4   5   6   7   8   9  11
#485 160  78  39  23  10   2   3   4   1

length(unlist(lapply(split(subjectHits(ol5000),queryHits(ol5000)),length)))
#[1] 805

allSD<-apply(do.call("rbind",lapply(res,function(x) x$y.orig)),1,sd)
allSD<-allSD[subjectHits(ol5000)]

allSD<-unlist(lapply(split(allSD,queryHits(ol5000)),function(x) names(x)[which.max(x)]))
names(allSD)<-names(geneCoords)[as.integer(names(allSD))]

head(allSD)
#55365|TMEM176A     9108|MTMR7       952|CD38     30812|SOX8      6863|TAC1
#  "cg02244695"   "cg12296772"   "cg26043257"   "cg05520409"   "cg09236284"
#     1750|DLX6
#  "cg04599026"

allCorrs<-matrix(nrow=length(allSD),ncol=3,dimnames=list(names(allSD),c("corr.raw","corr.adj","corrMePrePost")))

d1<-do.call("rbind",lapply(res,function(x) x$y.orig))
d2<-do.call("rbind",lapply(res,function(x) x$y.tum))

for(i in rownames(allCorrs)) {

      allCorrs[i,1]<-cor(d1[allSD[i],],gexTpm[i,samples_use])
      allCorrs[i,2]<-cor(d2[allSD[i],],gexTpm[i,samples_use])
      allCorrs[i,3]<-cor(d1[allSD[i],],d2[allSD[i],])

} ; rm(i,d1,d2)

##do plots of effects of correction
pdf(paste0(HOME,"/20191216_gexCorr_top5000_adjNonAdj.pdf"),width=10,height=10,useDingbats=F)

par(fig=c(0,.5,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=F)
plot(density(allCorrs[,3]),col=1,xlim=0:1,
  pch=16,cex=1.2,lwd=3,
  main="Beta correlation pre/post correction",
  xlab="correlation",
  ylab="Density",
  axes=F
)
axis(1,lwd=2,las=1,at=seq(0,1,.2),cex=1.2)
axis(2,lwd=2,las=1,cex=1.2)
text(0,14,"N=805 CpG-gene pairs",pos=4,cex=1.2)

par(fig=c(.5,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(allCorrs[,1],allCorrs[,2],col=1,xlim=c(-1,1),ylim=c(-1,1),
  pch=16,cex=1.2,lwd=3,
  main="Correlation Illumina beta to gene expression",
  xlab="unadjusted beta",
  ylab="adjusted beta",
  axes=F
)
axis(1,lwd=2,las=1,at=seq(-1,1,.5),cex=1.2)
axis(2,lwd=2,las=1,at=seq(-1,1,.5),cex=1.2)
text(-1,.75,paste0("mean abs difference = ",round(mean(abs(diff(t(allCorrs[,1:2])))),3)),pos=4,cex=1.2)

##add BRCA1 correction
i<-"cg09441966"
par(fig=c(0,.5,0,.5),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(res[[i]]$y.orig,gexTpm["672|BRCA1",samples_use],col=res[[i]]$groups,xlim=c(0,1),ylim=c(0,5),
  pch=16,cex=1.2,lwd=3,
  main="BRCA1 beta-gex correlation",
  xlab="unadjusted beta",
  ylab="log2(TPM+1) BRCA1",
  axes=F
)
legend("topright",legend=c("pop1","pop2"),col=1:2,bty="n",pch=16)
abline(lm(gexTpm["672|BRCA1",samples_use]~res[[i]]$y.orig),lwd=3,lty=2)
axis(1,lwd=2,las=1,at=seq(0,1,.25),cex=1.2)
axis(2,lwd=2,las=1,cex=1.2)
text(.25,4.5,paste0("adj R2 = ",round(summary(lm(gexTpm["672|BRCA1",samples_use]~res[[i]]$y.orig))$adj.r.squared,3)),pos=4,cex=1.2)

par(fig=c(.5,1,0,.5),font=2,font.axis=2,font.lab=2,font.sub=2,new=T)
plot(res[[i]]$y.tum,gexTpm["672|BRCA1",samples_use],col=res[[i]]$groups,xlim=c(0,1),ylim=c(0,5),
  pch=16,cex=1.2,lwd=3,
  main="BRCA1 beta-gex correlation",
  xlab="adjusted beta",
  ylab="log2(TPM+1) BRCA1",
  axes=F
)
legend("topright",legend=c("pop1","pop2"),col=1:2,bty="n",pch=16)
abline(lm(gexTpm["672|BRCA1",samples_use]~res[[i]]$y.tum),lwd=3,lty=2)
axis(1,lwd=2,las=1,at=seq(0,1,.25),cex=1.2)
axis(2,lwd=2,las=1,cex=1.2)
text(.25,4.5,paste0("adj R2 = ",round(summary(lm(gexTpm["672|BRCA1",samples_use]~res[[i]]$y.tum))$adj.r.squared,3)),pos=4,cex=1.2)
rm(i)
dev.off()

rm(p5000,pCoord,ol5000)

rownames(allCorrs)<-paste(rownames(allCorrs),allSD,sep="|")
rm(allSD)

save(allCorrs,file=paste0(HOME,"/20200116_object_gexCorr_top5000_adjNonAdj.RData"))

################################################################################
###redo plot for rand500 iterations

##track sens+spec+acc for all iter

load(file=paste0(HOME,"/20191214_top100percentBySd_basalVsLuminalSplitIn100randomSets_UsedProbeSets.RData"))

##do multicore
no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

clusterEvalQ(cl, {
  library("flexmix")
})

clusterSetRNGStream(cl, 20200918) ##will not make exactly replicable..

resMat2<-matrix(nrow=length(p_list),ncol=9)
colnames(resMat2)<-c("rawAcc","rawSens","rawSpec",
  "adjAcc","adjSens","adjSpec",
  "binAcc","binSens","binSpec")

refStat<-factor(1+(tnbcClass$PAM50_AIMS != "Basal"))
for( i in 1:length(p_list)) {
     cat(".")
     if(i%%10==0)  cat(" ",i,"\n")
     
    b<-parRapply(cl = cl, betaData[p_list[[i]],samples_use], adjustBeta,purity=fracTum2,snames=samples_use)

     ##unadjusted
     b1<-do.call("rbind",lapply(b,function(x) x$y.orig))
     b1<-b1[!apply(b1,1,function(x) any(is.na(x))),]
     b1<-factor(cutree( hclust( as.dist(1-cor(b1)),"ward.D"), k=2))
     r1<-confusionMatrix(b1,reference=refStat)
     r1<-c(r1$overall[1],r1$byClass[1:2])
     r1<-rbind(r1,1-r1)
     r1<-r1[which.max(r1[,1]),]
     resMat2[i,1:3]<-as.numeric(r1)
     ##adjusted
     b2<-do.call("rbind",lapply(b,function(x) x$y.tum))
     b2<-b2[!apply(b2,1,function(x) any(is.na(x))),]
     b2<-factor(cutree( hclust( as.dist(1-cor(b2)),"ward.D"), k=2))
     r1<-confusionMatrix(b2,reference=refStat)
     r1<-c(r1$overall[1],r1$byClass[1:2])
     r1<-rbind(r1,1-r1)
     r1<-r1[which.max(r1[,1]),]
     resMat2[i,4:6]<-as.numeric(r1)
     ##dichotomized
     b3<-do.call("rbind",lapply(b,function(x) x$y.orig > .3))
     b3<-b3[!apply(b3,1,function(x) any(is.na(x))),]
     b3<-factor(cutree( hclust( as.dist(1-cor(b3)),"ward.D"), k=2))
     r1<-confusionMatrix(b3,reference=refStat)
     r1<-c(r1$overall[1],r1$byClass[1:2])
     r1<-rbind(r1,1-r1)
     r1<-r1[which.max(r1[,1]),]
     resMat2[i,7:9]<-as.numeric(r1)
}     
     
rm(i,p_list,b,b1,b2,b3,r1,refStat,cl)

pdf(paste0(HOME,"/20200121_top100percentBySd_confusionStats_basalVsLuminalSplitIn100randomSets.pdf"),width=10,height=10,useDingbats=F)
par(fig=c(0,.5,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=F)

boxplot(resMat2[,4]-resMat2[,1],at=1,xlim=c(0.5,8.5),ylim=c(-1.,1.1),width=2,axes=F,lwd=2,
  main="Discrimination of PAM50 Basal vs Luminal split"
)
abline(h=0,lwd=3,col="lightgrey",lty=2)
boxplot(resMat2[,7]-resMat2[,1],at=2,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,5]-resMat2[,2],at=4,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,8]-resMat2[,2],at=5,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,6]-resMat2[,3],at=7,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,9]-resMat2[,3],at=8,add=T,width=2,axes=F,lwd=2)

axis(1,at=c(1:2,4:5,7:8),lwd=2,las=2,cex=1.2,
  labels=c("Adjusted",
  "Beta>0.3",
  "Adjusted",
  "Beta>0.3",
  "Adjusted",
  "Beta>0.3")
)
axis(2,at=round(seq(-1,1,length.out=5),2),lwd=2,las=1,cex=1.2)
mtext(side=2, "Relative to unadjusted data",font=2,line=2.5,cex=1.2)
lines(x=c(1,2),y=c(.85,.85),lwd=3)
text(1.5,.85,labels="Accuracy",pos=3)

lines(x=c(4,5),y=c(.85,.85),lwd=3)
text(4.5,.85,labels="Sensitivity",pos=3)

lines(x=c(7,8),y=c(.85,.85),lwd=3)
text(7.5,.85,labels="Specificity",pos=3)

##absolute terms
par(fig=c(.5,1,.5,1),font=2,font.axis=2,font.lab=2,font.sub=2,cex.lab=1.2,cex.lab=1.2,new=T)

boxplot(resMat2[,4],at=1,xlim=c(0.5,8.5),ylim=c(0,1),width=2,axes=F,lwd=2,
  main="Discrimination of PAM50 Basal vs Luminal split"
)
#abline(h=0,lwd=3,col="lightgrey",lty=2)
boxplot(resMat2[,7],at=2,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,5],at=4,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,8],at=5,add=T,width=2,axes=F,lwd=2)

boxplot(resMat2[,6],at=7,add=T,width=2,axes=F,lwd=2)
boxplot(resMat2[,9],at=8,add=T,width=2,axes=F,lwd=2)

axis(1,at=c(1:2,4:5,7:8),lwd=2,las=2,cex=1.2,
  labels=c("Adjusted",
  "Beta>0.3",
  "Adjusted",
  "Beta>0.3",
  "Adjusted",
  "Beta>0.3")
)
axis(2,at=seq(0,1,length.out=5),lwd=2,las=1,cex=1.2)
mtext(side=2, "Absolute level",font=2,line=2.5,cex=1.2)
lines(x=c(1,2),y=c(.05,.05),lwd=3)
text(1.5,.05,labels="Accuracy",pos=3)

lines(x=c(4,5),y=c(.05,.05),lwd=3)
text(4.5,.05,labels="Sensitivity",pos=3)

lines(x=c(7,8),y=c(.05,.05),lwd=3)
text(7.5,.05,labels="Specificity",pos=3)

dev.off()

##deltas - adj.
t.test(resMat2[,4]-resMat2[,1])
#         One Sample t-test

# data:  resMat2[, 4] - resMat2[, 1]
# t = 29.064, df = 99, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.2516462 0.2885240
# sample estimates:
# mean of x 
# 0.2700851 

t.test(resMat2[,5]-resMat2[,2])
#         One Sample t-test

# data:  resMat2[, 5] - resMat2[, 2]
# t = 19.755, df = 99, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.2648265 0.3239647
# sample estimates:
# mean of x 
# 0.2943956 

t.test(resMat2[,6]-resMat2[,3])
#         One Sample t-test

# data:  resMat2[, 6] - resMat2[, 3]
# t = 6.5142, df = 99, p-value = 3.065e-09
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.1297641 0.2434434
# sample estimates:
# mean of x 
# 0.1866038 


##deltas - beta>.3
t.test(resMat2[,7]-resMat2[,1])
#         One Sample t-test

# data:  resMat2[, 7] - resMat2[, 1]
# t = -4.3022, df = 99, p-value = 3.974e-05
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  -0.07094633 -0.02616005
# sample estimates:
#   mean of x 
# -0.04855319 

t.test(resMat2[,8]-resMat2[,2])
#         One Sample t-test

# data:  resMat2[, 8] - resMat2[, 2]
# t = -2.1794, df = 99, p-value = 0.03168
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  -0.070749732 -0.003316202
# sample estimates:
#   mean of x 
# -0.03703297 

t.test(resMat2[,9]-resMat2[,3])
#         One Sample t-test

# data:  resMat2[, 9] - resMat2[, 3]
# t = -2.2099, df = 99, p-value = 0.02942
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  -0.167227780 -0.008998635
# sample estimates:
#   mean of x 
# -0.08811321 

save(resMat2,file=paste0(HOME,"/20200121_top50percentBySd_basalVsLuminalSplitIn100randomSets_resultsConfusionMatrix.RData"))

################################################################################
##Temp save

#HOME<-"I:/data/adjustBetas"
#GIT<-"F:/gitProjects/adjustBetas"
save.image(paste0(HOME,"/20200121_tempSave_calibrateBetasPaper.RData"))

#load(paste0(HOME,"/20200121_tempSave_calibrateBetasPaper.RData"))

################################################################################
##Show additional examples of beta correciton for fig1





################################################################################
##Show wether CN-status influences presence of 3 groups



################################################################################
##Show how IM-group dissolves in corrected data













###END