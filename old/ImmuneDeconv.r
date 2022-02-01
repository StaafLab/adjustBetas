##home
HOME<-"I:/data/tcgaBrca"
MANIFEST<-"F:/gitProjects/getTcga/manifest"



load(file=paste0(HOME,"/","coreData_gexCnWesMeAtac_unique_matched_samples.RData"))

ls()
#  [1] "data.atac.core"      "data.cn.core"        "data.counts.core"   
#  [4] "data.fpkm.core"      "data.mutations.core" "data.seg.core"      
#  [7] "data.uq.core"        "HOME"                "MANIFEST"           
# [10] "mutations.core"      "samples.me.core"     "TUMOR_TYPE"        

##load betaFinal,sampleMat and probesKeep
load(file=paste0(HOME,"/me/","workspace_minfiNormInfiniumAdjFinalBetas.RData"))

ls()
#  [1] "betaFinal"           "data.atac.core"      "data.cn.core"       
#  [4] "data.counts.core"    "data.fpkm.core"      "data.mutations.core"
#  [7] "data.seg.core"       "data.uq.core"        "HOME"               
# [10] "MANIFEST"            "mutations.core"      "probesKeep"         
# [13] "sampleMat"           "samples.me.core"     "TUMOR_TYPE"         


##from TCGA pancan companion page
sample.purity<-read.table(paste0(MANIFEST,"/pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"),sep="\t",header=TRUE,as.is=TRUE)

str(sample.purity)
# 'data.frame':   10786 obs. of  10 variables:
#  $ array                    : chr  "TCGA-OR-A5J1-01" "TCGA-OR-A5J2-01" "TCGA-OR-A5J3-01" "TCGA-OR-A5J4-01" ...
#  $ sample                   : chr  "TCGA-OR-A5J1-01A-11D-A29H-01" "TCGA-OR-A5J2-01A-11D-A29H-01" "TCGA-OR-A5J3-01A-11D-A29H-01" "TCGA-OR-A5J4-01A-11D-A29H-01" ...
#  $ call.status              : chr  "called" "called" "called" "called" ...
#  $ purity                   : num  0.9 0.89 0.93 0.87 0.93 0.69 0.84 0.76 0.84 0.75 ...
#  $ ploidy                   : num  2 1.3 1.27 2.6 2.79 3.34 2.6 1.23 2.61 5.52 ...
#  $ Genome.doublings         : num  0 0 0 1 1 1 1 0 1 2 ...
#  $ Coverage.for.80..power   : num  9 6 5 12 12 17 12 7 12 32 ...
#  $ Cancer.DNA.fraction      : num  0.9 0.84 0.89 0.89 0.95 0.79 0.87 0.66 0.87 0.89 ...
#  $ Subclonal.genome.fraction: num  0.02 0.16 0.11 0.08 0.15 0.06 0 0.52 0.23 0.12 ...
#  $ solution                 : chr  "new" "new" "new" "new" ...

##use cancer DNA fraction as this is what should correspond best to what is seen on betas
str(samples.me.core)
# 'data.frame':   669 obs. of  2 variables:
#  $ tcga_id: chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ array  : chr  "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...

length(unique(sample.purity$array)) == length((sample.purity$array))
#[1] TRUE

length(unique(samples.me.core$tcga_id)) == length(unique(sub(".$","",samples.me.core$tcga_id)))
#[1] TRUE

length(intersect(sample.purity$array,sub(".$","",samples.me.core$tcga_id)))
#[1] 663

uids<-intersect(sample.purity$array,sub(".$","",samples.me.core$tcga_id))

all(sample.purity$array[match(uids,sample.purity$array)]==uids)
#[1] TRUE

sample.purity<-sample.purity[match(uids,sample.purity$array),]

##get ascat tumor
fracA<-sample.purity$purity

##get battenberg tumor
fracB<-sample.purity$Cancer.DNA.fraction



https://github.com/darneson/MethylResolver
install.packages('devtools')
devtools::install_github(repo = 'darneson/MethylResolver')
library(MethylResolver)

?MethylResolver
MethylSig2<-MethylSig
MethylSig2<-MethylSig[intersect(rownames(MethylSig),rownames(betaOrig)),]

str(MethylSig2)
# 'data.frame':   386 obs. of  11 variables:
#  $ Mon      : num  0.708 0.73 0.719 0.752 0.838 ...
#  $ Dendritic: num  0.566 0.606 0.578 0.643 0.661 ...
#  $ Macro    : num  0.572 0.57 0.599 0.629 0.663 ...
#  $ Neu      : num  0.342 0.609 0.363 0.479 0.729 ...
#  $ Eos      : num  0.353 0.36 0.235 0.384 0.555 ...
#  $ Treg     : num  0.191 0.173 0.125 0.167 0.16 ...
#  $ Tnaive   : num  0.1084 0.1333 0.0858 0.1256 0.1448 ...
#  $ Tmem     : num  0.241 0.183 0.147 0.171 0.187 ...
#  $ CD8      : num  0.123 0.145 0.0539 0.1363 0.1469 ...
#  $ NK       : num  0.206 0.22 0.157 0.277 0.278 ...
#  $ Bcell    : num  0.1236 0.1581 0.0967 0.2941 0.2692 ...

a1<-MethylResolver(betaOrig[,1:50],betaPrime=F,methylSig = MethylSig2)

a1<-read.delim("MethylResolver.txt")


#### EpiDISH

https://www.bioconductor.org/packages/release/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html
BiocManager::install("EpiDISH")
library("EpiDISH")


load(rgSetFnorm,file=paste0(HOME,"/me/","object_minfi_rgSetFnorm.RData"))
betaData <- getBeta(rgSetFnorm)

str( betaData )


centEpiFibIC.m2<-centEpiFibIC.m[intersect(rownames(centEpiFibIC.m),rownames(betaData)),]
betaData.m2<-betaData[intersect(rownames(centEpiFibIC.m),rownames(betaData)),]

out.l <- epidish(beta.m = betaData.m2, ref.m = centEpiFibIC.m2, method = "RPC") 




##locked new data
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167998


##Salas
https://github.com/immunomethylomics/FlowSorted.Blood.EPIC
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554


##Reinius
https://bioconductor.org/packages/release/data/experiment/html/FlowSorted.Blood.450k.html
FlowSorted.Blood.450k
