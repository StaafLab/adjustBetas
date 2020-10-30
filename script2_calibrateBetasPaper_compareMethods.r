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

##get Shamik's cell type deconvolution results
load(paste0(HOME,"/data/Cell_type_deconvolution/","Cell_type_fracs.RData"))

str(Cell_type_fractions)
# List of 4
#  $ EpiDish_centBloodSub: num [1:9, 1:235] 73.8 0 7.22 5.81 3.5 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:9] "Epi" "Fib" "B" "NK" ...
#   .. ..$ : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ EpiDish_centBloodDHS: num [1:9, 1:235] 73.8 0 7.15 0 3.73 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:9] "Epi" "Fib" "B" "NK" ...
#   .. ..$ : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ xCell               : num [1:64, 1:235] 0 3.101 0 0.165 0 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:64] "Adipocytes" "Astrocytes" "B-cells" "Basophils" ...
#   .. ..$ : chr [1:235] "PD31028a" "PD31029a" "PD31030a" "PD31031a" ...
#  $ CibersortX          :'data.frame':   235 obs. of  9 variables:
#   ..$ epithelial : num [1:235] 0.488 0.419 0.405 0.282 0.537 ...
#   ..$ macrophage : num [1:235] 0.1423 0.0904 0.1 0.2212 0.0416 ...
#   ..$ stroma     : num [1:235] 0.239 0.161 0.322 0.123 0.225 ...
#   ..$ Bcell      : num [1:235] 0.0359 0.1864 0.0554 0.2382 0.0629 ...
#   ..$ endothelial: num [1:235] 0.0915 0.1366 0.118 0.1234 0.1333 ...
#   ..$ Tcell      : num [1:235] 0.00422 0.00682 0 0.01198 0 ...
#   ..$ P-value    : num [1:235] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ Correlation: num [1:235] 0.636 0.815 0.675 0.721 0.63 ...
#   ..$ RMSE       : num [1:235] 0.79 0.706 0.759 0.716 0.806 ...

################################################################################
##check different methods

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
names(fracA)<-sampleSets$samplesAll

##get battenberg tumor
fracB<-clinAnno[sampleSets$samplesAll,"BATTENBERG_TUMOUR_FRAC"]
names(fracB)<-sampleSets$samplesAll

colnames(clinAnno)[grep("PLOIDY",colnames(clinAnno))]
#[1] "ASCAT_PLOIDY"      "BATTENBERG_PLOIDY"

##get ascat ploidy
ploiA<-clinAnno[sampleSets$samplesAll,"ASCAT_PLOIDY"]
names(ploiA)<-sampleSets$samplesAll

##get battenberg ploidy
ploiB<-clinAnno[sampleSets$samplesAll,"BATTENBERG_PLOIDY"]
names(ploiB)<-sampleSets$samplesAll


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
#rm(fracA,fracB)
names(fracTum)<-sampleSets$samplesAll

##get Shamik immune deconvolution data
shamik1<-as.data.frame(t(Cell_type_fractions[[1]]),stringsAsFactors=FALSE)

shamik2<-as.data.frame(t(Cell_type_fractions[[2]]),stringsAsFactors=FALSE)

shamik3<-as.data.frame(t(Cell_type_fractions[[3]]),stringsAsFactors=FALSE)

shamik4<-as.data.frame(Cell_type_fractions[[4]],stringsAsFactors=FALSE)

all.equal(rownames(shamik1),rownames(shamik2))
#[1] TRUE
all.equal(rownames(shamik1),rownames(shamik3))
#[1] TRUE
all.equal(rownames(shamik1),rownames(shamik4))
#[1] TRUE

all.equal(rownames(shamik1),sampleSets$samplesAll)
#[1] TRUE

str(shamik1)
# 'data.frame':   235 obs. of  9 variables:
#  $ Epi   : num  73.8 35.5 74.5 36.8 57.8 ...
#  $ Fib   : num  0 21.5 16.68 9.96 17.78 ...
#  $ B     : num  7.22 6.61 1.25 5.81 5 ...
#  $ NK    : num  5.81 12.54 2.58 14.07 4.27 ...
#  $ CD4T  : num  3.5 12.01 1.33 10.2 3.95 ...
#  $ CD8T  : num  0 0 0 1.75 0 ...
#  $ Mono  : num  5.31 7.68 2.7 9.79 11.16 ...
#  $ Neutro: num  3.986 4.127 0.971 11.579 0 ...
#  $ Eosino: num  0.369 0 0 0 0 ...

str(shamik2)
# 'data.frame':   235 obs. of  9 variables:
#  $ Epi   : num  73.8 35.5 74.5 36.8 57.8 ...
#  $ Fib   : num  0 21.5 16.68 9.96 17.78 ...
#  $ B     : num  7.15 1.55 0 0 2.58 ...
#  $ NK    : num  0 9.95 1.488 11.84 0.713 ...
#  $ CD4T  : num  3.73 13.11 2.59 6.32 10.04 ...
#  $ CD8T  : num  0 0 0 4.05 0 ...
#  $ Mono  : num  15.31 15.93 4.76 27.28 11.04 ...
#  $ Neutro: num  0 2.43 0 3.72 0 ...
#  $ Eosino: num  0 0 0 0 0 ...

str(shamik3)
# 'data.frame':   235 obs. of  64 variables:
#  $ Adipocytes                   : num  0 0.929 1.023 0 0 ...
#  $ Astrocytes                   : num  3.1 0 5.3 0 0 ...
#  $ B-cells                      : num  0 5.88 0 6.66 0 ...
#  $ Basophils                    : num  0.165 2.385 4.061 12.679 0 ...
#  $ CD4+ T-cells                 : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ CD4+ Tcm                     : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ CD4+ Tem                     : num  0 0 0 0 0 ...
#  $ CD4+ memory T-cells          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ CD4+ naive T-cells           : num  0 2.42 0 0 0 ...
#  $ CD8+ T-cells                 : num  0 0 0 6.23 0 ...
#  $ CD8+ Tcm                     : num  0 2.99 0 11.33 0 ...
#  $ CD8+ Tem                     : num  0 0 0 1.39 0 ...
#  $ CD8+ naive T-cells           : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ CLP                          : num  0 2.23 0 4.73 0 ...
#  $ CMP                          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Chondrocytes                 : num  0 2.14 6.01 0 1.46 ...
#  $ Class-switched memory B-cells: num  0 3.06 0 2.92 0 ...
#  $ DC                           : num  0.921 0.826 0 0.669 0 ...
#  $ Endothelial cells            : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Eosinophils                  : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Epithelial cells             : num  32.6 18.5 23.2 21.1 12.5 ...
#  $ Erythrocytes                 : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Fibroblasts                  : num  0 0 7.45 0 0 ...
#  $ GMP                          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ HSC                          : num  13.69 18.62 21.17 5.49 17.65 ...
#  $ Hepatocytes                  : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Keratinocytes                : num  10.719 2.503 5.257 4.832 0.779 ...
#  $ MEP                          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ MPP                          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ MSC                          : num  16.5 13.3 22.5 15.7 40.5 ...
#  $ Macrophages                  : num  0.826 3.321 0 9.405 0.992 ...
#  $ Macrophages M1               : num  0 3.25 0 5.36 0 ...
#  $ Macrophages M2               : num  0 0 0 3.03 2.34 ...
#  $ Mast cells                   : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Megakaryocytes               : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Melanocytes                  : num  0.417 0.866 0.661 0 0.456 ...
#  $ Memory B-cells               : num  0 1.76 0 1.8 0 ...
#  $ Mesangial cells              : num  2.9 0 0 0 0 ...
#  $ Monocytes                    : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Myocytes                     : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ NK cells                     : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ NKT                          : num  6.14 23.57 11.7 18.22 11.73 ...
#  $ Neurons                      : num  0.118 0 0 0 0.118 ...
#  $ Neutrophils                  : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Osteoblast                   : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Pericytes                    : num  6.27 3.45 7.68 2.04 0 ...
#  $ Plasma cells                 : num  0 2.668 0 1.456 0.929 ...
#  $ Platelets                    : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Preadipocytes                : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Sebocytes                    : num  7.847 0.779 3.597 3.125 0.756 ...
#  $ Skeletal muscle              : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Smooth muscle                : num  4.86 8.8 4.31 0 10.85 ...
#  $ Tgd cells                    : num  0.362 0 0 0.362 0 ...
#  $ Th1 cells                    : num  3.78 4.09 3.08 9.35 0 ...
#  $ Th2 cells                    : num  9.88 14.36 0 20.87 0 ...
#  $ Tregs                        : num  1.55 1.68 0 1.65 0 ...
#  $ aDC                          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ cDC                          : num  2.64 4.86 0 4.12 2.68 ...
#  $ iDC                          : num  27.8 18.8 0 22.4 0 ...
#  $ ly Endothelial cells         : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ mv Endothelial cells         : num  1.81 3.57 3.01 2.8 3.63 ...
#  $ naive B-cells                : num  0 0.425 0 0.582 0 ...
#  $ pDC                          : num  0 2.52 0 8.33 0 ...
#  $ pro B-cells                  : num  0 0 0 0 0 0 0 0 0 0 ...

str(shamik4)
# 'data.frame':   235 obs. of  9 variables:
#  $ epithelial : num  0.488 0.419 0.405 0.282 0.537 ...
#  $ macrophage : num  0.1423 0.0904 0.1 0.2212 0.0416 ...
#  $ stroma     : num  0.239 0.161 0.322 0.123 0.225 ...
#  $ Bcell      : num  0.0359 0.1864 0.0554 0.2382 0.0629 ...
#  $ endothelial: num  0.0915 0.1366 0.118 0.1234 0.1333 ...
#  $ Tcell      : num  0.00422 0.00682 0 0.01198 0 ...
#  $ P-value    : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Correlation: num  0.636 0.815 0.675 0.721 0.63 ...
#  $ RMSE       : num  0.79 0.706 0.759 0.716 0.806 ...

cor(shamik1,shamik2)
#                Epi         Fib           B           NK        CD4T        CD8T          Mono       Neutro      Eosino
# Epi     1.00000000 -0.12510891 -0.56718207 -0.641695112 -0.69911538 -0.58758735 -0.3002717388  0.119255575 -0.12221404
# Fib    -0.12510891  1.00000000 -0.25832666 -0.041449873 -0.28265236 -0.25913492 -0.2033031455 -0.282604524  0.05124926
# B      -0.63511345 -0.28595595  0.93882868  0.358376482  0.63452483  0.47056955 -0.0008136691 -0.140119387  0.21009969
# NK     -0.66165536 -0.05858821  0.27512336  0.873580582  0.50221606  0.28447082  0.3095582943 -0.147077137 -0.05497555
# CD4T   -0.66737475 -0.29035518  0.56136675  0.450251334  0.90105550  0.44967559  0.1056996310 -0.111828281  0.06182697
# CD8T   -0.45799009 -0.27823699  0.33601195  0.256937349  0.34095444  0.88193147  0.0405770328 -0.122074243  0.20429416
# Mono   -0.49278158 -0.11914403  0.01477799  0.283837259  0.25976093  0.14389693  0.8497342912  0.001009755 -0.07886912
# Neutro -0.01675777 -0.36792270 -0.04847744  0.003742788 -0.04674963 -0.08982678  0.2592046246  0.852807938 -0.23136051
# Eosino -0.19539121  0.06420647  0.17033701  0.199242564  0.15979487  0.12800193 -0.1286073450 -0.098878241  0.38313073

cor(fracA,shamik4$"epithelial")
#[1] 0.6535578
cor(fracB,shamik4$"epithelial")
#[1] 0.6527951
cor(fracTum,shamik4$"epithelial")
#[1] 0.6685273

cor(fracA,shamik1$Epi)
#[1] 0.7663427
cor(fracB,shamik1$Epi)
#[1] 0.7893982
cor(fracTum,shamik1$Epi)
#[1] 0.7961338

plot(shamik1$Epi,fracA)



cor(fracTum,shamik1$Epi)
#[1] 0.7961338

cor(fracTum,shamik3$"Epithelial cells")
#[1] 0.008694469

cor(fracTum,shamik4$epithelial)
#[1] 0.6685273

cor(shamik1$Epi,shamik4$epithelial)
#[1] 0.6401505



i<-"cg09441966"

plot(fracTum,betaNew[i,sampleSets$samplesAll])
dev.new()
plot(shamik1$Epi,betaNew[i,sampleSets$samplesAll])

plot(shamik1$Epi,betaNew[i,sampleSets$samplesAll])

which(shamik1$Epi >40 & betaNew[i,sampleSets$samplesAll]>.25 & betaNew[i,sampleSets$samplesAll]<.4)
PD31140a 
      76 


lm(betaNew[i,sampleSets$samplesAll]~shamik1)


fracTum["PD31140a"]
PD31140a 
   0.584 

betaNew[i,"PD31140a"]
[1] 0.308

shamik1["PD31140a",]
              Epi      Fib        B       NK     CD4T CD8T     Mono   Neutro Eosino
PD31140a 46.38691 20.45182 1.963657 7.751208 7.342916    0 9.696002 6.407489      0

which(shamik1$Epi >45 & betaNew[i,sampleSets$samplesAll]<.01)

fracTum["PD31028a"]
PD31028a 
   0.623 

shamik1["PD31028a",]
              Epi Fib        B      NK     CD4T CD8T     Mono  Neutro    Eosino
PD31028a 73.80454   0 7.224819 5.80749 3.500272    0 5.307645 3.98621 0.3690218


shamik4[102,]
         epithelial macrophage    stroma    Bcell endothelial     Tcell P-value Correlation      RMSE
PD31172a  0.1425544 0.09071978 0.1009835 0.425716   0.0903661 0.1496602       0   0.8949096 0.5589753

shamik1[102,]        
              Epi      Fib        B       NK     CD4T     CD8T     Mono    Neutro   Eosino
PD31172a 3.146022 5.506036 26.90015 8.495769 25.32735 24.70876 4.478195 0.3503148 1.087397

summary(lm(betaNew[getProm(symbol="CD4")$probes[5],sampleSets$samplesAll]~shamik1$CD4T+shamik1$Epi))

plot(shamik1$CD4T/100,betaNew[getProm(symbol="CD4")$probes[5],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))
plot(shamik1$CD8T/100,betaNew[getProm(symbol="CD4")$probes[5],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))
plot(shamik1$Epi/100,betaNew[getProm(symbol="CD4")$probes[5],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))

plot(shamik1$CD4T/100,betaNew[getProm(symbol="CD8A")$probes[10],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))
plot(shamik1$CD8T/100,betaNew[getProm(symbol="CD8A")$probes[10],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))
plot(shamik1$Epi/100,betaNew[getProm(symbol="CD8A")$probes[10],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))

plot(shamik4$Tcell,betaNew[getProm(symbol="CD8A")$probes[10],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))
plot(shamik4$Tcell,betaNew[getProm(symbol="CD8A")$probes[10],sampleSets$samplesAll],ylim=c(0,1),xlim=c(0,1))

boxplot(shamik1/100,las=2)

boxplot(shamik4[,1:6],las=2)

shamik4$Tcell


fracTum["PD31172a"]
PD31172a 
   0.465 

fracA["PD31172a"]
# PD31172a 
#     0.26 

fracB["PD31172a"]
#  PD31172a 
# 0.6697512 

pdf("ascat_battenberg_byPloidy.pdf")
plot(fracA,fracB,col=1+(ploiA> 2.5),xlab="ascat tumor%",ylab="battenberg tumor%",pch=16)
legend("bottomright",legend="ascat ploidy>2.5",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(fracA,fracB,col=1+(ploiB> 2.5),xlab="ascat tumor%",ylab="battenberg tumor%",pch=16)
legend("bottomright",legend="battenberg ploidy>2.5",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(fracA,fracB,col=1+( abs(ploiB-ploiA) >.5 ),xlab="ascat tumor%",ylab="battenberg tumor%",pch=16)
legend("bottomright",legend="abs ploidy difference >0.5",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(ploiA,ploiB,col=1+( abs(ploiB-ploiA) >.5 ),xlab="ascat ploidy",ylab="battenberg ploidy",pch=16)
legend("bottomright",legend="abs ploidy difference >0.5",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(ploiA,ploiB,col=1+(fracA> fracB),xlab="ascat ploidy",ylab="battenberg ploidy",pch=16)
legend("bottomright",legend="purity ascat > batteberg",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(fracA,shamik1$Epi/100,col=1+(fracA> shamik1$Epi/100),xlab="ascat %",ylab="shamik EpiDish",pch=16)
legend("bottomright",legend="purity ascat > epi",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(fracB,shamik1$Epi/100,col=1+(fracB> shamik1$Epi/100),xlab="battenberg %",ylab="shamik EpiDish",pch=16)
legend("bottomright",legend="purity battenberg > epi",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(fracA,shamik4$"epithelial",col=1+(fracA> shamik4$"epithelial"),xlab="ascat %",ylab="shamik CibersortX",pch=16)
legend("bottomright",legend="purity ascat > ciber",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

plot(fracB,shamik4$"epithelial",col=1+(fracB> shamik4$"epithelial"),xlab="battenberg %",ylab="shamik CibersortX",pch=16)
legend("bottomright",legend="purity battenberg > ciber",col=2,pch=16)
abline(a=0,b=1,lwd=2,lty=2)

dev.off()

##higher purities from sequencing based estimates
  ##epithelial vs tumor cell content



###END