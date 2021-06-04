#####======================================================================#####
### Correct TCGA BRCA beta values for infiltration
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

################################################################################
##Set home directory

##set/create own home directory below:

##work
HOME<-"~/hdd1/tcgaBrca"
MANIFEST<-"~/Documents/getTcga/manifest"
##home
HOME<-"I:/data/tcgaBrca"
MANIFEST<-"F:/gitProjects/getTcga/manifest"

##tumor type
TUMOR_TYPE<-"brca"

list.files(MANIFEST,full.names=T)
#[1] "/home/med-mai/Documents/getTcga/manifest/atac"  
#[2] "/home/med-mai/Documents/getTcga/manifest/brca"  
#[3] "/home/med-mai/Documents/getTcga/manifest/luad"  
#[4] "/home/med-mai/Documents/getTcga/manifest/lusc"  
#[5] "/home/med-mai/Documents/getTcga/manifest/pancan"

##create data directories
#dir.create(paste0(HOME,"/","me/norm"),recursive=TRUE)

################################################################################
##load required packages

if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer") }

library("RColorBrewer")

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap") }

library("pheatmap")

if(!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr") }

library("dplyr")

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges") }

library("GenomicRanges")

################################################################################
##Get data object

load(file=paste0(HOME,"/","finalWorkspace_atacCnGexMeWes_withAnnotations.RData"))

ls()
#  [1] "annoObj"    "betaAdj"    "betaNorm"   "betaOrig"   "dataAtac"   "dataCn"     "dataMut"    "dataSeg"    "gexCounts"  "gexFpkm"   
# [11] "gexUq"      "HOME"       "MANIFEST"   "sampleAnno" "TUMOR_TYPE"

################################################################################
##







q("no")
###END