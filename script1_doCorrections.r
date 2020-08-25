#####======================================================================#####
### Get TCGA/GDC data for tumor type of choice
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

################################################################################
##Set home directory

##downloads will require substantial amount of space, e.g. 1-2Tb for TCGA-BRCA

##set/create own home directory below:

##work
dir.create("~/hdd1/adjustBetas")
HOME<-"~/hdd1/adjustBetas"
GIT<-"~/Documents/adjustBetas"
##home
# dir.create("I:/data/adjustBetas")
# HOME<-"I:/data/adjustBetas"
# GIT<-"F:/gitProjects/adjustBetas"

paste0(HOME,"/data")


################################################################################
##load required packages

library(GenomicRanges)

library(RColorBrewer)

library(pheatmap)

library(magick)

library(ggalluvial)

#install.packages("flexmix")
library("flexmix")

################################################################################
##Define functions



################################################################################
##Load data manifests and parse

q("no")
###END