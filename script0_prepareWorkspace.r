######================================================######
### script 0 - prepare github workspace                  ###
######================================================######

cd ~/Documents

##The full workspace can be cloned from github using the command
'git clone https://github.com/StaafLab/adjustBetas.git'

##no need to run commands below if you use the gitHub-link

############################################################
##Main folder

dir.create("~/hdd1/adjustBetas")
HOME<-"~/hdd1/adjustBetas"
GIT<-"~/Documents/adjustBetas"
##home
# dir.create("I:/data/adjustBetas")
# HOME<-"I:/data/adjustBetas"
# GIT<-"F:/gitProjects/adjustBetas"

############################################################
##create subfolders in main folder

#data
dir.create(paste0(HOME,"/data"),recursive=T)

##Following data files required in "data"-folder 
#"GSE67919_Beta.RData" #methylation data for 96 Normal breast samples
#"GSE67919_Annotations.RData" #annotations for 96 Normal breast samples

#"20191021_workspace_tcgaAtacBrca.RData" #main workspace for JS TNBC set

#results
dir.create(paste0(HOME,"/results"),recursive=T)

q("no")
###END
