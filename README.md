# adjustBetas
This repository describes a simple approach for correcting for non-tumor infiltration in DNA methylation data derived using the Illumina Infinium 450/850K platform. The basic approach is described in our manuscript posted to [bioRxiv](https://doi.org/10.1101/2022.03.04.483052) in March 2022.

## The Issue
It is a well established fact that tumors are a heterogenous collection of cell types. In the context of DNA methylation analysis, true tumor methylation states may be confounded by the presence of infiltrating immune cells as well as non-malignant cells of the tissue in question, i.e. the tumor microenvironment (TME). Many approaches have been proposed to quantify tumor purity as well as the cellular composition of the TME. Some methods have also been proposed removing or reducing the influence of the TME on observed tumor DNA methylation estimates, similar to the aims of this work. Although the confounding of tumor DNA methylation estimates by non-tumor cells has been known for long, most analyses published to date have not addressed this issue. The aim of this work has been to devise a simple and widely applicable approach for producing tumor methylation estimates that more closely approximate the DNA methylation state of the malignant cellular compartment. To our knowledge there does not exist a single tool or approach that has gained universal acceptance as the way to accomplish this and although our approach will not be applicable in all situations, we believe it can be used to extract novel insights into e.g. the "pure" epigentic phenotypes of common human malignancies. 

## Approach
Our approach requires preexisting estimates of tumor purity derived from e.g. WGS/WES or the methylation data itself as well as a reasonably large cohort size. The main aim is to adjust tumor DNA methylation beta values for the influence of non-tumor cells. In cases where the tumor and non-tumor compartment have diametrically opposing methylation states, there is a linear influence of non-tumor methylation on observed beta values. To capture and correct for this, we iteratively apply the [FlexMix](https://cran.r-project.org/web/packages/flexmix) framework to each CpG on the platform to define (up to 3) sample populations with distinct linear relationships with tumor purity. Given a sufficiently large data set, inference of both the methylation state of the aggregate "normal" and tumor compartment can be made (see illustration from Staaf & Aine 2022 below).

![Figure 1 Staaf & Aine 2022](images/fig1.png?raw=true)

When applied to a large cohort of breast cancer tumors gathered by The Cancer Genome Atlas project (TCGA) our approach yields biologically and clinically sound results in terms of 1) the visual appearance of sample clustering, 2) a distribution of beta values that corresponds better to theretical expactiations, 3) an accurate approximation of the aggregate methylation state of the TME, and 4) a smaller influence of tumor purity (TME) on sample clustering.

![Figure 3 Staaf & Aine 2022](images/fig3.png?raw=true)

## Scripts 
The standalone function for adjusting Illumina beta values can be found in the file "function_correctBetas.r". The script "example_tcgaBrcaRun.r" illustrates the application of the described approach to the top 5000 most variable CpGs in a cohort of 630 breast cancer samples from TCGA with available data on the levels of 1) exome sequencing, 2) RNA-seq, 3) 450K methylation, 4) SNP6 copy number data. 

## Licence
You are free to use, modify or adapt the contents in any way you wish. 