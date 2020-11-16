

rm(list=ls())


dirproject      = "G:\\My Drive\\SOP\\methylation\\gsms\\exposure\\papers\\puberty\\paper\\JCPP\\revision\\code\\"

covariate_data  = read.csv(paste0(dirproject,"\\tanner.csv") ) 
modeloutcome    = "imtanner2_j4ndep_jsex" 

modelPCs        = 1 

celltypeMWAS    = c('CD03','CD14','CD15','CD19')


clincal_covariates = c("imtanner2","j4ndep","jsex","imtanner2_j4ndep","j4ndep_jsex","jage","jage2","race_aa","race_ai")
technical_covariates = c("PeakSQRT","AUC","AvgCpGcoverage","min_thres","peak2","peak3","peak4","peak5","ReadsUsedForCoverage","peak_height","skewness","NonCpg2CpGcoverageRatio","max_thres","NonCpGreadsPct","peak","r2_eb2")
bulk_celltype_covariates = c("CD03","CD14","CD15") # one cell-type deleted as celltypes sum to 1

functions_script = paste0(dirproject,"functions_celltype_MWAS.R")

### end settings

library(ramwas)
setwd(dirproject)
source(functions_script)

modelcovariates = c(clincal_covariates,technical_covariates,bulk_celltype_covariates)

param = ramwasParameters(
  dircoveragenorm = 'rw',
  covariates      = covariate_data,
  modelcovariates = modelcovariates,
  modeloutcome    = modeloutcome,
  modelPCs        = modelPCs,
  celltypeMWAS    = celltypeMWAS  )

param = parameterPreprocess(param)

if (!file.exists(paste0(param$dirpca,'\\eigenvalues.bmat'))) ramwas4PCA(param) 

# cell type MWAS
run_celltypeMWAS(param)
