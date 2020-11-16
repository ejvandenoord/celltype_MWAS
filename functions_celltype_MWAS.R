
run_celltypeMWAS = function( param  ) {
  
  # bulk MWAS
  ramwas5MWAS(param)
  qqPlotFast(getMWAS(param)$`p-value`)
  title(paste(param$modeloutcome, "_nCov", length(param$modelcovariates), "_nPC", param$modelPCs))
  
  # cell-type MWASs
  ctnames = param$celltypeMWAS
  cinames = paste0(ctnames,'_x_',param$modeloutcome)
  ccvar   = param$covariates[[param$modeloutcome]]
  
  for ( i in seq_along(ctnames) ) param$covariates[[cinames[i]]] = param$covariates[[ctnames[i]]] * ccvar
  
  for ( i in seq_along(ctnames) ) {
    dropct_i = which( ctnames[i] == ctnames   )
    message('Running cell-type MWAS for: ', cinames[dropct_i])
    pp = param
    pp$modeloutcome = cinames[dropct_i]
    pp$modelcovariates = unique(c(param$modelcovariates,ctnames,cinames[-dropct_i]))
    pp$modelhasconstant = FALSE 
    pp$dirmwas      = paste0(pp$dirmwas,"_",ctnames[dropct_i])
    
    ramwas5MWAS(pp)
    qqPlotFast(getMWAS(pp)$`p-value`)
    title(paste(cinames[dropct_i], "_nCov", length(modelcovariates), "_nPC", modelPCs))
    
  }
  
}

