
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

collect_results_old = function(param,analysis_label,collect_files) {
  
  # this RaMWAS comment adds the names of directories where RaMWAS will write results most relevant are param$dirpca and param$dirmwas
  param        = parameterPreprocess(param)
  
  results_dir = paste0(param$dirproject,"/",analysis_label,"/")
  dir.create(results_dir,showWarnings=F,recursive=T)
  
  collect_files = makefullpath(param,collect_files) # this function finds the exact location of the  files 
  
  file.copy( collect_files,results_dir) # this copies the files we need
  
  # this zips the files we need
  zip_file = paste0(results_dir,analysis_label,"_results.zip")
  zip::zip(zip_file,basename(collect_files), root=results_dir )
  message("Requested results archived in ",  zip_file )  
  
}


collect_results = function(param,analysis_label,collect_files) {
  
  # this RaMWAS comment adds the names of directories where RaMWAS will write results most relevant are param$dirpca and param$dirmwas
  param        = parameterPreprocess(param)
  
  results_dir = paste0(param$dirproject,"/",analysis_label,"/bulk")
  dir.create(results_dir,showWarnings=F,recursive=T)
  
  collect_files_orig = collect_files
  if ("Stats_and_pvalues.bmat" %in% collect_files) {
      loc_files = c("CpG_chromosome_names.txt","CpG_locations.bmat","CpG_locations.desc.txt","CpG_locations.nmscol.txt")
      collect_files = makefullpath(param,c(collect_files,loc_files)) # this function finds the exact location of the  files 
  }
    
  collect_files = makefullpath(param,collect_files)
  file.copy( collect_files,results_dir,overwrite = T) # this copies the files we need
  
  # this zips the files we need
  zip_file = paste0(param$dirproject,"/",analysis_label,"/bulk.zip")
  zip::zip(zip_file,basename(collect_files), root=results_dir )
  message("Bulk results archived in ",  zip_file )  
  
 
  # cell-type MWASs
  if ( !is.null(param$celltypeMWAS) ) {
    
    ctnames = param$celltypeMWAS
    for ( i in seq_along(ctnames) ) {  # i=3
      collect_files = collect_files_orig  
      results_dir = paste0(param$dirproject,"/",analysis_label,"/",ctnames[i])
      dir.create(results_dir,showWarnings=F,recursive=T)
    
      pp            = param
      pp$dirmwas    = paste0(pp$dirmwas,"_",ctnames[i])
      collect_files = makefullpath(pp,collect_files) # this function finds the exact location of the  files 
    
      file.copy( collect_files,results_dir,overwrite = T) # this copies the files we need
    
      # this zips the files we need
      zip_file = paste0(param$dirproject,"/",analysis_label,"/",ctnames[i],".zip")
      zip::zip(zip_file,basename(collect_files), root=results_dir )
      message(ctnames[i]," results archived in ",  zip_file )  
    }
  }
}



makefullpath = function(param,collect_files) {
  # param = pp
  dirmwas_files = list.files(param$dirmwas,full.names=TRUE)
  dirpca_files  = list.files(param$dirpca,full.names=TRUE)
  dircov_files  = list.files(paste0(param$dirproject,"/covariates/"),full.names=TRUE)
  dirmwas_files = list.files(param$dirmwas,full.names=TRUE)
  dircoveragenorm_files = list.files(param$dircoveragenorm,full.names=TRUE)
  
  for (i in seq_along(collect_files)) { # i=1
    
    sel = grepl(collect_files[i],basename(dirmwas_files) ) # starts looking here first
    if (any(sel)) collect_files[i] = dirmwas_files[sel]
    
    sel = grepl(collect_files[i],basename(dirpca_files) )
    if (any(sel)) collect_files[i] = dirpca_files[sel]
    
    sel = grepl(collect_files[i],basename(dircov_files) )
    if (any(sel)) collect_files[i] = dircov_files[sel]
    
    sel = grepl(collect_files[i],basename(dircoveragenorm_files) )
    if (any(sel)) collect_files[i] = dircoveragenorm_files[sel]
    
     
    
  }
  
  collect_files
}