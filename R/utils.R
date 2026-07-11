## other useful functions

check_reml = function(reml, model)
{
  if((reml!=0) & (reml!=1))
  {stop("reml must be either 0 or 1.")}
  
  if((reml==1) & (model!='NBLMM'))
  {
    reml <- 0
    warning("The value of reml is changed to zero because reml=1 is only supported for NBLMM in the current version.")
  }
  reml
}

check_conv = function(repml, conv, nb, vare, min, max, cutoff = 1e-8)
{
  if(conv==1)
  {
    if(vare[1]==max[1] | vare[2]==min[2])
    {
      conv = -60
    }else{
      if(is.nan(repml$loglik))
      {conv = -30}else{
        if(repml$iter==50)
        {
          conv = -20
        }else{
          if(repml$damp==11)
          {conv = -10}else{
            if(repml$damp==12)
            {conv = -40}
          }
        }
      }
    }
  }
  
  if(nb>1)
  {
    if(min(eigen(repml$var)$values) < cutoff)
    {conv = -25}
  }
  
  conv
}

#' Load sample_seurat data
#'
#' @title Load sample_seurat example data
#' @description Loads the sample_seurat example dataset for testing the scToNeb function.
#' This dataset contains a subset (1000 genes and 1000 cells) of the eight-pancreas scRNA-seq datasets.
#'
#' @return A Seurat object containing example single-cell data
#'
#' @details This function requires the Seurat package to be installed.
#' If Seurat is not available, the function will throw an error with installation instructions.
#'
#' @examples
#' \dontrun{
#' library(nebula)
#' sample_seurat <- load_sample_seurat()
#' re <- scToNeb(obj = sample_seurat, assay = "RNA", id = "replicate",
#'               pred = c("celltype", "tech"))
#' }
#'
#' @export
#' @seealso \code{\link{scToNeb}}
#' @keywords data
load_sample_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The 'sample_seurat' dataset requires the Seurat package. ",
         "Please install it with: install.packages('Seurat')")
  }

  data_file <- system.file("extdata", "sample_seurat.rda", package = "nebula")
  if (data_file == "" || !file.exists(data_file)) {
    stop("Could not find sample_seurat data file. ",
         "Please reinstall the nebula package.")
  }

  # Load the data into a temporary environment
  env <- new.env()
  load(data_file, envir = env)

  # Return the sample_seurat object
  if (exists("sample_seurat", envir = env)) {
    return(get("sample_seurat", envir = env))
  } else {
    stop("sample_seurat object not found in data file.")
  }
}