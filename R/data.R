#' An example data set for testing nebula
#'
#' @name sample_data
#' @title An example data set for testing nebula
#' @description A dataset containing a count matrix, subject IDs, a data frame of predictors and scaling factors.
#' @format A list of four objects:
#' \describe{
#'   \item{count}{A raw count matrix}
#'   \item{sid}{A vector of subject IDs}
#'   \item{pred}{A data frame of three predictors}
#'   \item{offset}{A vector of scaling factors}
#' }
#' @docType data
#' @keywords datasets
NULL

#' An example Seurat object for testing scToNeb
#'
#' @name sample_seurat
#' @title An example Seurat object for testing scToNeb
#' @description A Seurat object containing a subset (1000 genes and 1000 cells) of the eight-pancreas scRNA-seq datasets for testing the scToNeb function.
#' @format A Seurat object
#' @source \url{https://github.com/satijalab/seurat-data}
#' @docType data
#' @keywords datasets
NULL
