# Test cases for scToNeb function

test_that("scToNeb gives helpful error when Seurat not installed", {
  # Create a mock object with Seurat class but without Seurat installed
  # This tests the package availability check
  mock_obj <- structure(list(), class = "Seurat")

  # Temporarily hide Seurat from namespace if it exists
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Skip this test if Seurat is actually installed
    skip("Seurat is installed, cannot test missing package error")
  }

  expect_error(
    scToNeb(mock_obj),
    "Package 'Seurat' is required for Seurat objects"
  )
})

test_that("scToNeb gives helpful error when SingleCellExperiment not installed", {
  # Create a mock object with SingleCellExperiment class
  mock_obj <- structure(list(), class = "SingleCellExperiment")

  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    skip("SingleCellExperiment is installed, cannot test missing package error")
  }

  expect_error(
    scToNeb(mock_obj),
    "Package 'SingleCellExperiment' is required for SingleCellExperiment objects"
  )
})

test_that("scToNeb rejects unsupported object types", {
  # Test with a completely unsupported object type
  unsupported_obj <- list()

  expect_error(
    scToNeb(unsupported_obj),
    "Please provide either a SingleCellExperiment or Seurat object"
  )
})

# Test with actual Seurat object (only run if Seurat is installed)
test_that("scToNeb works with Seurat objects", {
  skip_if_not_installed("Seurat")

  # Load the sample data using the helper function
  sample_seurat <- load_sample_seurat()

  # Test basic functionality
  result <- scToNeb(obj = sample_seurat, assay = "RNA",
                    id = "replicate", pred = c("celltype", "tech"))

  # Check that result is a list
  expect_type(result, "list")

  # Check that required components exist
  expect_true(all(c("count", "pred") %in% names(result)))

  # Check that id is included when requested
  expect_true("id" %in% names(result))

  # Check that count matrix is present (either matrix or Matrix class)
  expect_true(is.matrix(result$count) || inherits(result$count, "Matrix"))
})

test_that("scToNeb works with Seurat object from latest Seurat version", {
  skip_if_not_installed("Seurat")

  library(Seurat)

  # Create a fresh Seurat object using modern Seurat API
  # This tests compatibility with the latest Seurat version

  # Create sample count matrix
  count_matrix <- matrix(
    rpois(200, lambda = 10),
    nrow = 20,
    ncol = 10
  )
  rownames(count_matrix) <- paste0("gene", 1:20)
  colnames(count_matrix) <- paste0("cell", 1:10)

  # Create metadata
  metadata <- data.frame(
    sample_id = rep(c("A", "B"), each = 5),
    celltype = rep(c("T", "B", "NK", "Mono", "DC"), 2),
    nCount_RNA = colSums(count_matrix),
    orig_ident = rep("sample1", 10),
    row.names = colnames(count_matrix)
  )

  # Create Seurat object using modern API
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    meta.data = metadata,
    project = "test_project"
  )

  # Test with the modern Seurat object
  result <- scToNeb(
    obj = seurat_obj,
    assay = "RNA",
    id = "sample_id",
    pred = c("celltype", "orig_ident"),
    offset = "nCount_RNA"
  )

  # Verify the result structure
  expect_type(result, "list")
  expect_true(all(c("count", "pred", "id", "offset") %in% names(result)))

  # Verify count matrix dimensions match input
  expect_equal(dim(result$count), dim(count_matrix))

  # Verify predictor data frame has correct rows
  expect_equal(nrow(result$pred), ncol(count_matrix))

  # Verify id is character vector
  expect_type(result$id, "character")

  # Verify offset is numeric/double
  expect_type(result$offset, "double")

  # Test without specifying assay (should use default)
  result_default <- scToNeb(
    obj = seurat_obj,
    id = "sample_id",
    pred = "celltype"
  )

  expect_true("count" %in% names(result_default))
  expect_true("pred" %in% names(result_default))
})

test_that("scToNeb handles Seurat object with sparse matrix", {
  skip_if_not_installed("Seurat")

  library(Seurat)

  # Create Seurat object with sparse matrix (common in modern Seurat)
  count_matrix <- Matrix::Matrix(
    rpois(300, lambda = 5),
    nrow = 30,
    ncol = 10,
    sparse = TRUE
  )
  rownames(count_matrix) <- paste0("gene", 1:30)
  colnames(count_matrix) <- paste0("cell", 1:10)

  metadata <- data.frame(
    batch = rep(c("batch1", "batch2"), 5),
    condition = rep(c("ctrl", "trt"), 5),
    row.names = colnames(count_matrix)
  )

  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    meta.data = metadata
  )

  result <- scToNeb(
    obj = seurat_obj,
    assay = "RNA",
    id = "batch",
    pred = "condition"
  )

  expect_type(result, "list")
  expect_true("count" %in% names(result))
  expect_true("pred" %in% names(result))
  expect_true("id" %in% names(result))
})

test_that("scToNeb handles offset parameter with Seurat", {
  skip_if_not_installed("Seurat")

  data("sample_seurat", package = "nebula")

  # Test with offset
  result <- scToNeb(obj = sample_seurat, assay = "RNA",
                    id = "replicate", pred = "celltype",
                    offset = "nCount_RNA")

  expect_true("offset" %in% names(result))
  expect_type(result$offset, "double")
})

test_that("scToNeb handles missing predictor variable gracefully", {
  skip_if_not_installed("Seurat")

  data("sample_seurat", package = "nebula")

  # Test with a non-existent predictor
  expect_error(
    scToNeb(obj = sample_seurat, assay = "RNA",
            id = "replicate", pred = "nonexistent_var"),
    "not available in the meta data"
  )
})

test_that("scToNeb handles missing id variable gracefully", {
  skip_if_not_installed("Seurat")

  data("sample_seurat", package = "nebula")

  # Test with a non-existent id
  expect_error(
    scToNeb(obj = sample_seurat, assay = "RNA",
            id = "nonexistent_id", pred = "celltype"),
    "not available in the meta data"
  )
})

test_that("scToNeb handles missing offset variable gracefully", {
  skip_if_not_installed("Seurat")

  data("sample_seurat", package = "nebula")

  # Test with a non-existent offset
  expect_error(
    scToNeb(obj = sample_seurat, assay = "RNA",
            id = "replicate", pred = "celltype",
            offset = "nonexistent_offset"),
    "not available in the meta data"
  )
})

test_that("scToNeb warns when no sample ID is provided", {
  skip_if_not_installed("Seurat")

  data("sample_seurat", package = "nebula")

  # Test warning when id is NULL
  expect_warning(
    scToNeb(obj = sample_seurat, assay = "RNA",
            pred = "celltype"),
    "No sample ID provided"
  )
})

test_that("scToNeb uses default assay when none specified", {
  skip_if_not_installed("Seurat")

  data("sample_seurat", package = "nebula")

  # Test with verbose=FALSE to suppress the default assay message
  result <- scToNeb(obj = sample_seurat, id = "replicate", pred = "celltype", verbose = FALSE)
  expect_true("count" %in% names(result))
  expect_true("pred" %in% names(result))
})

# Test with SingleCellExperiment object (only run if installed)
test_that("scToNeb works with SingleCellExperiment objects", {
  skip_if_not_installed("SingleCellExperiment")

  # Create a minimal SingleCellExperiment object for testing
  library(SingleCellExperiment)

  # Create sample data
  count_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  colnames(count_matrix) <- paste0("cell", 1:10)
  rownames(count_matrix) <- paste0("gene", 1:10)

  # Create metadata
  metadata <- data.frame(
    sample_id = rep(c("A", "B"), 5),
    celltype = rep(c("T", "B"), 5),
    row.names = colnames(count_matrix)
  )

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = count_matrix),
    colData = metadata
  )

  # Test basic functionality
  result <- scToNeb(obj = sce, id = "sample_id", pred = "celltype")

  expect_type(result, "list")
  expect_true(all(c("count", "pred", "id") %in% names(result)))
})

test_that("scToNeb handles offset with SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")

  library(SingleCellExperiment)

  count_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  colnames(count_matrix) <- paste0("cell", 1:10)
  rownames(count_matrix) <- paste0("gene", 1:10)

  metadata <- data.frame(
    sample_id = rep(c("A", "B"), 5),
    celltype = rep(c("T", "B"), 5),
    offset_factor = runif(10),
    row.names = colnames(count_matrix)
  )

  sce <- SingleCellExperiment(
    assays = list(counts = count_matrix),
    colData = metadata
  )

  result <- scToNeb(obj = sce, id = "sample_id", pred = "celltype",
                   offset = "offset_factor")

  expect_true("offset" %in% names(result))
  expect_type(result$offset, "double")
})

test_that("scToNeb handles missing variables in SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")

  library(SingleCellExperiment)

  count_matrix <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  colnames(count_matrix) <- paste0("cell", 1:10)

  metadata <- data.frame(
    sample_id = rep(c("A", "B"), 5),
    row.names = colnames(count_matrix)
  )

  sce <- SingleCellExperiment(
    assays = list(counts = count_matrix),
    colData = metadata
  )

  # Test with missing predictor
  expect_error(
    scToNeb(obj = sce, id = "sample_id", pred = "missing_var"),
    "not available in the meta data"
  )

  # Test with missing id
  expect_error(
    scToNeb(obj = sce, id = "missing_id", pred = "sample_id"),
    "not available in the meta data"
  )

  # Test with missing offset
  expect_error(
    scToNeb(obj = sce, id = "sample_id", pred = "sample_id",
            offset = "missing_offset"),
    "not available in the meta data"
  )
})
