
#' Connect to loom file
#'
#' @param loom_path A single-element character vector that indicates the path
#'  to the transcriptome loom file.
#'
#' @returns An object of class 'loom' that represents a connection to the loom
#'  file stored on disk.
#'
#' @examples
#' lpath <- "Path/to/loom"
#'
#' lfile <- read_loom(lpath)
read_loom <- function(loom_path) {
  tryCatch({loomR::connect(filename = loom_path, mode = "r", skip.validate = TRUE)},
           error = function(e) stop(e))
}

#' Convert loom file to Seuratobject
#'
#' @param loom_path A single-element character vector that indicates the path
#'  to the transcriptome loom file.
#' @param matrix_rowname_col A single element character vector that indicates
#'  the variable (i.e. column) of the row (i.e. gene) metadata in the loom file.
#'  This variable will subsequently be used as the feature names for the
#'  generated sparse matrix.
#' @param matrix_colname_col A single element character vector that indicates
#'  the variable (i.e. column) of the column (i.e. cell) metadata in the loom
#'  file. This variable will subsequently be used as the cell names for the
#'  generated sparse matrix.
#' @param seurat_assay_name A single element character vector representing the
#'  name of the assay in the Seuratobject.
#' @param seurat_min_cells Include only features in the seuratobject that are
#'  detected in at least this many cells.
#' @param seurat_names_field Argument used to get sample name from the cell
#'  names. Works in conjunction with the 'seurat_names_delim' argument.
#'   'seurat_names_field' is an integer selecting which substring
#'   (i.e. delimited by 'seurat_names_delim') should be selected as the sample
#'   name. The default is 1L (i.e. the first field of the cell name).
#' @param seurat_names_delim Argument used to get sample name from the cell
#'  names. Works in conjunction with the 'seurat_names_field' argument.
#'  'seurat_names_delim' is a single element character vector that represents
#'  the delimiter used to select the cell name substrings.
#' @param ... Additional arguments will be given to the
#'  SeuratObject::CreateSeuratObject() function call.
#'
#' @returns A Seuratobject populated with a sparse matrix representing the
#'  transcriptome data from the loom file located at 'loom_path'.
#' @export
#'
#' @examples
#' lpath <- "Path/to/loom"
#'
#' lseurat <- LoomAsSeurat(lpath,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   seurat_assay_name = "RNA",
#'   seurat_min_cells = 0,
#'   seurat_names_field = 1L,
#'   seurat_names_delim = "_")
LoomAsSeurat <- function(loom_path,
                         matrix_rowname_col = "Gene",
                         matrix_colname_col = "CellID",
                         seurat_assay_name = "RNA",
                         seurat_min_cells = 0,
                         seurat_names_field = 1L,
                         seurat_names_delim = "_",
                         ...) {

  ### Connect to loom file
  lfile <- withCallingHandlers(
    tryCatch(read_loom(loom_path), error = function(e) {
      stop(sprintf("Error opening loom file: %s\nOriginal error message:\n%s",
                   loom_path, conditionMessage(e)), call. = FALSE)
    }),
    warning = function(w) {
      message(sprintf("Warning while opening loom file: %s", loom_path))
      message("Original warning message:\n")
      message(conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Add a check for successful connection
  if (!inherits(lfile, "loom")) {
    stop(sprintf("Failed to connect to loom file at: %s. See above messages for details.", loom_path), call. = FALSE)
  }

  # Close loom files upon exit of function
  on.exit(lfile$close_all(), add = TRUE)

  ### Collect the metadata and matrix
  row_meta <- lfile$get.attribute.df(MARGIN = 1)
  col_meta <- lfile$get.attribute.df(MARGIN = 2)
  lmatrix <- lfile[["matrix"]]

  ### Convert to sparse matrix
  lmatrix_sparse <- t(Matrix::Matrix(lmatrix[,], sparse = TRUE))

  ### Perform checks
  if (is.null(row_meta[[matrix_rowname_col]])) {
    stop(sprintf("'%s' is not an existing column in the row metadata of the loom file at path '%s'", matrix_rowname_col, loom_path))
  }

  if (is.null(col_meta[[matrix_colname_col]])) {
    stop(sprintf("'%s' is not an existing column in the column metadata of the loom file at path '%s'", matrix_colname_col, loom_path))
  }

  if (dim(lmatrix_sparse)[1] != length(row_meta$Gene)) {
    stop(sprintf("The row number for the matrix and the row number for the row-metadata are not equal for the loom file at path '%s' has ", matrix_colname_col, loom_path))
  }

  if (dim(lmatrix_sparse)[2] != length(col_meta$CellID)) {
    stop(sprintf("The column number for the matrix and the row number for the column-metadata are not equal for the loom file at path '%s' has ", matrix_colname_col, loom_path))
  }

  ### Add row and column names to sparse matrix
  rownames(lmatrix_sparse) <- row_meta[[matrix_rowname_col]]
  colnames(lmatrix_sparse) <- col_meta[[matrix_colname_col]]

  ### Create the seurat object
  lseurat <- SeuratObject::CreateSeuratObject(counts = lmatrix_sparse,
                                assay = seurat_assay_name,
                                min.cells = seurat_min_cells,
                                names.field = seurat_names_field,
                                names.delim = seurat_names_delim,
                                ...)

  ### Add the cell metadata
  lseurat[[matrix_colname_col]] <- rownames(lseurat[[]])
  col_meta_merged <- merge(x = lseurat[[]], y = col_meta, by = matrix_colname_col)
  rownames(col_meta_merged) <- col_meta_merged[[matrix_colname_col]]
  # Reorder based on current seurat metadata, makes sure the same order of rows
  # is maintained
  col_meta_merged <- col_meta_merged[order(match(col_meta_merged[[matrix_colname_col]],
                                                 rownames(lseurat[[]]))),]
  # Add back the metadata
  lseurat <- SeuratObject::AddMetaData(lseurat, metadata = col_meta_merged)

  ### Add the row metadata
  row_meta_ordered <- row_meta[order(match(row_meta[[matrix_rowname_col]], rownames(lseurat[[seurat_assay_name]]))),]
  lseurat[[seurat_assay_name]] <- SeuratObject::AddMetaData(lseurat[[seurat_assay_name]], metadata = row_meta_ordered)

  return(lseurat)
}

#' Convert multiple loom files to a single merged Seuratobject
#'
#' @param loom_paths A character vector that indicates the paths to the
#'  transcriptome loom files.
#' @param matrix_rowname_col A single element character vector that indicates
#'  the variable (i.e. column) of the row (i.e. gene) metadata in the loom file.
#'  This variable will subsequently be used as the feature names for the
#'  generated sparse matrix.
#' @param matrix_colname_col A single element character vector that indicates
#'  the variable (i.e. column) of the column (i.e. cell) metadata in the loom
#'  file. This variable will subsequently be used as the cell names for the
#'  generated sparse matrix.
#' @param seurat_assay_name A single element character vector representing the
#'  name of the assay in the Seuratobject.
#' @param seurat_min_cells Include only features in the seuratobject that are
#'  detected in at least this many cells.
#' @param seurat_names_field Argument used to get sample name from the cell
#'  names. Works in conjunction with the 'seurat_names_delim' argument.
#'   'seurat_names_field' is an integer selecting which substring
#'   (i.e. delimited by 'seurat_names_delim') should be selected as the sample
#'   name. The default is 1L (i.e. the first field of the cell name).
#' @param seurat_names_delim Argument used to get sample name from the cell
#'  names. Works in conjunction with the 'seurat_names_field' argument.
#'  'seurat_names_delim' is a single element character vector that represents
#'  the delimiter used to select the cell name substrings.
#' @param ... Additional arguments will be given to the
#'  SeuratObject::CreateSeuratObject() function call.
#'
#' @returns A Seuratobject populated with a sparse matrix representing the
#'  transcriptome data from the loom files located at 'loom_paths'.
#'  The 'orig.idents' variable in the cell metadata is used to identify from
#'  which sample each cell came from in the merged object. This 'orig.idents'
#'  gets chosen based on the 'seurat_names_field' and 'seurat_names_delim'
#'  arguments.
#' @export
#'
#' @examples
#' lpaths <- c("Path/to/loomfile1",
#'             "Path/to/loomfile2",
#'             "Path/to/loomfile3")
#'
#' lseurat <- LoomAsSeurat(lpath,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   seurat_assay_name = "RNA",
#'   seurat_min_cells = 0,
#'   seurat_names_field = 1L,
#'   seurat_names_delim = "_")
MultiLoomAsSeurat <- function(loom_paths,
                              matrix_rowname_col = "Gene",
                              matrix_colname_col = "CellID",
                              seurat_assay_name = "RNA",
                              seurat_min_cells = 0,
                              seurat_names_field = 1L,
                              seurat_names_delim = "_",
                              ...) {

  # Check if file paths exist
  if (any(file.exists(loom_paths) == FALSE)) {
    non_existing_idx <- which(file.exists(loom_paths) == FALSE)
    stop(sprintf("The following files do not exist:\n%s", paste(loom_paths[non_existing_idx], collapse = "\n")))
  }

  loom_seurat_list <- lapply(loom_paths, function(lpath)
    LoomAsSeurat(lpath,
                 matrix_rowname_col = matrix_rowname_col,
                 matrix_colname_col = matrix_colname_col,
                 seurat_assay_name = seurat_assay_name,
                 seurat_min_cells = seurat_min_cells,
                 seurat_names_field = seurat_names_field,
                 seurat_names_delim = seurat_names_delim,
                 ...)
  )

  # Merge seurat objects if neccesary
  if (length(loom_seurat_list) > 1) {
    merged_seurat <- merge(x = loom_seurat_list[[1]], y = loom_seurat_list[[-1]])
    return(merged_seurat)
  } else if (length(loom_seurat_list) == 1) {
    return(loom_seurat_list[[1]])
  } else {
    stop("Loom_seurat_list is empty, meaning that the LoomAsSeurat function failed to produce any seurat object")
  }
}
