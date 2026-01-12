#' Create a custom "lseurat" object.
#'
#' This function acts as a constructor for a custom "lseurat" object,
#' which is a named list containing a Seurat object, a list of parameters,
#' and an optional gmm object.
#'
#' @param seurat_obj A valid SeuratObject. This is a required argument.
#' @param params A named list containing parameters. Defaults to an empty list.
#' @param gmm An optional gmm model object of class 'Mclust'. Defaults to NA.
#' @param loom_path A single-element character vector, representing the path to
#'  the loom file from which this object was generated.
#' @param id A single-element character vector. It is the unique ID chosen to
#'  represent this object (e.g. "i1_ES"). Defaults to NA.
#'
#' @return A named list object of class "lseurat".
#'
#' @export
#'
#' @examples
#' # Example usage with dummy objects that have the correct class.
#' # We'll create simple lists and assign the correct classes to them for the example.
#'
#' # Create a dummy SeuratObject
#' my_seurat <- list()
#' class(my_seurat) <- "Seurat"
#'
#' # Create a dummy Mclust object
#' my_gmm <- list(bic = -434.9987, components = 2)
#' class(my_gmm) <- "Mclust"
#'
#' # Define some parameters
#' all_params <- list(min_mapping_quality = 40, pdf_threshold = 0.0005)
#'
#' # Scenario 1: Create a lseurat object with all arguments provided
#' lseurat_obj_full <- create_lseurat_object(
#'   seurat_obj = my_seurat,
#'   params = all_params,
#'   gmm = my_gmm,
#'   loom_path = "/path/to/my/file.loom",
#'   id = "i5-TDT3-d17"
#' )
#'
#' # Scenario 2: Create a lseurat object without gmm, loom_path, or id
#' # (These arguments will use their default NA values)
#' lseurat_obj_basic <- create_lseurat_object(
#'   seurat_obj = my_seurat
#' )
#'
#' # Check the class of the new objects
#' class(lseurat_obj_full)
#' class(lseurat_obj_basic)
#'
#' # You can also check their structure
#' \dontrun{
#'   str(lseurat_obj_full)
#'   str(lseurat_obj_basic)
#' }
create_lseurat_object <- function(seurat_obj, params = list(), gmm = NA,
                                  loom_path = NA, id = NA) {

  # Function to check if object is a single-element character vector or NA
  is_1length_character <- function(x) {
    ((is.atomic(x) && length(x) == 1 && inherits(x, "character")) || identical(x, NA))
  }

  # Ensure 'seurat_obj' is provided and has the correct class
  if (missing(seurat_obj)) {
    stop("Argument 'seurat_obj' is required.")
  }
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Argument 'seurat_obj' must be an object of class 'Seurat'.")
  }

  # Ensure 'params' is a list
  if (!is.list(params)) {
    warning("Argument 'params' is not a list and has been coerced to one.")
    params <- as.list(params)
  }

  # Ensure 'gmm' is either NA or an Mclust object
  if (!identical(gmm, NA) && !inherits(gmm, "Mclust")) {
    warning(paste0("Argument 'gmm' is not NA and does not have the expected ",
                   "class 'Mclust'.\n'gmm' has been set to NA."))
    gmm <- NA
  }

  # Ensure 'loom_path' is a single-element character vector
  if (!is_1length_character(loom_path)) {
    warning(paste0("While creating lseurat object.\nInvalid 'loom_path' ",
                   "argument. 'loom_path should be either a 1-length atomic vector of ",
                   "type 'character' or NA.\nSetting 'loom_path' to NA"))
    loom_path <- NA
  }

  # Ensure 'id' is a single-element character vector
  if (!is_1length_character(id)) {
    warning(paste0("While creating lseurat object.\nInvalid 'id' ",
                   "argument. 'id' should be either a 1-length atomic vector of ",
                   "type 'character' or NA.\nSetting 'id' to NA"))
    id <- NA
  }

  # Combine the components into a single named list
  lseurat_obj <- list(
    "seurat" = seurat_obj,
    "params" = params,
    "gmm" = gmm,
    "loom_path" = loom_path,
    "id" = id
  )

  # Assign class attribute
  class(lseurat_obj) <- "lseurat"

  # Return object
  return(lseurat_obj)
}

#' Function for calling cells using gaussian mixture models
#'
#' @param counts When 'counts' is an atomic integer vector, the values should
#'  represent total read counts per cell (e.g. cDNA or gDNA counts). A matrix or
#'  dataframe can also be given, in which case the rows correspond with cells
#'  and the columns represent separate modalities (e.g. RNA and DamID). Thus,
#'  cell class membership can be determined based on multiple modalities.
#' @param k An atomic vector of type integer. A gmm is estimated for each
#'  value of 'k', where this value represents the number of components in
#'  the gmm. The cell classifications will be based on the model with the lowest
#'  BIC. To prevent testing of multiple components and to select this yourself,
#'  provide a single integer. Default 'k' is '1:6'.
#' @param log_convert Logical value indicating if the counts should be
#'  log-transformed before model estimation.
#' @param verbose Logical value indicating if the function should message the
#'  chosen number of components and the corresponding bic.
#'
#' @returns A two-element list. The first element is an atomic integer vector,
#'  where the integers represent class membership to one of the predicted
#'  classes (e.g. 1: 'no cell', 2: 'real cell'). The second element represents
#'  the 'Mclust' object that was used for the class membership prediction.
#'
#' @export
#'
#' @importFrom mclust mclustBIC
gmm_cell_caller <- function(counts, k = 1:6,
                            log_convert = TRUE,

                            verbose = TRUE) {
  if (log_convert) {
    # Convert counts to log scale
    counts = log1p(counts)
  }
  gmm_model <- mclust::Mclust(counts, G = k)

  if (verbose) {
    message(sprintf("Provided %s observations for cell calling",
                    as.character(gmm_model$n)))
    message(sprintf("optimal gmm model selected had %s components",
                    as.character(gmm_model$G)))
    message(sprintf("BIC optimal model: %s", as.character(gmm_model$bic)))
  }
  return(list(mclust::predict.Mclust(gmm_model)$classification, gmm_model))
}

#' Add Odd-barcode–derived metadata to a Seurat object
#'
#' sciT datasets often consist of cells from different treatment conditions.
#' Using Odd barcodes embedded in the full cell barcode, cells can be stratified
#' by treatment or condition. This function extracts the Odd barcode from the
#' full barcode stored in the Seurat metadata, matches it to a user-provided
#' metadata table, and adds both the Odd barcode and the corresponding condition
#' as new metadata columns to the Seurat object contained within an \code{lseurat}
#' object.
#'
#' @param lseurat An object of class \code{lseurat} containing a Seurat object in
#'   \code{lseurat$seurat}.
#' @param Odd_barcode_md_file A comma-separated file linking Odd barcode IDs to
#'   metadata (e.g. treatment or condition). Must contain at least two columns:
#'   one with Odd barcode identifiers and one with the corresponding condition.
#' @param full_barcode_col Character string specifying the column name in the
#'   Seurat metadata that contains the full cell barcode.
#'   Defaults to \code{"BC"}.
#' @param condition_col Character string specifying the name of the metadata
#'   column to store the condition labels in the Seurat object.
#'   Defaults to \code{"condition"}.
#' @param Odd_barcode_col Character string specifying the name of the metadata
#'   column to store the extracted Odd barcodes in the Seurat object.
#'   Defaults to \code{"Odd_barcode"}.
#' @param Odd_barcode_position Integer specifying the position of the Odd barcode
#'   after splitting the full barcode on underscores (\code{"_"}).
#'   Defaults to \code{2}.
#'
#' @details
#' The full barcode is split using underscores, and the element at
#' \code{Odd_barcode_position} is interpreted as the Odd barcode. These Odd
#' barcodes are then matched to the metadata file using exact matching.
#'
#' @returns
#' The input \code{lseurat} object with updated Seurat metadata. Two new columns
#' are added to \code{lseurat$seurat@meta.data}: one containing the extracted Odd
#' barcodes and one containing the corresponding condition labels.
#'
#' @export
Add_odd_barcode_metadata <- function(lseurat,
                                     Odd_barcode_md_file,
                                     full_barcode_col = "BC",
                                     condition_col = "condition",
                                     Odd_barcode_col = "Odd_barcode",
                                     Odd_barcode_position = 2) {

  # Input checks

  if (!is.list(lseurat) || is.null(lseurat$seurat)) {
    stop("'lseurat' must be an object containing a Seurat object in lseurat$seurat.")
  }

  if (!inherits(lseurat$seurat, "Seurat")) {
    stop("lseurat$seurat must be a valid Seurat object.")
  }

  if (!is.character(Odd_barcode_md_file) || length(Odd_barcode_md_file) != 1) {
    stop("'Odd_barcode_md_file' must be a single character string.")
  }

  if (!file.exists(Odd_barcode_md_file)) {
    stop("Odd barcode metadata file does not exist: ", Odd_barcode_md_file)
  }

  if (!is.numeric(Odd_barcode_position) || length(Odd_barcode_position) != 1 ||
      Odd_barcode_position < 1) {
    stop("'Odd_barcode_position' must be a single positive integer.")
  }

  # Load data

  Odd_barcode_md <- utils::read.csv(Odd_barcode_md_file)

  required_cols <- c("Odd_barcode", "Condition")
  missing_cols <- setdiff(required_cols, colnames(Odd_barcode_md))

  if (length(missing_cols) > 0) {
    stop(
      "Odd barcode metadata file is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  seurat_md <- lseurat$seurat@meta.data

  if (!full_barcode_col %in% colnames(seurat_md)) {
    stop("Column '", full_barcode_col, "' not found in Seurat metadata.")
  }

  # Extract Odd barcodes

  Odd_bcs <- sapply(
    strsplit(seurat_md[[full_barcode_col]], split = "_"),
    "[",
    Odd_barcode_position
  )

  # Map to condition

  condition <- Odd_barcode_md$Condition[
    match(Odd_bcs, Odd_barcode_md$Odd_barcode)
  ]

  # Warn if some cells could not be matched to a condition
  if (any(is.na(condition))) {
    warning(sum(is.na(condition)), " cells could not be matched to a condition.")
  }

  # Add metadata

  seurat_md[[Odd_barcode_col]] <- Odd_bcs
  seurat_md[[condition_col]] <- condition

  lseurat$seurat@meta.data <- seurat_md

  return(lseurat)
}

#' Convert a loom file to a Seurat object
#'
#' This function loads a transcriptome loom file and converts it into a
#' Seurat object wrapped inside an \code{lseurat} container. It supports
#' optional Gaussian mixture model (GMM)–based cell calling, automatic
#' resolution of duplicate feature names, and annotation of cells using
#' Odd-barcode–linked metadata.
#'
#' @section Loom file input:
#' Arguments: \code{loom_path}, \code{matrix_rowname_col},
#' \code{matrix_colname_col}.
#'
#' Control how the loom file is read and how features and cells are identified.
#'
#' @param loom_path A single-element character vector specifying the path to
#'   the transcriptome loom file.
#' @param matrix_rowname_col A single-element character vector indicating the
#'   column in the row (gene) metadata of the loom file that should be used as
#'   feature names for the count matrix.
#' @param matrix_colname_col A single-element character vector indicating the
#'   column in the column (cell) metadata of the loom file that should be used
#'   as cell names for the count matrix.
#'
#' @section Seurat object construction:
#' Arguments: \code{seurat_assay_name}, \code{seurat_min_cells},
#' \code{seurat_names_field}, \code{seurat_names_delim}, \code{...}.
#'
#' Control creation, filtering, and naming of the Seurat object.
#'
#' @param seurat_assay_name A single-element character vector specifying the
#'   name of the assay in the created Seurat object.
#' @param seurat_min_cells Include only features detected in at least this many
#'   cells.
#' @param seurat_names_field Integer specifying which substring of the cell
#'   name (after splitting by \code{seurat_names_delim}) is used as the sample
#'   name.
#' @param seurat_names_delim A single-element character vector specifying the
#'   delimiter used to split cell names when extracting sample names.
#' @param ... Additional arguments passed to
#'   \code{SeuratObject::CreateSeuratObject()}.
#'
#' @section Duplicate feature handling:
#' Arguments: \code{resolve_duplicates}.
#'
#' Control how duplicate feature names are detected and resolved.
#'
#' @param resolve_duplicates Logical value indicating whether duplicate
#'   feature names should be resolved automatically. If \code{TRUE}, only the
#'   first occurrence of each duplicated feature is retained. If \code{FALSE},
#'   the function errors when duplicates are detected. Default is \code{FALSE}.
#'
#' @section GMM-based cell calling:
#' Arguments: \code{gmm_cell_calling}, \code{gmm_k}, \code{gmm_verbose}.
#'
#' Control Gaussian mixture model–based classification of cells.
#'
#' @param gmm_cell_calling Logical value indicating whether GMM-based cell
#'   calling should be performed. Default is \code{FALSE}.
#' @param gmm_k Integer vector specifying the number(s) of mixture components
#'   to test. The model with the lowest BIC is selected. Default is \code{1:6}.
#' @param gmm_verbose Logical value indicating whether information about the
#'   selected GMM (BIC and number of components) should be printed.
#'
#' @section Odd-barcode metadata annotation:
#' Arguments: \code{Add_Odd_bc_md}, \code{Odd_barcode_md_file},
#' \code{full_barcode_col}, \code{condition_col},
#' \code{Odd_barcode_col}, \code{Odd_barcode_position}.
#'
#' Control annotation of cells using Odd-barcode–linked metadata.
#'
#' @param Add_Odd_bc_md Logical value indicating whether Odd-barcode metadata
#'   annotation should be performed. Default is \code{FALSE}.
#' @param Odd_barcode_md_file A comma-separated file linking Odd barcode IDs to
#'   metadata (e.g. treatment or condition). Must contain at least two columns:
#'   one with Odd barcode identifiers and one with the corresponding condition.
#' @param full_barcode_col Character string specifying the column in the Seurat
#'   metadata that contains the full cell barcode.
#' @param condition_col Character string specifying the name of the metadata
#'   column in which to store condition labels.
#' @param Odd_barcode_col Character string specifying the name of the metadata
#'   column in which to store extracted Odd barcodes.
#' @param Odd_barcode_position Integer specifying the position of the Odd
#'   barcode after splitting the full barcode on underscores (\code{"_"}).
#'
#' @returns
#' An object of class \code{lseurat}. This is a named list containing:
#' \describe{
#'   \item{seurat}{A Seurat object containing the count matrix and metadata
#'   derived from the loom file.}
#'   \item{params}{A named list of parameters supplied to
#'   \code{LoomAsSeurat()}.}
#'   \item{gmm}{The \code{Mclust} object used for cell calling, or \code{NA} if
#'   GMM-based cell calling was not performed.}
#' }
#'
#' @export
#'
#' @examples
#' # Example 1: basic loom loading
#' lpath <- system.file("extdata", "i1_ES_subsample.loom",
#'   package = "LoadSciLooms")
#' lseurat <- LoomAsSeurat(
#'   lpath,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   resolve_duplicates = FALSE,
#'   seurat_assay_name = "RNA",
#'   seurat_min_cells = 0,
#'   seurat_names_field = 1L,
#'   seurat_names_delim = "_"
#' )
#'
#' # Example 2: loom loading with Odd-barcode metadata annotation
#' loom_file <- system.file("extdata", "i31_GSK126_subsample.loom",
#'   package = "LoadSciLooms")
#' odd_md_file <- system.file("extdata", "Odd_barcode_md.csv",
#'   package = "LoadSciLooms")
#' lseurat_odd <- LoomAsSeurat(
#'   loom_path = loom_file,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   resolve_duplicates = FALSE,
#'   Add_Odd_bc_md = TRUE,
#'   Odd_barcode_md_file = odd_md_file,
#'   full_barcode_col = "BC"
#' )
#'
#' # Resulting Seurat metadata will contain Odd-barcode metadata
#' head(lseurat_odd$seurat[[c("Odd_barcode", "condition")]])
LoomAsSeurat <- function(loom_path,
                         matrix_rowname_col = "Gene",
                         matrix_colname_col = "CellID",
                         gmm_cell_calling = FALSE,
                         gmm_k = 1:6,
                         gmm_verbose = TRUE,
                         resolve_duplicates = FALSE,
                         Add_Odd_bc_md = FALSE,
                         Odd_barcode_md_file = NULL,
                         full_barcode_col = "BC",
                         condition_col = "condition",
                         Odd_barcode_col = "Odd_barcode",
                         Odd_barcode_position = 2,
                         seurat_assay_name = "RNA",
                         seurat_min_cells = 0,
                         seurat_names_field = 1L,
                         seurat_names_delim = "_",
                         ...) {

  # Capture all passed arguments
  all_args <- c(as.list(environment()), list(...))

  # Capture loom_path's value immediately to avoid promise issues in handlers
  loom_path_val <- loom_path

  ### Check types
  if (!((is.logical(gmm_cell_calling)) & (length(gmm_cell_calling) == 1))) {
    stop("gmm_cell_calling is not of type logical or is not of length one.")
  }

  if (!((is.logical(resolve_duplicates)) & (length(resolve_duplicates) == 1))) {
    stop("resolve_duplicates is not of type logical or is not of length one.")
  }

  if (!((is.logical(Add_Odd_bc_md)) & (length(Add_Odd_bc_md) == 1))) {
    stop("Add_Odd_bc_md is not of type logical or is not of length one.")
  }

  ### Connect to loom file
  lfile <- withCallingHandlers(
    tryCatch(rhdf5::H5Fopen(loom_path_val, flags = "H5F_ACC_RDONLY"),
             error = function(e) {
               stop(sprintf("Error opening loom file: %s\nOriginal error message:\n%s",
                            loom_path_val, conditionMessage(e)), call. = FALSE)
             }),
    warning = function(w) {
      message(sprintf("Warning while opening loom file: %s", loom_path_val))
      message("Original warning message:\n")
      message(conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Add a check for successful connection
  if (!inherits(lfile, "H5IdComponent")) {
    stop(sprintf("Failed to connect to loom file at: %s.\nfile connection is not of type 'hdf5'",
                 loom_path_val), call. = FALSE)
  }

  # Close loom files upon exit of function
  on.exit(rhdf5::h5closeAll(), add = TRUE)

  ### Perform checks
  if (is.null(lfile$row_attrs[[matrix_rowname_col]])) {
    stop(sprintf("'%s' is not an existing column in the row metadata of the loom file at path '%s'", matrix_rowname_col, loom_path_val))
  }

  if (is.null(lfile$col_attrs[[matrix_colname_col]])) {
    stop(sprintf("'%s' is not an existing column in the column metadata of the loom file at path '%s'", matrix_colname_col, loom_path_val))
  }

  mtx_dims <- dim(lfile$matrix)

  if (mtx_dims[2] != length(lfile$row_attrs[[matrix_rowname_col]])) {
    stop(sprintf("The row number for the matrix and the row number for the row-metadata are not equal for the loom file at path '%s'.", loom_path_val))
  }

  if (mtx_dims[1] != length(lfile$col_attrs[[matrix_colname_col]])) {
    stop(sprintf("The column number for the matrix and the row number for the column-metadata are not equal for the loom file at path '%s'.", loom_path_val))
  }

  ### Collect the metadata
  row_meta <- as.data.frame(lfile$row_attrs)
  col_meta <- as.data.frame(lfile$col_attrs)

  ### Convert to sparse matrix. Transposing last is more memory efficient.
  lmatrix_sparse <- Matrix::t(Matrix::Matrix(lfile$matrix, sparse = TRUE))

  ### Add row and column names to sparse matrix
  rownames(lmatrix_sparse) <- row_meta[[matrix_rowname_col]]
  colnames(lmatrix_sparse) <- col_meta[[matrix_colname_col]]

  ### Create the seurat object, keeping in mind potential duplicate features

  # Create index pointing to duplicate rownames
  duplicate_index <- (duplicated(rownames(lmatrix_sparse)) == FALSE)

  # Check if there are any duplicates
  if (any(duplicate_index == FALSE)) {
    if (resolve_duplicates == FALSE) {
      stop(sprintf(paste0("The feature column '%s' in '%s' contains duplicate values.\n",
                          "This column is used to set the feature names in the resulting Seurat object.\n",
                          "To fix this, either set 'resolve_duplicates = TRUE' to automatically resolve duplicates,\n",
                          "or provide a different column via the 'matrix_rowname_col' argument so that feature names are unique."),
                   matrix_rowname_col,
                   loom_path_val
      ))
    } else {
      # Create seurat subsetted to remove duplicates
      lseurat <- SeuratObject::CreateSeuratObject(counts = lmatrix_sparse[duplicate_index, ],
                                                  assay = seurat_assay_name,
                                                  min.cells = seurat_min_cells,
                                                  names.field = seurat_names_field,
                                                  names.delim = seurat_names_delim,
                                                  ...)
      # Remove duplicates from row_meta
      row_meta <- row_meta[duplicate_index, ]
    }
  } else {
    # No duplicates, so create the seurat object without subsetting
    lseurat <- SeuratObject::CreateSeuratObject(counts = lmatrix_sparse,
                                                assay = seurat_assay_name,
                                                min.cells = seurat_min_cells,
                                                names.field = seurat_names_field,
                                                names.delim = seurat_names_delim,
                                                ...)
  }

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

  ### Perform gmm cell calling and return lseurat object
  if (gmm_cell_calling) {
    gmm_calls <- gmm_cell_caller(lseurat$nCount_RNA, k = gmm_k,
                                 log_convert = TRUE, verbose = gmm_verbose)
    # Add gmm cell calls to seurat metadata
    lseurat$gmm <- gmm_calls[[1]]

    # Create lseurat object with gmm object
    lseurat <- create_lseurat_object(lseurat, all_args, gmm_calls[[2]])
  } else {
    # Create lseurat object without gmm object
    lseurat <- create_lseurat_object(lseurat, all_args, NA)
  }

  ### Perform Odd-barcode linked metadata annotation
  if (Add_Odd_bc_md) {
    lseurat <- Add_odd_barcode_metadata(lseurat = lseurat,
                                        Odd_barcode_md_file = Odd_barcode_md_file,
                                        full_barcode_col = full_barcode_col,
                                        condition_col = condition_col,
                                        Odd_barcode_col = Odd_barcode_col,
                                        Odd_barcode_position = Odd_barcode_position)
  }
  return(lseurat)
}

#' Convert multiple loom files to a merged Seurat object
#'
#' This function converts multiple transcriptome loom files into individual
#' Seurat objects using \code{\link{LoomAsSeurat}} and optionally merges them
#' into a single Seurat object. Each loom file is tracked by a user-supplied
#' identifier, which is stored in the cell-level metadata of the resulting
#' Seurat object(s).
#'
#' @section Loom file input:
#' Arguments: \code{loom_paths}, \code{matrix_rowname_col},
#' \code{matrix_colname_col}.
#'
#' Control how multiple loom files are specified and how features and cells
#' are identified.
#'
#' @param loom_paths A named character vector specifying paths to transcriptome
#'   loom files. The names are used as sample identifiers and will be stored in
#'   the cell metadata (column \code{loom}) of the resulting Seurat object(s).
#' @param matrix_rowname_col A single-element character vector indicating the
#'   column in the row (gene) metadata of the loom files that should be used as
#'   feature names for the count matrix.
#' @param matrix_colname_col A single-element character vector indicating the
#'   column in the column (cell) metadata of the loom files that should be used
#'   as cell names for the count matrix.
#'
#' @section Seurat object construction:
#' Arguments: \code{seurat_assay_name}, \code{seurat_min_cells},
#' \code{seurat_names_field}, \code{seurat_names_delim}, \code{...}.
#'
#' Control creation, filtering, and naming of Seurat objects derived from each
#' loom file.
#'
#' @param seurat_assay_name A single-element character vector specifying the
#'   name of the assay in the created Seurat object(s).
#' @param seurat_min_cells Include only features detected in at least this many
#'   cells.
#' @param seurat_names_field Integer specifying which substring of the cell
#'   name (after splitting by \code{seurat_names_delim}) is used as the sample
#'   name.
#' @param seurat_names_delim A single-element character vector specifying the
#'   delimiter used to split cell names when extracting sample names.
#' @param ... Additional arguments passed to
#'   \code{SeuratObject::CreateSeuratObject()}.
#'
#' @section Duplicate feature handling:
#' Arguments: \code{resolve_duplicates}.
#'
#' Control how duplicate feature names are detected and resolved.
#'
#' @param resolve_duplicates Logical value indicating whether duplicate
#'   feature names should be resolved automatically. If \code{TRUE}, only the
#'   first occurrence of each duplicated feature is retained. If \code{FALSE},
#'   the function errors when duplicates are detected. Default is \code{FALSE}.
#'
#' @section GMM-based cell calling:
#' Arguments: \code{gmm_cell_calling}, \code{gmm_k}, \code{gmm_verbose}.
#'
#' Control Gaussian mixture model–based classification of cells for each loom
#' file.
#'
#' @param gmm_cell_calling Logical value indicating whether GMM-based cell
#'   calling should be performed for each loom file. Default is \code{FALSE}.
#' @param gmm_k Integer vector specifying the number(s) of mixture components
#'   to test. The model with the lowest BIC is selected. Default is \code{1:6}.
#' @param gmm_verbose Logical value indicating whether information about the
#'   selected GMM (BIC and number of components) should be printed.
#'
#' @section Odd-barcode metadata annotation:
#' Arguments: \code{Add_Odd_bc_md}, \code{Odd_barcode_md_file},
#' \code{full_barcode_col}, \code{condition_col},
#' \code{Odd_barcode_col}, \code{Odd_barcode_position}.
#'
#' Control annotation of cells using Odd-barcode–linked metadata.
#'
#' @param Add_Odd_bc_md Logical value indicating whether Odd-barcode metadata
#'   annotation should be performed. Default is \code{FALSE}.
#' @param Odd_barcode_md_file A comma-separated file linking Odd barcode IDs to
#'   metadata (e.g. treatment or condition). Must contain at least two columns:
#'   one with Odd barcode identifiers and one with the corresponding condition.
#' @param full_barcode_col Character string specifying the column in the Seurat
#'   metadata that contains the full cell barcode.
#' @param condition_col Character string specifying the name of the metadata
#'   column in which to store condition labels.
#' @param Odd_barcode_col Character string specifying the name of the metadata
#'   column in which to store extracted Odd barcodes.
#' @param Odd_barcode_position Integer specifying the position of the Odd
#'   barcode after splitting the full barcode on underscores (\code{"_"}).
#'
#' @section Merging behavior:
#' Arguments: \code{merge_seurat}, \code{remove_unmerged}.
#'
#' Control whether individual Seurat objects are merged and whether unmerged
#' objects are retained.
#'
#' @param merge_seurat Logical value indicating whether individual Seurat
#'   objects should be merged into a single Seurat object.
#' @param remove_unmerged Logical value indicating whether the individual
#'   Seurat objects should be removed after merging. This argument is ignored
#'   if \code{merge_seurat = FALSE}.
#'
#' @returns
#' A named list containing:
#' \describe{
#'   \item{merged_seurat}{A merged Seurat object if
#'   \code{merge_seurat = TRUE}, otherwise \code{NA}.}
#'   \item{<sample_id>}{One element per loom file, each an object of class
#'   \code{lseurat} containing the Seurat object, parameters, and optional GMM
#'   results.}
#' }
#'
#' @seealso
#' \code{\link{LoomAsSeurat}}, \code{\link{Add_odd_barcode_metadata}}
#'
#' @export
#'
#' @examples
#' ## Example 1: Regular loading and merging of multiple loom files
#' lpaths <- c(
#'   "i1_ES"  = system.file("extdata", "i1_ES_subsample.loom",
#'                          package = "LoadSciLooms"),
#'   "i3_d6"  = system.file("extdata", "i3_d6_subsample.loom",
#'                          package = "LoadSciLooms"),
#'   "i5_d17" = system.file("extdata", "i5_d17_subsample.loom",
#'                          package = "LoadSciLooms")
#' )
#' lseurat_list <- MultiLoomAsSeurat(
#'   lpaths,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   seurat_assay_name = "RNA",
#'   seurat_min_cells = 0,
#'   seurat_names_field = 1L,
#'   seurat_names_delim = "_"
#' )
#'
#' ## Example 2: loom loading with Odd-barcode metadata annotation
#' lpaths <- c(
#'   "dummy_sample1" = system.file("extdata", "i31_GSK126_subsample.loom",
#'                                 package = "LoadSciLooms"),
#'   "dummy_sample2" = system.file("extdata", "i31_GSK126_subsample.loom",
#'                                 package = "LoadSciLooms"),
#'   "dummy_sample3" = system.file("extdata", "i31_GSK126_subsample.loom",
#'                                 package = "LoadSciLooms")
#' )
#' odd_md_file <- system.file("extdata", "Odd_barcode_md.csv",
#'   package = "LoadSciLooms")
#' lseurat_list_odd <- MultiLoomAsSeurat(
#'   lpaths,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   Add_Odd_bc_md = TRUE,
#'   Odd_barcode_md_file = odd_md_file,
#'   remove_unmerged = TRUE,
#'   full_barcode_col = "BC"
#'   )
#'
#' # Resulting merged-Seurat metadata will contain Odd-barcode metadata
#' head(lseurat_list_odd$merged_seurat[[c("Odd_barcode", "condition")]])
MultiLoomAsSeurat <- function(loom_paths,
                              matrix_rowname_col = "Gene",
                              matrix_colname_col = "CellID",
                              gmm_cell_calling = FALSE,
                              gmm_k = 1:6,
                              gmm_verbose = TRUE,
                              merge_seurat = TRUE,
                              remove_unmerged = TRUE,
                              resolve_duplicates = FALSE,
                              Add_Odd_bc_md = FALSE,
                              Odd_barcode_md_file = NULL,
                              full_barcode_col = "BC",
                              condition_col = "condition",
                              Odd_barcode_col = "Odd_barcode",
                              Odd_barcode_position = 2,
                              seurat_assay_name = "RNA",
                              seurat_min_cells = 0,
                              seurat_names_field = 1L,
                              seurat_names_delim = "_",
                              ...) {

  # Function to check if it's a named atomic vector
  is_named_atomic_vector <- function(x) {
    is.atomic(x) && !is.null(names(x)) && length(x) > 0
  }

  # Function to check if argument is logical of length one
  is_1length_logical <- function(x) {
    return((is.logical(x)) & (length(x) == 1))
  }

  # Check if loom_paths is named atomic vector
  if (!(is_named_atomic_vector(loom_paths))) {
    stop("Loom_paths is not a named vector")
  }

  # Check if arguments are logical
  if (!is_1length_logical(merge_seurat)) {
    stop("'merge_seurat' is not of type logical or is not of length one.")
  }

  if (!is_1length_logical(remove_unmerged)) {
    stop("'remove_unmerged' is not of type logical or is not of length one.")
  }

  path_names <- names(loom_paths)

  # Check if all paths have non-empty names
  if (any(path_names == "")) {
    idx_noname <- which(path_names == "")
    path_noname <- paste(loom_paths[idx_noname], collapse = "\n")
    stop(sprintf(paste0("The following paths in the 'loom_paths' vector have no ",
                        "identifiers (i.e. names):\n%s"), path_noname))
  }

  # Check if names are unique
  if (any(duplicated(path_names))) {
    idx_notuniq <- which(path_names %in% path_names[anyDuplicated(path_names)])
    path_notuniq <- paste(loom_paths[idx_notuniq], collapse = "\n")
    stop(sprintf(paste0("The following paths in the 'loom_paths' vector have non-unique ",
                        "identifiers (i.e. names):\n%s"), path_notuniq))
  }

  # Check if file paths exist
  if (any(file.exists(loom_paths) == FALSE)) {
    non_existing_idx <- which(file.exists(loom_paths) == FALSE)
    stop(sprintf("The following files do not exist:\n%s", paste(loom_paths[non_existing_idx], collapse = "\n")))
  }

  lseurat_obj_list <- lapply(seq_along(loom_paths), function(lp_idx) {
    message(sprintf("Converting %s to seurat", names(loom_paths)[lp_idx]))
    lseurat <- LoomAsSeurat(loom_path = loom_paths[lp_idx],
                 matrix_rowname_col = matrix_rowname_col,
                 matrix_colname_col = matrix_colname_col,
                 gmm_cell_calling = gmm_cell_calling,
                 gmm_k = gmm_k,
                 gmm_verbose = gmm_verbose,
                 resolve_duplicates = resolve_duplicates,
                 Add_Odd_bc_md = Add_Odd_bc_md,
                 Odd_barcode_md_file = Odd_barcode_md_file,
                 full_barcode_col = full_barcode_col,
                 condition_col = condition_col,
                 Odd_barcode_col = Odd_barcode_col,
                 Odd_barcode_position = Odd_barcode_position,
                 seurat_assay_name = seurat_assay_name,
                 seurat_min_cells = seurat_min_cells,
                 seurat_names_field = seurat_names_field,
                 seurat_names_delim = seurat_names_delim,
                 ...)
    # Add column with loom path name to metadata
    lseurat$seurat$loom <- names(loom_paths)[lp_idx]

    # Add additional data to object
    lseurat$loom_path <- loom_paths[lp_idx]
    lseurat$id <- names(loom_paths)[lp_idx]

    return(lseurat)
  })
  names(lseurat_obj_list) <- path_names

  # Merge seurat objects if option is set
  if (merge_seurat) {
    # Retrieve seuratobjects from lseurat object list
    loom_seurat_list <- lapply(lseurat_obj_list, `[[`, "seurat")

    if (length(loom_seurat_list) > 1) {
      merged_seurat <- merge(x = loom_seurat_list[[1]], y = loom_seurat_list[-1])
    } else if (length(loom_seurat_list) == 1) {
      merged_seurat <- loom_seurat_list[[1]]
    } else {
      stop("Loom_seurat_list is empty, meaning that the LoomAsSeurat function failed to produce any seurat object")
    }

    # Add merged seurat to lseurat object list
    lseurat_obj_list$merged_seurat <- merged_seurat

    # Remove the individual seurat objects if specified
    if (remove_unmerged) {
      lseurat_obj_list <- lapply(lseurat_obj_list, function(i) {
        if (inherits(i, "lseurat")) {
          i$seurat <- NA
        }
        return(i)
      })
    }
  } else {
    lseurat_obj_list$merged_seurat <- NA
  }
  return(lseurat_obj_list)
}
