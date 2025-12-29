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
#' class(my_seurat) <- "SeuratObject"
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
#' str(lseurat_obj_full)
#' str(lseurat_obj_basic)
#'
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
  return(list(predict(gmm_model)$classification, gmm_model))
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
#' @param gmm_cell_calling A logical value indicating if a gaussian mixture
#'  model (gmm) is used to predict which cells are likely real cell. If set to
#'  TRUE, The 'mclust' R package is used for this purpose. Default is FALSE.
#' @param gmm_k An atomic vector of type integer. A gmm is estimated for each
#'  value of 'gmm_k', where this value represents the number of components in
#'  the gmm. The cell classifications will be based on the model with the lowest
#'  BIC. To prevent testing of multiple components and to select this yourself,
#'  provide a single integer. Default 'gmm_k' is '1:6'.
#' @param gmm_verbose A logical value. For TRUE, will report the BIC and the
#'  number of components for the chosen gmm. for FALSE, no messages will be
#'  printed about the chosen model.
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
#' @returns An object of class 'lseurat'. This class represents a named list
#'  containing three elements, namely 'seurat', 'params', 'gmm'. 'seurat'
#'  represents the created seuratobject containing the count data
#'  (in sparse matrix format) and metadata from the loom file. 'params' is a
#'  named list storing the parameters given to LoomAsSeurat function. 'gmm' is
#'  the 'Mclust' object which was used to generate the cell class memberships.
#'
#' @export
#'
#' @examples
#' lpath <- system.file("extdata", "i1_ES_subsample.loom",
#'                       package = "LoadSciLooms")
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
                         gmm_cell_calling = FALSE,
                         gmm_k = 1:6,
                         gmm_verbose = TRUE,
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

  ### Perform gmm cell calling and return lseurat object
  if (gmm_cell_calling) {
    gmm_calls <- gmm_cell_caller(lseurat$nCount_RNA, k = gmm_k,
                                 log_convert = TRUE, verbose = gmm_verbose)
    # Add gmm cell calls to seurat metadata
    lseurat$gmm <- gmm_calls[[1]]

    # Create lseurat object with gmm object
    return(create_lseurat_object(lseurat, all_args, gmm_calls[[2]]))
  } else {
    # Create lseurat object without gmm object
    return(create_lseurat_object(lseurat, all_args, NA))
  }
}

#' Convert multiple loom files to a single merged Seuratobject
#'
#' @param loom_paths A named character vector that indicates the paths to the
#'  transcriptome loom files. The Name attribute will be used to indicate which
#'  cell came from which loom file.
#' @param matrix_rowname_col A single element character vector that indicates
#'  the variable (i.e. column) of the row (i.e. gene) metadata in the loom file.
#'  This variable will subsequently be used as the feature names for the
#'  generated sparse matrix.
#' @param matrix_colname_col A single element character vector that indicates
#'  the variable (i.e. column) of the column (i.e. cell) metadata in the loom
#'  file. This variable will subsequently be used as the cell names for the
#'  generated sparse matrix.
#' @param gmm_cell_calling A logical value indicating if a gaussian mixture
#'  model (gmm) is used to predict which cells are likely real cell. If set to
#'  TRUE, The 'mclust' R package is used for this purpose. Default is FALSE.
#' @param gmm_k An atomic vector of type integer. A gmm is estimated for each
#'  value of 'gmm_k', where this value represents the number of components in
#'  the gmm. The cell classifications will be based on the model with the lowest
#'  BIC. Default 'gmm_k' is '1:6'
#' @param gmm_verbose A logical value. For TRUE, will report the BIC and the
#'  number of components for the chosen gmm. for FALSE, no messages will be
#'  printed about the chosen model.
#' @param merge_seurat A logical value. If TRUE, the individual seuratobjects
#'  for each loom file will be merged. The merged seuratobject will then be
#'  added to the returned list.
#' @param remove_unmerged A logical value. If TRUE, the individual seuratobjects
#'  for each loom file will be removed after merging. If FALSE, the individual
#'  seuratobjects will be kept in the returned list, even after merging.
#'  This parameter is not used when 'merge_seurat' is set to FALSE.
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
#' @returns A list-like object containing a merged seuratobject and objects
#'  of type 'lseurat' that were produced by the LoomasSeurat function.
#'
#' A Seuratobject populated with a sparse matrix representing the
#'  transcriptome data from the loom files located at 'loom_paths'.
#'  The names registered in the 'loom_paths' variable are used to identify from
#'  which sample each cell came from in the merged object. This information can
#'  be found in the 'loom' column of the seurat cell metadata.
#'
#' @export
#'
#' @examples
#' # Get named vector of loom file paths
#' lpaths <- c("i1_ES" = system.file("extdata", "i1_ES_subsample.loom",
#'                                    package = "LoadSciLooms"),
#'             "i3_d6" = system.file("extdata", "i3_d6_subsample.loom",
#'                                    package = "LoadSciLooms"),
#'             "i5_d17" = system.file("extdata", "i5_d17_subsample.loom",
#'                                     package = "LoadSciLooms"))
#'
#' # Load loom files as seurat
#' lseurat <- MultiLoomAsSeurat(lpaths,
#'   matrix_rowname_col = "Gene",
#'   matrix_colname_col = "CellID",
#'   seurat_assay_name = "RNA",
#'   seurat_min_cells = 0,
#'   seurat_names_field = 1L,
#'   seurat_names_delim = "_")
MultiLoomAsSeurat <- function(loom_paths,
                              matrix_rowname_col = "Gene",
                              matrix_colname_col = "CellID",
                              gmm_cell_calling = FALSE,
                              gmm_k = 1:6,
                              gmm_verbose = TRUE,
                              merge_seurat = TRUE,
                              remove_unmerged = TRUE,
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
        if (class(i) == "lseurat") {
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
