# Convert multiple loom files to a single merged Seuratobject

Convert multiple loom files to a single merged Seuratobject

## Usage

``` r
MultiLoomAsSeurat(
  loom_paths,
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
  ...
)
```

## Arguments

- loom_paths:

  A named character vector that indicates the paths to the transcriptome
  loom files. The Name attribute will be used to indicate which cell
  came from which loom file.

- matrix_rowname_col:

  A single element character vector that indicates the variable (i.e.
  column) of the row (i.e. gene) metadata in the loom file. This
  variable will subsequently be used as the feature names for the
  generated sparse matrix.

- matrix_colname_col:

  A single element character vector that indicates the variable (i.e.
  column) of the column (i.e. cell) metadata in the loom file. This
  variable will subsequently be used as the cell names for the generated
  sparse matrix.

- gmm_cell_calling:

  A logical value indicating if a gaussian mixture model (gmm) is used
  to predict which cells are likely real cell. If set to TRUE, The
  'mclust' R package is used for this purpose. Default is FALSE.

- gmm_k:

  An atomic vector of type integer. A gmm is estimated for each value of
  'gmm_k', where this value represents the number of components in the
  gmm. The cell classifications will be based on the model with the
  lowest BIC. Default 'gmm_k' is '1:6'

- gmm_verbose:

  A logical value. For TRUE, will report the BIC and the number of
  components for the chosen gmm. for FALSE, no messages will be printed
  about the chosen model.

- merge_seurat:

  A logical value. If TRUE, the individual seuratobjects for each loom
  file will be merged. The merged seuratobject will then be added to the
  returned list.

- remove_unmerged:

  A logical value. If TRUE, the individual seuratobjects for each loom
  file will be removed after merging. If FALSE, the individual
  seuratobjects will be kept in the returned list, even after merging.
  This parameter is not used when 'merge_seurat' is set to FALSE.

- seurat_assay_name:

  A single element character vector representing the name of the assay
  in the Seuratobject.

- seurat_min_cells:

  Include only features in the seuratobject that are detected in at
  least this many cells.

- seurat_names_field:

  Argument used to get sample name from the cell names. Works in
  conjunction with the 'seurat_names_delim' argument.
  'seurat_names_field' is an integer selecting which substring (i.e.
  delimited by 'seurat_names_delim') should be selected as the sample
  name. The default is 1L (i.e. the first field of the cell name).

- seurat_names_delim:

  Argument used to get sample name from the cell names. Works in
  conjunction with the 'seurat_names_field' argument.
  'seurat_names_delim' is a single element character vector that
  represents the delimiter used to select the cell name substrings.

- ...:

  Additional arguments will be given to the
  SeuratObject::CreateSeuratObject() function call.

## Value

A list-like object containing a merged seuratobject and objects of type
'lseurat' that were produced by the LoomasSeurat function.

A Seuratobject populated with a sparse matrix representing the
transcriptome data from the loom files located at 'loom_paths'. The
names registered in the 'loom_paths' variable are used to identify from
which sample each cell came from in the merged object. This information
can be found in the 'loom' column of the seurat cell metadata.

## Examples

``` r
# Get named vector of loom file paths
lpaths <- c("i1_ES" = system.file("extdata", "i1_ES_subsample.loom",
                                   package = "LoadSciLooms"),
            "i3_d6" = system.file("extdata", "i3_d6_subsample.loom",
                                   package = "LoadSciLooms"),
            "i5_d17" = system.file("extdata", "i5_d17_subsample.loom",
                                    package = "LoadSciLooms"))

# Load loom files as seurat
lseurat <- MultiLoomAsSeurat(lpaths,
  matrix_rowname_col = "Gene",
  matrix_colname_col = "CellID",
  seurat_assay_name = "RNA",
  seurat_min_cells = 0,
  seurat_names_field = 1L,
  seurat_names_delim = "_")
#> Converting i1_ES to seurat
#> Converting i3_d6 to seurat
#> Converting i5_d17 to seurat
```
