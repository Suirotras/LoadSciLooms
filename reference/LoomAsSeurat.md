# Convert a loom file to a Seurat object

This function loads a transcriptome loom file and converts it into a
Seurat object wrapped inside an `lseurat` container. It supports
optional Gaussian mixture model (GMM)–based cell calling, automatic
resolution of duplicate feature names, and annotation of cells using
Odd-barcode–linked metadata.

## Usage

``` r
LoomAsSeurat(
  loom_path,
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
  ...
)
```

## Arguments

- loom_path:

  A single-element character vector specifying the path to the
  transcriptome loom file.

- matrix_rowname_col:

  A single-element character vector indicating the column in the row
  (gene) metadata of the loom file that should be used as feature names
  for the count matrix.

- matrix_colname_col:

  A single-element character vector indicating the column in the column
  (cell) metadata of the loom file that should be used as cell names for
  the count matrix.

- gmm_cell_calling:

  Logical value indicating whether GMM-based cell calling should be
  performed. Default is `FALSE`.

- gmm_k:

  Integer vector specifying the number(s) of mixture components to test.
  The model with the lowest BIC is selected. Default is `1:6`.

- gmm_verbose:

  Logical value indicating whether information about the selected GMM
  (BIC and number of components) should be printed.

- resolve_duplicates:

  Logical value indicating whether duplicate feature names should be
  resolved automatically. If `TRUE`, only the first occurrence of each
  duplicated feature is retained. If `FALSE`, the function errors when
  duplicates are detected. Default is `FALSE`.

- Add_Odd_bc_md:

  Logical value indicating whether Odd-barcode metadata annotation
  should be performed. Default is `FALSE`.

- Odd_barcode_md_file:

  A comma-separated file linking Odd barcode IDs to metadata (e.g.
  treatment or condition). Must contain at least two columns: one with
  Odd barcode identifiers and one with the corresponding condition.

- full_barcode_col:

  Character string specifying the column in the Seurat metadata that
  contains the full cell barcode.

- condition_col:

  Character string specifying the name of the metadata column in which
  to store condition labels.

- Odd_barcode_col:

  Character string specifying the name of the metadata column in which
  to store extracted Odd barcodes.

- Odd_barcode_position:

  Integer specifying the position of the Odd barcode after splitting the
  full barcode on underscores (`"_"`).

- seurat_assay_name:

  A single-element character vector specifying the name of the assay in
  the created Seurat object.

- seurat_min_cells:

  Include only features detected in at least this many cells.

- seurat_names_field:

  Integer specifying which substring of the cell name (after splitting
  by `seurat_names_delim`) is used as the sample name.

- seurat_names_delim:

  A single-element character vector specifying the delimiter used to
  split cell names when extracting sample names.

- ...:

  Additional arguments passed to
  [`SeuratObject::CreateSeuratObject()`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html).

## Value

An object of class `lseurat`. This is a named list containing:

- seurat:

  A Seurat object containing the count matrix and metadata derived from
  the loom file.

- params:

  A named list of parameters supplied to `LoomAsSeurat()`.

- gmm:

  The `Mclust` object used for cell calling, or `NA` if GMM-based cell
  calling was not performed.

## Loom file input

Arguments: `loom_path`, `matrix_rowname_col`, `matrix_colname_col`.

Control how the loom file is read and how features and cells are
identified.

## Seurat object construction

Arguments: `seurat_assay_name`, `seurat_min_cells`,
`seurat_names_field`, `seurat_names_delim`, `...`.

Control creation, filtering, and naming of the Seurat object.

## Duplicate feature handling

Arguments: `resolve_duplicates`.

Control how duplicate feature names are detected and resolved.

## GMM-based cell calling

Arguments: `gmm_cell_calling`, `gmm_k`, `gmm_verbose`.

Control Gaussian mixture model–based classification of cells.

## Odd-barcode metadata annotation

Arguments: `Add_Odd_bc_md`, `Odd_barcode_md_file`, `full_barcode_col`,
`condition_col`, `Odd_barcode_col`, `Odd_barcode_position`.

Control annotation of cells using Odd-barcode–linked metadata.

## Examples

``` r
# Example 1: basic loom loading
lpath <- system.file("extdata", "i1_ES_subsample.loom",
  package = "LoadSciLooms")
lseurat <- LoomAsSeurat(
  lpath,
  matrix_rowname_col = "Gene",
  matrix_colname_col = "CellID",
  resolve_duplicates = FALSE,
  seurat_assay_name = "RNA",
  seurat_min_cells = 0,
  seurat_names_field = 1L,
  seurat_names_delim = "_"
)

# Example 2: loom loading with Odd-barcode metadata annotation
loom_file <- system.file("extdata", "i31_GSK126_subsample.loom",
  package = "LoadSciLooms")
odd_md_file <- system.file("extdata", "Odd_barcode_md.csv",
  package = "LoadSciLooms")
lseurat_odd <- LoomAsSeurat(
  loom_path = loom_file,
  matrix_rowname_col = "Gene",
  matrix_colname_col = "CellID",
  resolve_duplicates = FALSE,
  Add_Odd_bc_md = TRUE,
  Odd_barcode_md_file = odd_md_file,
  full_barcode_col = "BC"
)

# Resulting Seurat metadata will contain Odd-barcode metadata
head(lseurat_odd$seurat[[c("Odd_barcode", "condition")]])
#>                                  Odd_barcode   condition
#> i31-FV-GSK126-1K_multiome_100859    Odd2Bo60 Treatment_B
#> i31-FV-GSK126-1K_multiome_101271    Odd2Bo88 Treatment_C
#> i31-FV-GSK126-1K_multiome_102301    Odd2Bo62 Treatment_B
#> i31-FV-GSK126-1K_multiome_102545    Odd2Bo18        DMSO
#> i31-FV-GSK126-1K_multiome_103754    Odd2Bo75 Treatment_C
#> i31-FV-GSK126-1K_multiome_104        Odd2Bo9        DMSO
```
