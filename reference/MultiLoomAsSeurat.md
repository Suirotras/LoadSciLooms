# Convert multiple loom files to a merged Seurat object

This function converts multiple transcriptome loom files into individual
Seurat objects using
[`LoomAsSeurat`](https://suirotras.github.io/LoadSciLooms/reference/LoomAsSeurat.md)
and optionally merges them into a single Seurat object. Each loom file
is tracked by a user-supplied identifier, which is stored in the
cell-level metadata of the resulting Seurat object(s).

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

- loom_paths:

  A named character vector specifying paths to transcriptome loom files.
  The names are used as sample identifiers and will be stored in the
  cell metadata (column `loom`) of the resulting Seurat object(s).

- matrix_rowname_col:

  A single-element character vector indicating the column in the row
  (gene) metadata of the loom files that should be used as feature names
  for the count matrix.

- matrix_colname_col:

  A single-element character vector indicating the column in the column
  (cell) metadata of the loom files that should be used as cell names
  for the count matrix.

- gmm_cell_calling:

  Logical value indicating whether GMM-based cell calling should be
  performed for each loom file. Default is `FALSE`.

- gmm_k:

  Integer vector specifying the number(s) of mixture components to test.
  The model with the lowest BIC is selected. Default is `1:6`.

- gmm_verbose:

  Logical value indicating whether information about the selected GMM
  (BIC and number of components) should be printed.

- merge_seurat:

  Logical value indicating whether individual Seurat objects should be
  merged into a single Seurat object.

- remove_unmerged:

  Logical value indicating whether the individual Seurat objects should
  be removed after merging. This argument is ignored if
  `merge_seurat = FALSE`.

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
  the created Seurat object(s).

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

A named list containing:

- merged_seurat:

  A merged Seurat object if `merge_seurat = TRUE`, otherwise `NA`.

- \<sample_id\>:

  One element per loom file, each an object of class `lseurat`
  containing the Seurat object, parameters, and optional GMM results.

## Loom file input

Arguments: `loom_paths`, `matrix_rowname_col`, `matrix_colname_col`.

Control how multiple loom files are specified and how features and cells
are identified.

## Seurat object construction

Arguments: `seurat_assay_name`, `seurat_min_cells`,
`seurat_names_field`, `seurat_names_delim`, `...`.

Control creation, filtering, and naming of Seurat objects derived from
each loom file.

## Duplicate feature handling

Arguments: `resolve_duplicates`.

Control how duplicate feature names are detected and resolved.

## GMM-based cell calling

Arguments: `gmm_cell_calling`, `gmm_k`, `gmm_verbose`.

Control Gaussian mixture model–based classification of cells for each
loom file.

## Odd-barcode metadata annotation

Arguments: `Add_Odd_bc_md`, `Odd_barcode_md_file`, `full_barcode_col`,
`condition_col`, `Odd_barcode_col`, `Odd_barcode_position`.

Control annotation of cells using Odd-barcode–linked metadata.

## Merging behavior

Arguments: `merge_seurat`, `remove_unmerged`.

Control whether individual Seurat objects are merged and whether
unmerged objects are retained.

## See also

[`LoomAsSeurat`](https://suirotras.github.io/LoadSciLooms/reference/LoomAsSeurat.md),
[`Add_odd_barcode_metadata`](https://suirotras.github.io/LoadSciLooms/reference/Add_odd_barcode_metadata.md)

## Examples

``` r
## Example 1: Regular loading and merging of multiple loom files
lpaths <- c(
  "i1_ES"  = system.file("extdata", "i1_ES_subsample.loom",
                         package = "LoadSciLooms"),
  "i3_d6"  = system.file("extdata", "i3_d6_subsample.loom",
                         package = "LoadSciLooms"),
  "i5_d17" = system.file("extdata", "i5_d17_subsample.loom",
                         package = "LoadSciLooms")
)
lseurat_list <- MultiLoomAsSeurat(
  lpaths,
  matrix_rowname_col = "Gene",
  matrix_colname_col = "CellID",
  seurat_assay_name = "RNA",
  seurat_min_cells = 0,
  seurat_names_field = 1L,
  seurat_names_delim = "_"
)
#> Converting i1_ES to seurat
#> Converting i3_d6 to seurat
#> Converting i5_d17 to seurat

## Example 2: loom loading with Odd-barcode metadata annotation
lpaths <- c(
  "dummy_sample1" = system.file("extdata", "i31_GSK126_subsample.loom",
                                package = "LoadSciLooms"),
  "dummy_sample2" = system.file("extdata", "i31_GSK126_subsample.loom",
                                package = "LoadSciLooms"),
  "dummy_sample3" = system.file("extdata", "i31_GSK126_subsample.loom",
                                package = "LoadSciLooms")
)
odd_md_file <- system.file("extdata", "Odd_barcode_md.csv",
  package = "LoadSciLooms")
lseurat_list_odd <- MultiLoomAsSeurat(
  lpaths,
  matrix_rowname_col = "Gene",
  matrix_colname_col = "CellID",
  Add_Odd_bc_md = TRUE,
  Odd_barcode_md_file = odd_md_file,
  remove_unmerged = TRUE,
  full_barcode_col = "BC"
  )
#> Converting dummy_sample1 to seurat
#> Converting dummy_sample2 to seurat
#> Converting dummy_sample3 to seurat
#> Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.

# Resulting merged-Seurat metadata will contain Odd-barcode metadata
head(lseurat_list_odd$merged_seurat[[c("Odd_barcode", "condition")]])
#>                                    Odd_barcode   condition
#> i31-FV-GSK126-1K_multiome_100859_1    Odd2Bo60 Treatment_B
#> i31-FV-GSK126-1K_multiome_101271_1    Odd2Bo88 Treatment_C
#> i31-FV-GSK126-1K_multiome_102301_1    Odd2Bo62 Treatment_B
#> i31-FV-GSK126-1K_multiome_102545_1    Odd2Bo18        DMSO
#> i31-FV-GSK126-1K_multiome_103754_1    Odd2Bo75 Treatment_C
#> i31-FV-GSK126-1K_multiome_104_1        Odd2Bo9        DMSO
```
