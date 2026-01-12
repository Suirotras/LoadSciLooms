# Add Odd-barcode–derived metadata to a Seurat object

sciT datasets often consist of cells from different treatment
conditions. Using Odd barcodes embedded in the full cell barcode, cells
can be stratified by treatment or condition. This function extracts the
Odd barcode from the full barcode stored in the Seurat metadata, matches
it to a user-provided metadata table, and adds both the Odd barcode and
the corresponding condition as new metadata columns to the Seurat object
contained within an `lseurat` object.

## Usage

``` r
Add_odd_barcode_metadata(
  lseurat,
  Odd_barcode_md_file,
  full_barcode_col = "BC",
  condition_col = "condition",
  Odd_barcode_col = "Odd_barcode",
  Odd_barcode_position = 2
)
```

## Arguments

- lseurat:

  An object of class `lseurat` containing a Seurat object in
  `lseurat$seurat`.

- Odd_barcode_md_file:

  A comma-separated file linking Odd barcode IDs to metadata (e.g.
  treatment or condition). Must contain at least two columns: one with
  Odd barcode identifiers and one with the corresponding condition.

- full_barcode_col:

  Character string specifying the column name in the Seurat metadata
  that contains the full cell barcode. Defaults to `"BC"`.

- condition_col:

  Character string specifying the name of the metadata column to store
  the condition labels in the Seurat object. Defaults to `"condition"`.

- Odd_barcode_col:

  Character string specifying the name of the metadata column to store
  the extracted Odd barcodes in the Seurat object. Defaults to
  `"Odd_barcode"`.

- Odd_barcode_position:

  Integer specifying the position of the Odd barcode after splitting the
  full barcode on underscores (`"_"`). Defaults to `2`.

## Value

The input `lseurat` object with updated Seurat metadata. Two new columns
are added to `lseurat$seurat@meta.data`: one containing the extracted
Odd barcodes and one containing the corresponding condition labels.

## Details

The full barcode is split using underscores, and the element at
`Odd_barcode_position` is interpreted as the Odd barcode. These Odd
barcodes are then matched to the metadata file using exact matching.
