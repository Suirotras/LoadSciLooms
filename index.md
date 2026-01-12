# LoadSciLooms

*A companion R package for the sciT-snakemake workflow*

## Overview

**LoadSciLooms** is an R package designed as a **companion to the
[`SciT_snakemake`](https://github.com/Suirotras/SciT_snakemake)
pipeline**.

While **SciT_snakemake** performs the heavy lifting of preprocessing
single-cell combinatorial indexing transcriptomics (**sciT**) data
(demultiplexing, alignment, counting, QC, and loom generation),
**LoadSciLooms** provides the **R-side tooling** needed to:

- Load sciT-generated **`.loom` files** into R
- Convert them into **Seurat objects**
- Optionally perform **Gaussian Mixture Model (GMM)**–based cell calling
- Merge multiple sciT libraries while preserving provenance and metadata
- Keep analysis parameters and model objects bundled with the data

In short:

> **SciT_snakemake → `.loom` files → LoadSciLooms → Seurat-based
> analysis**

------------------------------------------------------------------------

## Relationship to `SciT_snakemake`

| SciT_snakemake          | LoadSciLooms                |
|-------------------------|-----------------------------|
| Snakemake pipeline      | R package                   |
| FASTQ → BAM → counts    | `.loom` → Seurat            |
| Alignment, QC, counting | Loading, filtering, merging |
| Produces `.loom` files  | Consumes `.loom` files      |

Thus, this package can be used to load sciT-snakemake-processed data
into R for further analysis.

## Typical workflow

    FASTQ files
    ↓
    SciT_snakemake
    ↓
    Transcriptome .loom files
    ↓
    LoadSciLooms (this package)
    ↓
    Seurat objects
    ↓
    Downstream analysis (clustering, annotation, visualization)

## Installation

``` r
library(devtools)
devtools::install_github("Suirotras/LoadSciLooms", ref = "master")
```

## Required R dependencies

- SeuratObject
- rhdf5
- Matrix
- mclust

These are automatically installed when missing.

## The lseurat object

All major functions in this package return an object of class lseurat.

An lseurat object is a lightweight container that keeps data,
parameters, and models together.

### Structure

| Slot      | Description                                    |
|-----------|------------------------------------------------|
| seurat    | A Seurat object containing counts and metadata |
| params    | Parameters used to create the object           |
| gmm       | Mclust model used for cell calling             |
| loom_path | Path to the source .loom file                  |
| id        | Sample / library identifier                    |

## Loading a sciT loom file into R

After running SciT_snakemake, load a transcriptome loom file:

``` r
# Load as seuratobject
Loom_path <- "path/to/libraries/i31-sample1-1K/transcriptome_se.loom"
lseurat <- LoomAsSeurat(Loom_path, matrix_rowname_col = "Gene",
                        resolve_duplicates = TRUE,
                        gmm_cell_calling = FALSE)
```

This creates:

- A sparse transcriptome matrix
- A Seurat object with cell and feature metadata
- An `lseurat` container capturing multiple parameters

## Optional: GMM-based cell calling

SciT data often contains a mixture of real cells and background
barcodes. LoadSciLooms optionally performs Gaussian Mixture Model–based
cell calling using total transcript counts.

``` r
lseurat <- LoomAsSeurat(
  loom_path = "transcriptome.loom",
  gmm_cell_calling = TRUE,
  gmm_k = 1:6,
  gmm_verbose = TRUE
)
```

Results:

- Cell classifications stored in `lseurat$seurat$gmm`
- Fitted model stored in `lseurat$gmm`

## Loading and merging multiple sciT libraries

For experiments processed as multiple sciT libraries:

``` r
loom_paths <- c(
  "i31-sample1-1K" = "path/to/i31-sample1-1K.loom",
  "i31-sample2-10K" = "path/to/i31-sample2-10K.loom"
)

lseurat_list <- MultiLoomAsSeurat(
  loom_paths = loom_paths,
  merge_seurat = TRUE,
  remove_unmerged = TRUE
)
```

Output

- One `lseurat` object per input loom
- A merged Seurat object at `lseurat_list$merged_seurat`
- A `loom` metadata column indicating sample of origin

## Duplicate gene handling

SciT loom files from the sciT-snakemake workflow have both gene symbols
and ENSEMBL gene IDs.

When using the gene symbols to create the seurat object, duplicate gene
symbols become problematic. Since Seurat does not allow duplicate
feature names:

- Default behavior: error
- Optional behavior: automatically resolve duplicates

``` r
# Load as seuratobject
Loom_path <- "path/to/libraries/i31-sample1-1K/transcriptome_se.loom"
lseurat <- LoomAsSeurat(Loom_path, matrix_rowname_col = "name",
                        resolve_duplicates = TRUE,
                        gmm_cell_calling = FALSE)
```

Duplicate features are resolved by keeping the first occurrence.

## When should I use this package?

Use **LoadSciLooms** if:

- You ran **SciT_snakemake**
- You have transcriptome `.loom` files.
- You want to analyze data in Seurat.

## Related repository

- SciT_snakemake (<https://github.com/Suirotras/SciT_snakemake>)
