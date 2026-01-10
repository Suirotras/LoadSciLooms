# Create a custom "lseurat" object.

This function acts as a constructor for a custom "lseurat" object, which
is a named list containing a Seurat object, a list of parameters, and an
optional gmm object.

## Usage

``` r
create_lseurat_object(
  seurat_obj,
  params = list(),
  gmm = NA,
  loom_path = NA,
  id = NA
)
```

## Arguments

- seurat_obj:

  A valid SeuratObject. This is a required argument.

- params:

  A named list containing parameters. Defaults to an empty list.

- gmm:

  An optional gmm model object of class 'Mclust'. Defaults to NA.

- loom_path:

  A single-element character vector, representing the path to the loom
  file from which this object was generated.

- id:

  A single-element character vector. It is the unique ID chosen to
  represent this object (e.g. "i1_ES"). Defaults to NA.

## Value

A named list object of class "lseurat".

## Examples

``` r
# Example usage with dummy objects that have the correct class.
# We'll create simple lists and assign the correct classes to them for the example.

# Create a dummy SeuratObject
my_seurat <- list()
class(my_seurat) <- "Seurat"

# Create a dummy Mclust object
my_gmm <- list(bic = -434.9987, components = 2)
class(my_gmm) <- "Mclust"

# Define some parameters
all_params <- list(min_mapping_quality = 40, pdf_threshold = 0.0005)

# Scenario 1: Create a lseurat object with all arguments provided
lseurat_obj_full <- create_lseurat_object(
  seurat_obj = my_seurat,
  params = all_params,
  gmm = my_gmm,
  loom_path = "/path/to/my/file.loom",
  id = "i5-TDT3-d17"
)

# Scenario 2: Create a lseurat object without gmm, loom_path, or id
# (These arguments will use their default NA values)
lseurat_obj_basic <- create_lseurat_object(
  seurat_obj = my_seurat
)

# Check the class of the new objects
class(lseurat_obj_full)
#> [1] "lseurat"
class(lseurat_obj_basic)
#> [1] "lseurat"

# You can also check their structure
if (FALSE) { # \dontrun{
  str(lseurat_obj_full)
  str(lseurat_obj_basic)
} # }
```
