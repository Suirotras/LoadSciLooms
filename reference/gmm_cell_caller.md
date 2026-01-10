# Function for calling cells using gaussian mixture models

Function for calling cells using gaussian mixture models

## Usage

``` r
gmm_cell_caller(counts, k = 1:6, log_convert = TRUE, verbose = TRUE)
```

## Arguments

- counts:

  When 'counts' is an atomic integer vector, the values should represent
  total read counts per cell (e.g. cDNA or gDNA counts). A matrix or
  dataframe can also be given, in which case the rows correspond with
  cells and the columns represent separate modalities (e.g. RNA and
  DamID). Thus, cell class membership can be determined based on
  multiple modalities.

- k:

  An atomic vector of type integer. A gmm is estimated for each value of
  'k', where this value represents the number of components in the gmm.
  The cell classifications will be based on the model with the lowest
  BIC. To prevent testing of multiple components and to select this
  yourself, provide a single integer. Default 'k' is '1:6'.

- log_convert:

  Logical value indicating if the counts should be log-transformed
  before model estimation.

- verbose:

  Logical value indicating if the function should message the chosen
  number of components and the corresponding bic.

## Value

A two-element list. The first element is an atomic integer vector, where
the integers represent class membership to one of the predicted classes
(e.g. 1: 'no cell', 2: 'real cell'). The second element represents the
'Mclust' object that was used for the class membership prediction.
