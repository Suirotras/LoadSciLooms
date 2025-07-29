### Test variations of the transpose function
### Conclusion: No indexing parantheses
###             First conversion to matrix, than transpose.
###             This combination is the fastest and takes the least memory
# Transpose after sparse conversion with parentheses
lpath <- "/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/transcriptome_se_35x_inflated.loom"
lfile <- rhdf5::H5Fopen(lpath, flags = "H5F_ACC_RDONLY")
lmatrix_sparse <- Matrix::t(Matrix::Matrix(lfile$matrix[,], sparse = TRUE))

# Transpose before sparse conversion with parentheses
lpath <- "/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/transcriptome_se_35x_inflated.loom"
lfile <- rhdf5::H5Fopen(lpath, flags = "H5F_ACC_RDONLY")
lmatrix_sparse <- Matrix::Matrix(t(lfile$matrix[,]), sparse = TRUE)

# Transpose after sparse conversion with parentheses
lpath <- "/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/transcriptome_se_35x_inflated.loom"
lfile <- rhdf5::H5Fopen(lpath, flags = "H5F_ACC_RDONLY")
lmatrix_sparse <- Matrix::t(Matrix::Matrix(lfile$matrix, sparse = TRUE))

# Transpose before sparse conversion with parentheses
lpath <- "/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/transcriptome_se_35x_inflated.loom"
lfile <- rhdf5::H5Fopen(lpath, flags = "H5F_ACC_RDONLY")
lmatrix_sparse <- Matrix::Matrix(t(lfile$matrix), sparse = TRUE)

### Test Loom loading function
### Conclusion: Made sure that dimensions of matrix is only retrieved once from
### the loom file during checks, instead of twice.
devtools::load_all()
lpath <- "/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/transcriptome_se_35x_inflated.loom"
lseurat <- LoomAsSeurat(lpath)

### Test Multi-Loom loading function
### Conclusion: Not really anything to improve. It mostly just uses
### LoomAsSeurat, so improvements thare apply here as well.
devtools::load_all()
lpaths <- c("i1_ES_infl" = "/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/transcriptome_se_35x_inflated.loom",
            "i3_d6" = system.file("extdata", "i3_d6_subsample.loom",
                                  package = "LoadSciLooms"),
            "i5_d17" = system.file("extdata", "i5_d17_subsample.loom",
                                   package = "LoadSciLooms"))

lseurat <- MultiLoomAsSeurat(lpaths)



