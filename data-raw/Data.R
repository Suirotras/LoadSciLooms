### This script generates the subsampled exmpla loom files

subset_loom_file <- function(full_loom_path, subsampled_loom_path,
                             cell_subset, gene_subset) {
  # Create hdf5 file
  rhdf5::H5Fcreate(subsampled_loom_path, flags = "H5F_ACC_EXCL")

  full_loom <- rhdf5::H5Fopen(full_loom_path,
                              flags = "H5F_ACC_RDONLY")

  # Get the group names
  full_loom_content <- rhdf5::h5ls(full_loom)
  full_loom_groups <- full_loom_content$name[full_loom_content$otype == "H5I_GROUP"]

  # Create groups
  for (group in full_loom_groups) {
    rhdf5::h5createGroup(subsampled_loom_path, group)
  }

  # Write the single cell matrix
  rhdf5::h5write(obj = full_loom$matrix[cell_subset, gene_subset],
                 file = subsampled_loom_path,
                 name = "matrix")

  # Add the metadata to the new hdf5 file
  for (group in full_loom_groups) {

    expr_string <- paste0("names(full_loom$", group, ")")
    group_attrs_names <- eval(parse(text = expr_string))

    # If not metadata found in group, continue to next iteration
    if (length(group_attrs_names) == 0) {
      message(paste0("skipping group '", group, "', as it is empty"))
      next
    }

    for (grp_attr in group_attrs_names) {

      expr_str_grp <- paste0("full_loom$\"", group, "/", grp_attr, "\"")
      group_attr <- eval(parse(text = expr_str_grp))

      # Subset the attribute the same way as the matrix
      if (group == "col_attrs") {
        group_attr <- group_attr[cell_subset]
      } else if (group == "row_attrs") {
        group_attr <- group_attr[gene_subset]
      }

      # Write the metadata from the attribute to the hdf5 file
      rhdf5::h5write(obj = group_attr,
                     file = subsampled_loom_path,
                     name = paste0(group, "/", grp_attr))
    }
  }
  # Close all remaining hdf5 connections
  rhdf5::h5closeAll()
}

# Create raw data directory
dir.create("inst/extdata/", recursive = TRUE)

subset_loom_file(full_loom_path = paste0("/Users/jari/Hubrecht_data/projects/",
                                         "NPC/KIN9199_TDT3/results/",
                                         "DamID_snakemake/",
                                         "libraries/i1-TDT3-ES/",
                                         "transcriptome_se.loom"),
                 subsampled_loom_path = "inst/extdata/i1_ES_subsample.loom",
                 cell_subset = 1:150, gene_subset = 1:5000)

subset_loom_file(full_loom_path = paste0("/Users/jari/Hubrecht_data/projects/",
                                         "NPC/KIN9199_TDT3/results/",
                                         "DamID_snakemake/",
                                         "libraries/i3-TDT3-d6/",
                                         "transcriptome_se.loom"),
                 subsampled_loom_path = "inst/extdata/i3_d6_subsample.loom",
                 cell_subset = 1:150, gene_subset = 1:5000)

subset_loom_file(full_loom_path = paste0("/Users/jari/Hubrecht_data/projects/",
                                         "NPC/KIN9199_TDT3/results/",
                                         "DamID_snakemake/",
                                         "libraries/i5-TDT3-d17/",
                                         "transcriptome_se.loom"),
                 subsampled_loom_path = "inst/extdata/i5_d17_subsample.loom",
                 cell_subset = 1:150, gene_subset = 1:5000)

subset_loom_file(full_loom_path = paste0("/Users/jari/Hubrecht_data/projects/",
                                         "SciT_snakemake/KIN10808/data/",
                                         "sciT_snakemake_results/",
                                         "libraries/i31-FV-GSK126-1K/",
                                         "transcriptome_se.loom"),
                 subsampled_loom_path = "inst/extdata/i31_GSK126_subsample.loom",
                 cell_subset = 1:150, gene_subset = 1:5000)
