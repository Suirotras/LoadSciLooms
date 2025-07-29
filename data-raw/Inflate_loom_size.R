### This script inflates the size of a loom file by a certain multiplication
### factor

inflate_loom_file <- function(full_loom_path, inflated_loom_path,
                             mlt_factor = 35) {
  # Create hdf5 file
  rhdf5::H5Fcreate(inflated_loom_path, flags = "H5F_ACC_EXCL")

  full_loom <- rhdf5::H5Fopen(full_loom_path,
                              flags = "H5F_ACC_RDONLY")

  # Get the group names
  full_loom_content <- rhdf5::h5ls(full_loom)
  full_loom_groups <- full_loom_content$name[full_loom_content$otype == "H5I_GROUP"]

  # Create groups
  for (group in full_loom_groups) {
    rhdf5::h5createGroup(inflated_loom_path, group)
  }

  # Retrieve matrix from file
  lmatrix <- full_loom$matrix

  # Create expression to generate the inflated matrix
  inf_mtx_expr <- parse(text = paste0("rbind(", paste(rep("lmatrix", mlt_factor), collapse = ", "), ")"))

  # Write the inflated single cell matrix
  rhdf5::h5write(obj = eval(inf_mtx_expr),
                 file = inflated_loom_path,
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

      if (group == "col_attrs") {
        # If attribute is 'CellID', make sure there are no duplicate names
        if (grp_attr == "CellID") {
          # Inflate the column attributes
          group_attr <- rep(group_attr, times = mlt_factor)
          group_attr <- make.unique(group_attr)
        } else {
          # Inflate the column attributes
          group_attr <- rep(group_attr, times = mlt_factor)
        }
      }

      # Write the metadata from the attribute to the hdf5 file
      rhdf5::h5write(obj = group_attr,
                     file = inflated_loom_path,
                     name = paste0(group, "/", grp_attr))
    }
  }
  # Close all remaining hdf5 connections
  rhdf5::h5closeAll()
}

# Create raw data directory
dir.create("/Users/jari/Hubrecht_data/projects/Miscellaneous/LoadSciLooms/libraries/i1-TDT3-ES/",
           recursive = TRUE)

inflate_loom_file(full_loom_path = paste0("/Users/jari/Hubrecht_data/projects/",
                                         "NPC/KIN9199_TDT3/results/",
                                         "DamID_snakemake/",
                                         "libraries/i1-TDT3-ES/",
                                         "transcriptome_se.loom"),
                 inflated_loom_path = paste0("/Users/jari/Hubrecht_data/projects/",
                                             "Miscellaneous/LoadSciLooms/libraries/",
                                             "i1-TDT3-ES/transcriptome_se_35x_inflated.loom"),
                 mlt_factor = 35)
