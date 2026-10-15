# Run from mixed_paternity_final; export the exact inputs of the final model.
full <- readRDS("models/common_scaling/model_full_brouwer_griffith.rds")
matched_data <- full$data
matched_phylo_mat <- full$data2$phylo_mat
stopifnot(nrow(matched_data) == 533L,
          length(unique(matched_data$Study_ID)) == 481L,
          length(unique(matched_data$species_nonphylo)) == 326L,
          !is.null(matched_phylo_mat),
          all(as.character(matched_data$species_phylo_brms) %in%
                rownames(matched_phylo_mat)))
dir.create("data", showWarnings = FALSE)
saveRDS(matched_data, "data/data_null_matched_533.rds")
saveRDS(matched_phylo_mat, "data/phylo_mat_matched_533.rds")
saved_data <- readRDS("data/data_null_matched_533.rds")
# The terms environment is re-created by deserialisation; verify all values.
stopifnot(isTRUE(all.equal(saved_data, full$data)),
          all(vapply(names(saved_data), function(n)
            identical(saved_data[[n]], full$data[[n]]), logical(1))),
          identical(readRDS("data/phylo_mat_matched_533.rds"), full$data2$phylo_mat))
cat("Verified: 533 observations; 481 study identities; 326 species.\n")
print(tools::md5sum(c("data/data_null_matched_533.rds",
                     "data/phylo_mat_matched_533.rds")))
