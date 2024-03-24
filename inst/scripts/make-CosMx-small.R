###########################################################################
#                CosMx demo dataset lung P9S1, fov 10, 11                 #
###########################################################################
cosmx_lung_p9s1_path <- "/Users/estelladong/Desktop/SampleData/Raw/Nanostring_CosMX/Lung"
counts <- read.csv(file.path(cosmx_lung_p9s1_path, "Lung5_Rep1_exprMat_file.csv")) # 91992   982
meta <- read.csv(file.path(cosmx_lung_p9s1_path, "Lung5_Rep1_metadata_file.csv"))  # 91972    20

library(dplyr)
counts_test <- counts %>%
  dplyr::filter(fov %in% c(10)) %>%
  slice(1:10) %>%
  dplyr::filter(cell_ID != 0) # take 9 cells

counts_test <- counts_test[, 1:10] # take 10 genes

meta_test <- meta %>%
  right_join(counts_test[, c("fov", "cell_ID")], by = c("fov", "cell_ID")) %>%
  select("fov", "cell_ID", "Area", "CenterX_local_px", "CenterY_local_px", "CenterX_global_px", "CenterY_global_px")

cosmx_lung_p9s1_demo_path <- "/Users/estelladong/Desktop/SpatialExperimentIO/inst/extdata/CosMx_small"
write.csv(counts_test, file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_exprMat_file.csv"), row.names = FALSE)
write.csv(meta_test, file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_metadata_file.csv"), row.names = FALSE)

#devtools::install_github("estellad/SpatialExperimentReader")
library(SpatialExperimentReader)
cos_spe <- readCosmxSXE(cosmx_lung_p9s1_demo_path, return_type = "SPE"); dim(cos_spe)
library(SpatialExperiment)
counts(cos_spe)[1:5, 1:5]

#remotes::install_github("estellad/ggspavis", ref = "estellad")
library(ggspavis)
cos_spe$sum <- colSums(counts(cos_spe))
plotSpots(cos_spe, annotate = "sum", in_tissue = NULL)


# # Whole_Image -------------------------------------------------------------
# cos_all_spe <- readCosmxSPE(cosmx_lung_p9s1_path)
# cos_all_spe$sum <- colSums(counts(cos_all_spe))
# plotSpots(cos_all_spe, annotate = "sum", in_tissue = NULL)
# 
# cos_sub <- cos_all_spe[, cos_all_spe$fov %in% c(10, 11)]
# cos_sub$libsize <- colSums(counts(cos_sub))
# plotSpots(cos_sub, annotate = "sum", in_tissue = NULL)
# plotSpots(cos_sub, annotate = "libsize", in_tissue = NULL)



