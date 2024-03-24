###########################################################################
#     MERSCOPE demo dataset ovarian P2S2, coord x > 10000, y > 8000       #
###########################################################################
mer_ovarian_p2s2_path <- "/Users/estelladong/Desktop/SampleData/Raw/Vizgen_MERSCOPE/Ovarian/P2/S1"

# # Whole_Image -------------------------------------------------------------
# mer_all_spe <- readMerscopeSPE(mer_ovarian_p2s2_path)
# mer_all_spe$sum <- colSums(counts(mer_all_spe))
# plotSpots(mer_all_spe, annotate = "sum", in_tissue = NULL)
# 
# mer_sub <- mer_all_spe[, spatialCoords(mer_all_spe)[, 1] > 10000 & spatialCoords(mer_all_spe)[, 2] > 8000]; dim(mer_sub)
# mer_sub$libsize <- colSums(counts(mer_sub))
# plotSpots(mer_sub, annotate = "sum", in_tissue = NULL)
# plotSpots(mer_sub, annotate = "libsize", in_tissue = NULL)
# 
# 
# head(spatialCoords(mer_all_spe))
# table(spatialCoords(mer_all_spe)[, 1] > 10000 & spatialCoords(mer_all_spe)[, 2] > 8000) # this will be the criteria
# 
# 
# library(viridis)
# df <- cbind.data.frame(data.frame(spatialCoords(mer_all_spe)),
#                        var = mer_all_spe$sum)
# p <- ggplot(df, aes(x = center_x, y = center_y, color = var)) +
#   geom_point(size = 0.1) + 
#   theme_bw() + 
#   scale_color_viridis(option = "magma", name = "cor") + 
#   theme(legend.title.align = 0.5, plot.title = element_text(hjust = 0.5)) 
# 
# source("~/Desktop/plotSpots_ED.R")
# plotSpots_ED(mer_sub, annotate = "libsize", in_tissue = NULL, show_axis = TRUE)



# Small image -------------------------------------------------------------
mer_ovarian_p2s2_path <- "/Users/estelladong/Desktop/SampleData/Raw/Vizgen_MERSCOPE/Ovarian/P2/S1"
counts <- read.csv(file.path(mer_ovarian_p2s2_path, "OvarianP2S1_cell_by_gene.csv")) # 91992   982
meta <- read.csv(file.path(mer_ovarian_p2s2_path, "OvarianP2S1_cell_metadata.csv"))  # 91972    20

library(dplyr)
meta_test <- meta %>%
  dplyr::filter(center_x > 12500 & center_y > 10000) %>%
  select("X", "fov", "volume", "center_x", "center_y")

counts_test <- counts %>%
  mutate(X = cell) %>%
  right_join(meta_test[, c("X", "fov")], by = "X") %>%
  select(-c(X, fov)) 

counts_test <- counts_test[, 1:10]

mer_ovarian_p2s2_demo_path <- "/Users/estelladong/Desktop/SpatialExperimentIO/inst/extdata/MERSCOPE_small"
write.csv(counts_test, file.path(mer_ovarian_p2s2_demo_path, "ovarian_p2s2_cell_by_gene.csv"), row.names = FALSE)
write.csv(meta_test, file.path(mer_ovarian_p2s2_demo_path, "ovarian_p2s2_cell_metadata.csv"), row.names = FALSE)

# ReadSPE
spe <- readMerscopeSXE(mer_ovarian_p2s2_demo_path, return_type = "SPE")
plotSpots(spe, annotate = "volume", in_tissue = NULL)

