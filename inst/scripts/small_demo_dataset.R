###########################################################################
#                CosMx demo dataset lung P9S1, fov 10, 11                 #
###########################################################################

cosmx_lung_p9s1_path <- "~/Desktop/SampleData/Raw/Nanostring_CosMx/Lung/P9/S1"
counts <- read.csv(file.path(cosmx_lung_p9s1_path, "cell_by_gene.csv")) # 91992   982
meta <- read.csv(file.path(cosmx_lung_p9s1_path, "cell_metadata.csv"))  # 91972    20

library(dplyr)
counts_test <- counts %>%
  dplyr::filter(fov %in% c(10, 11)) %>%
  dplyr::filter(cell_ID != 0)# take 2 fov, around 1/10 of the original # of cells

meta_test <- meta %>%
  right_join(counts_test[, c("fov", "cell_ID")], by = c("fov", "cell_ID")) %>%
  select("fov", "cell_ID", "Area", "CenterX_local_px", "CenterY_local_px", "CenterX_global_px", "CenterY_global_px")

cosmx_lung_p9s1_demo_path <- "~/Desktop/SampleData/Raw/Nanostring_CosMx/Lung/P9/test_small_demo"
write.csv(counts_test, file.path(cosmx_lung_p9s1_demo_path, "small_lung_p9s1_exprMat_file.csv"), row.names = FALSE)
write.csv(meta_test, file.path(cosmx_lung_p9s1_demo_path, "small_lung_p9s1_metadata_file.csv"), row.names = FALSE)

#devtools::install_github("estellad/SpatialExperimentReader")
library(SpatialExperimentReader)
cos_spe <- readCosmxSPE(cosmx_lung_p9s1_demo_path); dim(cos_spe)
library(SpatialExperiment)
counts(cos_spe)[1:5, 1:5]

#remotes::install_github("estellad/ggspavis", ref = "estellad")
library(ggspavis)
cos_spe$sum <- colSums(counts(cos_spe))
plotSpots(cos_spe, annotate = "sum", in_tissue = NULL)


# Whole_Image -------------------------------------------------------------
cos_all_spe <- readCosmxSPE(cosmx_lung_p9s1_path)
cos_all_spe$sum <- colSums(counts(cos_all_spe))
plotSpots(cos_all_spe, annotate = "sum", in_tissue = NULL)

cos_sub <- cos_all_spe[, cos_all_spe$fov %in% c(10, 11)]
cos_sub$libsize <- colSums(counts(cos_sub))
plotSpots(cos_sub, annotate = "sum", in_tissue = NULL)
plotSpots(cos_sub, annotate = "libsize", in_tissue = NULL)



###########################################################################
#     MERSCOPE demo dataset ovarian P2S2, coord x > 10000, y > 8000       #
###########################################################################
mer_ovarian_p2s2_path <- "~/Desktop/SampleData/Raw/Vizgen_MERSCOPE/Ovarian/P2/S2"

# Whole_Image -------------------------------------------------------------
mer_all_spe <- readMerscopeSPE(mer_ovarian_p2s2_path)
mer_all_spe$sum <- colSums(counts(mer_all_spe))
plotSpots(mer_all_spe, annotate = "sum", in_tissue = NULL)

mer_sub <- mer_all_spe[, spatialCoords(mer_all_spe)[, 1] > 10000 & spatialCoords(mer_all_spe)[, 2] > 8000]; dim(mer_sub)
mer_sub$libsize <- colSums(counts(mer_sub))
plotSpots(mer_sub, annotate = "sum", in_tissue = NULL)
plotSpots(mer_sub, annotate = "libsize", in_tissue = NULL)


head(spatialCoords(mer_all_spe))
table(spatialCoords(mer_all_spe)[, 1] > 10000 & spatialCoords(mer_all_spe)[, 2] > 8000) # this will be the criteria


library(viridis)
df <- cbind.data.frame(data.frame(spatialCoords(mer_all_spe)),
                       var = mer_all_spe$sum)
p <- ggplot(df, aes(x = center_x, y = center_y, color = var)) +
  geom_point(size = 0.1) + 
  theme_bw() + 
  scale_color_viridis(option = "magma", name = "cor") + 
  theme(legend.title.align = 0.5, plot.title = element_text(hjust = 0.5)) 

source("~/Desktop/plotSpots_ED.R")
plotSpots_ED(mer_sub, annotate = "libsize", in_tissue = NULL, show_axis = TRUE)



# Small image -------------------------------------------------------------
mer_ovarian_p2s2_path <- "~/Desktop/SampleData/Raw/Vizgen_MERSCOPE/Ovarian/P2/S2"
counts <- read.csv(file.path(mer_ovarian_p2s2_path, "OvarianP2S2_cell_by_gene.csv")) # 91992   982
meta <- read.csv(file.path(mer_ovarian_p2s2_path, "OvarianP2S2_cell_metadata.csv"))  # 91972    20

library(dplyr)
meta_test <- meta %>%
  dplyr::filter(center_x > 10000 & center_y > 8000) %>%
  select("X", "fov", "volume", "center_x", "center_y")

counts_test <- counts %>%
  mutate(X = cell) %>%
  right_join(meta_test[, c("X", "fov")], by = "X") %>%
  select(-c(X, fov))

mer_ovarian_p2s2_demo_path <- "~/Desktop/SampleData/Raw/Vizgen_MERSCOPE/Ovarian/P2/p2s2_demo"
write.csv(counts_test, file.path(mer_ovarian_p2s2_demo_path, "small_ovarian_p2s2_cell_by_gene.csv"), row.names = FALSE)
write.csv(meta_test, file.path(mer_ovarian_p2s2_demo_path, "small_ovarian_p2s2_cell_metadata.csv"), row.names = FALSE)



















