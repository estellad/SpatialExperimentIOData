# Small image -------------------------------------------------------------
starmap_mousebrain_path <- "/Users/estelladong/Desktop/SampleData/Raw/STARmap_PLUS"
counts <- read.csv(file.path(starmap_mousebrain_path, "well05raw_expression_pd.csv")) # 91992   982
meta <- read.csv(file.path(starmap_mousebrain_path, "well05_spatial.csv"))  # 91972    20

library(dplyr)
meta_test <- meta %>%
  slice(1:10)  # take 9 cells 

cells <- meta_test$NAME[meta_test$NAME != "TYPE"]

counts_test <- counts %>%
  slice(1:8) # take 8 genes

counts_test <- counts_test[, colnames(counts_test) %in% c("GENE", cells)]

starmap_mousebrain_demo_path <- "/Users/estelladong/Desktop/SpatialExperimentIO/inst/extdata/STARmapPLUS_small"
write.csv(counts_test, file.path(starmap_mousebrain_demo_path, "mockraw_expression_pd.csv"), row.names = FALSE)
write.csv(meta_test, file.path(starmap_mousebrain_demo_path, "mock_spatial.csv"), row.names = FALSE)

# ReadSPE
spe <- readStarmapplusSXE(starmap_mousebrain_demo_path, return_type = "SPE")
plotSpots(spe, annotate = "Main_molecular_cell_type", in_tissue = NULL)
