# Small everything

library(rhdf5)
## Xenium sample 1 
# h5 
h5file <- "/home/estelladong/Desktop/SampleData/archive/Xenium_Preprint_Data/Xenium/outs/cell_feature_matrix.h5"
h5ls(file = h5file)

h5_barcodes = h5read(h5file, "/matrix/barcodes")
h5_barcodes_new <- as(head(h5_barcodes), class(h5_barcodes)) #1 2 3 4 5 6 # take only 6 cells

h5_data = h5read(h5file, "/matrix/data")
# matrix(c(0, 2:6, 0, 8:10, 0, 0, 0, 14:15, 0, 17:18, 0, 20, 0, 2:3, 0), nrow = 4, ncol = 6) # something like this
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices
h5_data_new <- c(2, 3, 4, 5, 6, 8, 9, 10, 14, 15, 17, 18, 20, 2, 3) # int 1 1 1 3 2 4  # Nonzero UMI counts in column-major order


## Features
h5_all_tag_keys = h5read(h5file, "/matrix/features/_all_tag_keys")
h5_all_tag_keys_new <- h5_all_tag_keys # "genome"

h5_feature_type = h5read(h5file, "/matrix/features/feature_type")
h5_feature_type_new <- head(h5_feature_type, 4) # "Gene Expression" # take only 4 genes

h5_genome = h5read(h5file, "/matrix/features/genome")
h5_genome_new <- head(h5_genome, 4) # "Unknown" # take only 4 genes

h5_id = h5read(h5file, "/matrix/features/id")
h5_id_new <- head(h5_id, 4) # "ENSG000121270" # take only 4 genes

h5_name = h5read(h5file, "/matrix/features/name")
h5_name_new <- head(h5_name, 4) # "ABCC11" # take only 4 genes

# matrix
h5_indices = h5read(h5file, "/matrix/indices")
h5_indices_new <- h5_data_new <- c(1, 2, 3, 0, 1, 3, 0, 1, 1, 2, 0, 1, 3, 1, 2) # int 33 36 67  # Zero-based row index of corresponding element in `data`

h5_indptr = h5read(h5file, "/matrix/indptr")
h5_indptr_new <- c(0, 3, 6, 8, 10, 13, 15) # int 33 36 67 # take 6 cells, # Zero-based index into data / indices of the start of each column, i.e., the data corresponding to each barcode sequence
# https://stackoverflow.com/questions/52299420/scipy-csr-matrix-understand-indptr

h5_shape = h5read(h5file, "/matrix/shape")
h5_shape_new <- h5_shape # 541 167780
h5_shape_new[1] <- 4
h5_shape_new[2] <- 6

#               group          name       otype  dclass      dim
# 0                 /        matrix   H5I_GROUP                 
# 1           /matrix      barcodes H5I_DATASET  STRING   167780
# 2           /matrix          data H5I_DATASET INTEGER 10640759
# 3           /matrix      features   H5I_GROUP                 
# 4  /matrix/features _all_tag_keys H5I_DATASET  STRING        1
# 5  /matrix/features  feature_type H5I_DATASET  STRING      541
# 6  /matrix/features        genome H5I_DATASET  STRING      541
# 7  /matrix/features            id H5I_DATASET  STRING      541
# 8  /matrix/features          name H5I_DATASET  STRING      541
# 9           /matrix       indices H5I_DATASET INTEGER 10640759
# 10          /matrix        indptr H5I_DATASET INTEGER   167781
# 11          /matrix         shape H5I_DATASET INTEGER        2

smallh5 <- "/home/estelladong/Desktop/SampleData/SpatialExperimentIOData/Xenium/cell_feature_matrix.h5"

h5createFile(smallh5)
h5createGroup(smallh5,"matrix")
h5createGroup(smallh5,"matrix/features")
h5ls(smallh5)

h5write(h5_barcodes_new, smallh5,"matrix/barcodes")
h5write(h5_data_new, smallh5,"matrix/data")

h5write(h5_all_tag_keys_new, smallh5,"matrix/features/_all_tag_keys")
h5write(h5_feature_type_new, smallh5,"matrix/features/feature_type")
h5write(h5_genome_new, smallh5,"matrix/features/genome")
h5write(h5_id_new, smallh5,"matrix/features/id")
h5write(h5_name_new, smallh5,"matrix/features/name")

h5write(h5_indices_new, smallh5,"matrix/indices")
h5write(h5_indptr_new, smallh5,"matrix/indptr")
h5write(h5_shape_new, smallh5,"matrix/shape")
h5ls(smallh5)


# Now try to read it as SCE
library(SingleCellExperiment)
# chromh5 <- "/home/estelladong/Desktop/SampleData/archive/Xenium_Preprint_Data/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer/count_sample_filtered_feature_bc_matrix.h5"

sce <- DropletUtils::read10xCounts(h5file, type = "HDF5", col.names = TRUE)
# sce <- DropletUtils::read10xCounts(chromh5, type = "HDF5", col.names = TRUE)
scetest <- SingleCellExperiment::SingleCellExperiment(assays = assays(sce))
counts(scetest)[34:38,]

sce <- DropletUtils::read10xCounts(smallh5, type = "HDF5", col.names = TRUE)
scetest <- SingleCellExperiment::SingleCellExperiment(assays = assays(sce))

data_check <- counts(sce)

read.csv()

############
library(dplyr)
# Read un .gz ed cell.csv
cells <- read.csv("~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/cells.csv")
cells <- cells %>% slice(1:6)

write.csv(cells, "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/cells.csv")

################
library(SpatialExperimentIO)
spe <- readXeniumSXE("~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/", return_type = "SPE")
plotSpots(spe, annotate = "total_counts", in_tissue = NULL)






