#' Obtain the MERSCOPE data
#'
#' Obtain a small demo of MERSCOPE ovarian patient 2 slice 2 data from: 
#' https://info.vizgen.com/ffpe-showcase
#'
#' @details
#' * cell_by_gene.csv 
#' https://console.cloud.google.com/storage/browser/_details/vz-ffpe-showcase/HumanOvarianCancerPatient2Slice2/cell_by_gene.csv?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
#' * cell_metadata.csv
#' https://console.cloud.google.com/storage/browser/_details/vz-ffpe-showcase/HumanOvarianCancerPatient2Slice2/cell_metadata.csv?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
#'
#' @author Estella Yixing Dong
#' 
#' @references
#' TODO: how to deal with not publication?
#'
#' @return A \linkS4class{SpatialExperiment} object with a matrix of UMI counts 
#' and spatial coordinates
#' 
#' @export
#' 
#' @importFrom ...
mockMerscope <- function(){

}