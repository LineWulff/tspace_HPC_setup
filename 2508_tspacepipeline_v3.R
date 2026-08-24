### Rscript for running imputation and tspace trajectory ####
## Line Wulff, 25.08.27

.libPaths("/home/line.winthereik/software/miniforge3/envs/tspace_env_v2/lib/R/library")
# "/global/software/r/4.3.1/lib64/R/library",
lib_path <- "/home/line.winthereik/software/miniforge3/envs/tspace_env_v2/lib/R/library"
# for magicbatch
python_path <- system("which python3", intern = TRUE) 

##libraries
library(SeuratObject)
library(Signac)
library(Seurat)
library(magicBatch)
library(doSNOW)
library(doParallel)
source("/work/mccoy_lab/0000_generalscripts/tSp_helperfunctions.R")
source("/work/mccoy_lab/0000_generalscripts/tSpace.R")


dir <- getwd()
print(paste0("tSpace script started at: ",Sys.time()))
print(paste0("I am here: ",dir,""))
control_df <- read.csv(paste(dir,"/control_df.csv",sep=""), header = T)
control_df <- control_df[,c(1,2)]
colnames(control_df) <- c("parameter","value")
print(control_df)

#create output folder
if (file.exists("output")){print("Directory already exist.")} else {
    dir.create(file.path(dir, "output"))}

#variables from control csv file
obj_name <- as.character(control_df[control_df$parameter=="Seurat_obj",]$value)
T_impute <- control_df[control_df$parameter=="impute",]$value
imp_feat <- as.character(control_df[control_df$parameter=="magic_features",]$value)
assay <- as.character(control_df[control_df$parameter=="assay",]$value)
Pc_incl_max <- as.numeric(as.character(control_df[control_df$parameter=="pca_incl",]$value))
if (assay=="RNA"){Pc_incl_min <- 1} else if (assay %in% c("ATAC","ATAC/RNA") ){Pc_incl_min <- 2}
save_imputation <- control_df[control_df$parameter=="save_imputation",]$value
run_tspace <- control_df[control_df$parameter=="run_tspace",]$value
tsp_feat <- control_df[control_df$parameter=="tspace_features",]$value

#### read in Seurat data ####
obj <- readRDS(paste(dir, "/", obj_name, ".rds", sep=""))
#print(obj)

# find top expressed features if data is scATACseq
#if (assay=="ATAC"){
#  print("This is an ATAC seq assay, reclculating top exp. peaks")
#obj <- RunTFIDF(obj)
#obj <- FindTopFeatures(obj, min.cutoff = "q99", assay = "ATAC")
#length(VariableFeatures(obj)) # 10428, top %5 %
#obj <- RunSVD(obj, n = Pc_incl_max)
#}

##### magic imputation on normalized data ####
if (T_impute==TRUE){
#extract data from Seurat object
#decide which genes are to be imputated and set matrix based on input data
if (assay=="RNA"){
  if (imp_feat=="ALL"){
    var_feat <- rownames(obj[[obj@active.assay]])}
  else if (imp_feat=="variable"){
    var_feat <- VariableFeatures(obj)}
  magic_data <- FetchData(obj,vars = var_feat, assay="RNA")
  magic_data <- as.matrix(magic_data)}
else if (assay=="ATAC/RNA") {
    var_feat <- rownames(obj[["integrated"]])
    #magic_data <- FetchData(obj,vars = var_feat, assay="integrated")
    magic_data <- as.matrix(t(obj[["integrated"]]@data))}
else if (assay=="ATAC") {
  if (imp_feat=="variable"){
  obj <- FindVariableFeatures(object = obj, nfeatures = 1000)
  var_feat <- VariableFeatures(obj)}
  else if (imp_feat=="ALL"){
    var_feat <- rownames(obj[[obj@active.assay]])}
  magic_data <- FetchData(obj,vars = var_feat, assay="RNA")
  magic_data <- as.matrix(magic_data)}

if (assay=="RNA"){aff_magic_data <- obj@reductions$pca@cell.embeddings[,Pc_incl_min:Pc_incl_max]}
else if (assay=="ATAC"){aff_magic_data <- obj@reductions$lsi@cell.embeddings[,Pc_incl_min:Pc_incl_max]}
else if (assay=="ATAC/RNA"){aff_magic_data <- obj@reductions$lsi@cell.embeddings[,Pc_incl_min:Pc_incl_max]}

#magic imputation
magic_data <- magicBatch(data=magic_data, mar_mat_input=aff_magic_data, t_param=c(2,6), python_command = python_path)
#extract imputed matrix
imputed_t2 <- magic_data[[1]]$t2
imputed_t6 <- magic_data[[1]]$t6
print("imputation done, 5x5 matrix of t2:"); print(imputed_t2[1:5,1:5])

#add imputation data to Seurat obj
if (save_imputation==TRUE){
obj[["imputed_t2"]] <- CreateAssayObject(data = t(imputed_t2))
obj[["imputed_t6"]] <- CreateAssayObject(data = t(imputed_t6))
saveRDS(obj, file = paste(dir,"/output/",obj_name,"_imputed.rds",sep=""))}
} else {
print("Saving imputation was not requested.")} 

print(paste0("Imputation finished at: ",Sys.time()))

#### tSPACE ####
if (run_tspace==TRUE){
  print("You have asked to run tspce.")
#tspace_df <- FALSE
#which input to use, imputed variable features or adjusted pca
if (tsp_feat=="variable"){
  tspace_df <- imputed_t6[,VariableFeatures(obj)]
  dimred <- "both"
} else if (tsp_feat=="pca"){
  tspace_df <- obj@reductions$pca@cell.embeddings[,Pc_incl_min:Pc_incl_max]
  dimred <- "umap"
} else if (tsp_feat=="lsi"){
  print("You are running tspace on lsi coordinates.")
  tspace_df <- obj@reductions$lsi@cell.embeddings[,Pc_incl_min:Pc_incl_max]
  dimred <- "umap"
} else {print("You have not assigned a matrix for the tSpace calculation")}
  print("Here is a subset of your data input:")
  print(tspace_df[1:5,1:5])
class(tspace_df)
tspace_out <- tSpace(tspace_df, K = 20, graph = 5,
       trajectories = 20, wp = 20, ground_truth = F,
       weights = "exponential", dr = dimred, seed = 8, core_no = 20,
       env_path=lib_path)
#save files
print(paste0("tSpace finished running at: ",Sys.time()))
print("I made it here.")
print(paste(dir,"/output/",obj_name,"_tspacefile.rds",sep=""))
tspace_out
saveRDS(tspace_out,file = paste(dir,"/output/",obj_name,"_tspacefile.rds",sep=""))
tspace_out2 <- tspace_out$ts_file
saveRDS(tspace_out2,file = paste(dir,"/output/",obj_name,"_tspace_tsfile.rds",sep=""))
#} #tspace_df cannot be empty
} #end if run_tspace is true

print("All done!")
print(paste0("tSpace script finished at: ",Sys.time()))



