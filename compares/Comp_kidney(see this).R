# Kidney comparisons

# Codes of Garnett, CellAssign, and SingleR about kidney dataset.
# Codes are similar for all datasets. Thus we only shows codes for one dataset.
# For packages, see 'functions/packages.R'.

### 1 input data
ann <- read.table("annotations_droplet.csv",sep=",",header=TRUE)
ann[,2] <- as.factor(ann[,2])
ann_cluster <- ann[,2] %>% as.numeric
ann[,4] <- ann_cluster

GM_kidney <- read.csv("cellmatch_kidney.csv")
barcodes_kidney <- read.table("droplet\\Kidney-10X_P4_5\\barcodes.tsv", sep = "-")[,1]
genes_kidney <- read.table("droplet\\Kidney-10X_P4_5\\genes.tsv")[,1]
matrix_kidney <- Matrix::readMM("droplet\\Kidney-10X_P4_5\\matrix.mtx")
barcodes_kidney <- paste("10X_P4_5", barcodes_kidney, sep="_")


barcodes_kidney2 <- read.table("droplet\\Kidney-10X_P4_6\\barcodes.tsv", sep = "-")[,1]
matrix_kidney2 <- Matrix::readMM("droplet\\Kidney-10X_P4_6\\matrix.mtx")
barcodes_kidney2 <- paste("10X_P4_6", barcodes_kidney2, sep="_")

barcodes_kidney3 <- read.table("droplet\\Kidney-10X_P7_5\\barcodes.tsv", sep = "-")[,1]
matrix_kidney3 <- Matrix::readMM("droplet\\Kidney-10X_P7_5\\matrix.mtx")
barcodes_kidney3 <- paste("10X_P7_5", barcodes_kidney3, sep="_")


B_kidney <- c(barcodes_kidney,barcodes_kidney2,barcodes_kidney3)
genes_kidney <- toupper(genes_kidney)


Matrix_kidney <- cbind(matrix_kidney, matrix_kidney2, matrix_kidney3)
Matrix_kidney <- Seurat::LogNormalize(data = Matrix_kidney, scale.factor = 10000, verbose = F)
na_index <- which(apply(Matrix_kidney,1,sum) == 0)
Matrix_kidney <- Matrix_kidney[-na_index,]
genes_kidney <- genes_kidney[-na_index]
DFtest_kidney <- Matrix_kidney
rownames(DFtest_kidney) <- genes_kidney

Y_kidney <- ann[match(B_kidney, ann[,1]),][,4]
na.Y_kidney <- which(is.na(Y_kidney))
Y_kidney <- Y_kidney[-na.Y_kidney]
DFtest_kidney <- DFtest_kidney[,-na.Y_kidney]
colnames(DFtest_kidney) <- B_kidney[-na.Y_kidney]
class_kidney <- read.csv("ann_kidney.csv")
Y_kidney.raw <- class_kidney[match(Y_kidney, class_kidney[,1]),2]
DFtest_kidney <- as.matrix(DFtest_kidney)
DFtest_kidney <- as(DFtest_kidney,'dgCMatrix')
rm(Matrix_kidney,matrix_kidney,matrix_kidney2,matrix_kidney3)
rm(barcodes_kidney,barcodes_kidney2,barcodes_kidney3)


### 2 marker files

markers_kidney <- 'markers_kidney.txt'
bigGM_kidney <- read.csv('bigGM_kidney.csv', header = T,check.names = FALSE)
rownames(bigGM_kidney) <- bigGM_kidney[,1]
bigGM_kidney <- bigGM_kidney[,-1]
colnames(bigGM_kidney) <- unique(GM_kidney[,2])
bigGM_kidney <- bigGM_kidney[which(rownames(bigGM_kidney) %in% genes_kidney),]

### Garnett
pdata_kidney <- data.frame(FACS_type = Y_kidney.raw)
rownames(pdata_kidney) <- colnames(DFtest_kidney)
fdata_kidney <- data.frame(gene_short_name = rownames(DFtest_kidney))
rownames(DFtest_kidney) <- str_to_title(rownames(DFtest_kidney))
rownames(fdata_kidney) <- rownames(DFtest_kidney)
cds_kidney <- new_cell_data_set(DFtest_kidney,
                                cell_metadata = pdata_kidney,
                                gene_metadata = fdata_kidney)
rownames(DFtest_kidney) <- toupper(rownames(DFtest_kidney))
set.seed(1)
kidney_classifier <- train_cell_classifier(cds = cds_kidney,
                                           marker_file = markers_kidney,
                                           db = org.Mm.eg.db,
                                           cds_gene_id_type = "SYMBOL",
                                           num_unknown = 50,
                                           marker_file_gene_id_type = "SYMBOL")
set.seed(1)
cds1_kidney <- classify_cells(cds_kidney, kidney_classifier,
                              db = org.Mm.eg.db,
                              cluster_extend = TRUE,
                              cds_gene_id_type = "SYMBOL")

length(which(pData(cds1_kidney)$cluster_ext_type == pData(cds1_kidney)$FACS_type))



### cellassign
sce_kidney <- SingleCellExperiment(assays = list(counts = DFtest_kidney),
                                   colData = Y_kidney.raw,
                                   rowData = rownames(DFtest_kidney))
gene_index_kidney <- which(rownames(DFtest_kidney) %in% rownames(bigGM_kidney))
# The reason to use edgeR::calcNormFactors, see \url{https://github.com/Irrationone/cellassign/issues/23}.
sce_kidney_Size_factors <- edgeR::calcNormFactors(counts(sce_kidney))

sce_kidney <- sce_kidney[gene_index_kidney,]



library(reticulate)
# GPU options changed to accelerate cellassign
repl_python()
import tensorflow as tf
config = tf.compat.v1.ConfigProto(allow_soft_placement=True)
config.gpu_options.allow_growth = True
config.gpu_options.per_process_gpu_memory_fraction = 0.8
sess = tf.compat.v1.Session(config=config)
tf.compat.v1.keras.backend.set_session(sess)
exit

Y_kidney_cellassign <- rep(0, length(Y_kidney.raw))
set.seed(1)
cuttest_kidney <- cut(sample(1:length(Y_kidney.raw)),breaks = 5,labels = F)

# If CellAssign is done jointly, use this.
# sce_kidney_fit <- cellassign(exprs_obj = sce_kidney, 
#                              marker_gene_info = as.matrix(bigGM_kidney[which(rownames(bigGM_kidney) %in% rownames(sce_kidney)),]),
#                              s = sce_kidney_Size_factors, 
#                              shrinkage = TRUE,
#                              verbose = TRUE)

# Do separately.
for(i in 1:5){
  na_index <- which(apply(DFtest_kidney[gene_index_kidney,which(cuttest_kidney == i)],1,sum) <10)
  set.seed(1)
  sce_kidney_fit <- cellassign(exprs_obj = sce_kidney[-na_index,which(cuttest_kidney == i)], 
                               marker_gene_info = as.matrix(bigGM_kidney[which(rownames(bigGM_kidney) %in% rownames(sce_kidney[-na_index,])),]), 
                               s = sce_kidney_Size_factors[which(cuttest_kidney == i)], 
                               shrinkage = TRUE,
                               verbose = TRUE)
  gc()
  Y_kidney_cellassign[which(cuttest_kidney == i)] <- celltypes(sce_kidney_fit)
}
length(which(Y_kidney.raw == Y_kidney_cellassign))




### singleR
rownames(bigGM_kidney) <- rownames(bigGM_kidney)
celltypes_kidney <- colnames(bigGM_kidney)
cellmaintypes_kidney <- colnames(bigGM_kidney)
kidney_ref <- list(name = 'CellMatch',data = as.matrix(bigGM_kidney),
                   types = celltypes_kidney,main_types = cellmaintypes_kidney)
kidney_ref$de.genes <- CreateVariableGeneSet(as.matrix(bigGM_kidney),celltypes_kidney,200)
kidney_ref$de.genes.main <- CreateVariableGeneSet(as.matrix(bigGM_kidney),cellmaintypes_kidney,300)
rm(celltypes_kidney,cellmaintypes_kidney)

kidney_singleR <- CreateSinglerObject(counts = as.matrix(DFtest_kidney),
                                      project.name = 'kidney',
                                      ref.list = list(kidney_ref))

kidney_singleR <- kidney_singleR$singler
kidney_singleR1 <- kidney_singleR[[1]]
kidney_singleR1 <- kidney_singleR1$SingleR.single.main
Y_kidney_singleR <- NULL
Y_kidney_singleR$cellmatch <- kidney_singleR1$labels1
length(which(Y_kidney_singleR$cellmatch == Y_kidney.raw))