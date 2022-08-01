# kidney

# Codes are similar for all datasets. Thus we only shows codes for one dataset.
# For packages, see 'functions/packages.R'.

# input data

### 1 

ann <- read.table("annotations_droplet.csv",sep=",",header=TRUE)
ann[,2] <- as.factor(ann[,2])
ann_cluster <- ann[,2] %>% as.numeric
ann[,4] <- ann_cluster

GM_kidney <- read.csv("droplet\\Kidney-10X_P4_5\\cellmatch_kidney.csv")
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
index_kidney <- which(genes_kidney %in% GM_kidney[,1])

Matrix_kidney <- cbind(matrix_kidney, matrix_kidney2, matrix_kidney3)
Matrix_kidney <- Seurat::LogNormalize(data = Matrix_kidney, scale.factor = 10000, verbose = F)
Matrix_kidney <- Matrix_kidney[index_kidney, ]
Matrix_kidney <- as.matrix(t(Matrix_kidney))


DF_kidney <- data.frame(Matrix_kidney)
row.names(DF_kidney) <- B_kidney
names(DF_kidney) <- genes_kidney[index_kidney]

Y_kidney <- ann[match(rownames(DF_kidney), ann[,1]),][,4]
na.Y_kidney <- which(is.na(Y_kidney))
Y_kidney <- Y_kidney[-na.Y_kidney]
DF_kidney <- DF_kidney[-na.Y_kidney,]
rm(Matrix_kidney)
na_index <- which(apply(DF_kidney,2,sum) == 0)
DF_kidney <- DF_kidney[,-na_index]

class_kidney <- read.csv("ann_kidney.csv")
Y_kidney.raw <- class_kidney[match(Y_kidney, class_kidney[,1]),2]

# The proposed method: ACAM. Codes below are the same for all datasets.
# For Wu dataset, it is randomly split by the following function.
# set.seed(1)
# cut_wu <- cut(sample(1:length(Y_wu)),breaks = 5,labels = F)

### 2 
# For details, see 'functions/clustering.R'.
SELF_kidney <- individual_clustering(inputTags = t(DF_kidney), mt_filter = TRUE,
                                     percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
                                     nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, 
                                     save.results = FALSE, SEED = 123)
# It takes time to run clustering function. 
# Louvain results are saved in 'results' file for convenience.
kidney.comb <- CommunityMining(SELF_kidney)


min_num <- 10
comb_vector_kidney <- rep(0,length(unique(kidney.comb)))
for(i in 1:length((kidney.comb))){
  comb_vector_kidney[kidney.comb[i]] <- comb_vector_kidney[kidney.comb[i]] + 1
}

Y_min_in_kidney <- which(comb_vector_kidney >= min_num)
Y_min_out_kidney <- which(comb_vector_kidney < min_num)
Ycomb_kidney <- kidney.comb
Ycomb_kidney[which(Ycomb_kidney %in% Y_min_out_kidney)] <- 0
Ycomb_kidney_in <- which(Ycomb_kidney != 0)
Ycomb_kidney_out <- which(Ycomb_kidney == 0)


### 3
DFpca_kidney <- stats::prcomp(DF_kidney, rank = 50)$x
set.seed(1)
umap_kidney <- uwot::umap(DFpca_kidney, n_components = 10)


kidney_import <- colnames(DF_kidney) %>% toupper
label_kidney <- NULL
Yclass_kidney <- NULL
DFclass_kidney <- NULL
M_kidney <- NULL
X_kidney <- NULL
IM_kidney <- NULL
type_kidney <- NULL
info_kidney <- NULL
gain_kidney <- NULL
Gain_kidney <- NULL
Ycol_kidney <- NULL
for(i in 1:length(Y_min_in_kidney)){
  Yclass_kidney[[i]] <- rep(0, length(Ycomb_kidney))
  Yclass_kidney[[i]][which(Ycomb_kidney == Y_min_in_kidney[i])] <- 1
  
  Yclass_kidney[[i]] <- Yclass_kidney[[i]] %>% as.numeric
  partial_kidney <- ceiling(length(Ycomb_kidney) / length(which(Ycomb_kidney == Y_min_in_kidney[i])) ) 
  DFclass_kidney <- cbind(Yclass_kidney[[i]], DF_kidney)
  
  if(partial_kidney > 2){
    for(k in 1:(partial_kidney- 2)){
      DFclass_kidney <- rbind(DFclass_kidney, DFclass_kidney[which(Ycomb_kidney == Y_min_in_kidney[i]),])
    }
  }
  M_kidney[[i]] <- xgb.DMatrix(data = DFclass_kidney[,-1] %>% as.matrix,label = DFclass_kidney[,1])
  X_kidney[[i]] <- xgboost(M_kidney[[i]], 
                           max.depth = 1, 
                           nround = 50, 
                           objective = 'binary:hinge', 
                           eval_metric = "auc")
  IM_kidney[[i]] <- xgb.importance(kidney_import, model = X_kidney[[i]])
  IM_kidney[[i]][,2] <- as.vector(IM_kidney[[i]][,2])
  for(k in 1:nrow(IM_kidney[[i]])){
    this_gene <- IM_kidney[[i]][k,1]
    this_mean <- mean(DF_kidney[which(Ycomb_kidney == Y_min_in_kidney[i]),
                                which(as.factor(this_gene) == dimnames(DF_kidney)[[2]])])
    other_mean <- mean(DF_kidney[which(Ycomb_kidney != Y_min_in_kidney[i]),
                                 which(as.factor(this_gene) == dimnames(DF_kidney)[[2]])])
    if(other_mean > this_mean){
      IM_kidney[[i]][k,2] <- 0
    }
  }
  info_kidney[[i]] <- 0
  for(k in 1:nrow(IM_kidney[[i]])){
    info_kidney[[i]] <- c(info_kidney[[i]], which(as.data.frame(GM_kidney)[,1]  == as.data.frame(IM_kidney[[i]])[k,1]))
  }
  info_kidney[[i]] <- info_kidney[[i]][-1]
  
  type_kidney[[i]] <- GM_kidney[info_kidney[[i]],2]
  
  gain_kidney <- NULL
  mid_gene <- NULL
  for(j in 1:length(unique(type_kidney[[i]]))){
    mid_gene <- which(GM_kidney[,2] == unique(type_kidney[[i]])[j])
    mid_index <- which(as.data.frame(IM_kidney[[i]])[,1] %in% as.vector(GM_kidney[mid_gene,1]))
    gain_kidney[j] <- sum(IM_kidney[[i]][mid_index,2])
    
  }
  Gain_kidney[[i]] <- data.frame(unique(type_kidney[[i]]), gain_kidney)
  Gain_kidney[[i]] <- Gain_kidney[[i]][order(Gain_kidney[[i]][,2],decreasing = T),]
  label_kidney[i] <- Gain_kidney[[i]][1,1]
  rm(gain_kidney)
}

### 4 

train.kidney <- cbind(Ycomb_kidney,umap_kidney)[which(Ycomb_kidney !=0), ] %>% as.data.frame()
test.kidney <- cbind(Ycomb_kidney,umap_kidney)[which(Ycomb_kidney ==0), ] %>% as.data.frame()
set.seed(1)
knn_kidney <- kknn(Ycomb_kidney ~., 
                   train = train.kidney,
                   test = test.kidney, 
                   k = 1,
                   kernel = 'rectangular',
)

Ynew_kidney <- Ycomb_kidney
Ynew_kidney[which(Ycomb_kidney == 0)] <- knn_kidney$fitted.values

label2_kidney <- class_kidney[match(label_kidney,class_kidney[,2]),1]


Yfinal_kidney <- rep(0,length(Y_kidney))
for(i in 1:length(Y_min_in_kidney)){
  Yfinal_kidney[which(Ynew_kidney == Y_min_in_kidney[i])] <- label_kidney[i]
}

### results
length(which(Y_kidney.raw == Yfinal_kidney))