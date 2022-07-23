# Wu

#This dataset is randomly split into 5 subsets to be analyzed.


### 1
Wu_1 <- read.table("GSE103976_MeA_AllCells_DGE.txt", 
                   header = T, na.strings = c("NA"), nrows = 10000)

Wu_2 <- read.table("GSE103976_MeA_AllCells_DGE.txt", 
                   header = F,na.strings = c("NA"), skip = 10000)

Wu_2 <- Wu_2[,-1]
dimnames(Wu_2)[[2]] <- dimnames(Wu_1)[[2]]
Y_wu <- strsplit(dimnames(Wu_2)[[2]],'_') %>% unlist %>% matrix(ncol = 3, byrow = T)
Y_wu <- Y_wu[,3]
unique(Y_wu)
Y_wu.raw <- Y_wu
Y_wu.raw[which(Y_wu.raw == 'OPC')] <- "Precursor Cell"
Y_wu.raw[which(Y_wu.raw == 'AS')] <- "Astrocyte"
Y_wu.raw[which(Y_wu.raw == 'OL')] <- "Oligodendrocyte"
Y_wu.raw[which(Y_wu.raw == 'N')] <- "Neuron"
Y_wu.raw[which(Y_wu.raw == 'MG')] <- "Microglial Cell"
Y_wu.raw[which(Y_wu.raw == 'EN')] <- "Endothelial Cell"
Y_wu.raw[which(Y_wu.raw == 'MU')] <- "Mural Cell"

Y_wu <- as.factor(Y_wu.raw) %>% as.numeric

class_wu <- data.frame(cluster=c(1:7),
                       type = c("Astrocyte","Endothelial Cell","Microglial Cell","Mural Cell",
                                "Neuron","Oligodendrocyte","Precursor Cell"))
#write.csv(class_wu, 'ann_wu.csv')

GM_wu <- read.csv("cellmatch_wu.csv")
Wu_3 <- rbind(Wu_1,Wu_2)
genes_wu <- dimnames(Wu_3)[[1]]
genes_wu <- toupper(genes_wu)
index_wu <- which(genes_wu %in% GM_wu[,1])
dimnames(DF_wu)[[2]] <- toupper(dimnames(DF_wu)[[2]])
DF_wu <- Wu_3[index_wu,]
rm(Wu_1,Wu_2,Wu_3)
DF_wu <- Seurat::LogNormalize(data = DF_wu, 
                              scale.factor = 10000, 
                              verbose = F)
DF_wu <- t(DF_wu)

#random split
set.seed(1)
cut_wu <- cut(sample(1:length(Y_wu)),breaks = 5,labels = F)

DF_wu1 <- DF_wu[which(cut_wu == 1),]
Y_wu1 <- Y_wu[which(cut_wu == 1)]
Y_wu1.raw <- Y_wu.raw[which(cut_wu == 1)]
na_index <- which(apply(DF_wu1,2,sum) == 0)
DF_wu1 <- DF_wu1[,-na_index]
dimnames(DF_wu1)[[2]] <- toupper(dimnames(DF_wu1)[[2]])

DF_wu2 <- DF_wu[which(cut_wu == 2),]
Y_wu2 <- Y_wu[which(cut_wu == 2)]
Y_wu2.raw <- Y_wu.raw[which(cut_wu == 2)]
#na_index <- which(apply(DF_wu2,2,sum) == 0)
#DF_wu2 <- DF_wu2[,-na_index]
dimnames(DF_wu2)[[2]] <- toupper(dimnames(DF_wu2)[[2]])

DF_wu3 <- DF_wu[which(cut_wu == 3),]
Y_wu3 <- Y_wu[which(cut_wu == 3)]
Y_wu3.raw <- Y_wu.raw[which(cut_wu == 3)]
na_index <- which(apply(DF_wu3,2,sum) == 0)
DF_wu3 <- DF_wu3[,-na_index]
dimnames(DF_wu3)[[2]] <- toupper(dimnames(DF_wu3)[[2]])

DF_wu4 <- DF_wu[which(cut_wu == 4),]
Y_wu4 <- Y_wu[which(cut_wu == 4)]
Y_wu4.raw <- Y_wu.raw[which(cut_wu == 4)]
na_index <- which(apply(DF_wu4,2,sum) == 0)
DF_wu4 <- DF_wu4[,-na_index]
dimnames(DF_wu4)[[2]] <- toupper(dimnames(DF_wu4)[[2]])

DF_wu5 <- DF_wu[which(cut_wu == 5),]
Y_wu5 <- Y_wu[which(cut_wu == 5)]
Y_wu5.raw <- Y_wu.raw[which(cut_wu == 5)]
#na_index <- which(apply(DF_wu5,2,sum) == 0)
#DF_wu5 <- DF_wu5[,-na_index]
dimnames(DF_wu5)[[2]] <- toupper(dimnames(DF_wu5)[[2]])


### 2
SELF_wu1 <- self_clustering(inputTags = t(DF_wu1), k_fixed = 15, SEED = 123)
wu1.comb <- CommunityMining(SELF_wu1)

SELF_wu2 <- self_clustering(inputTags = t(DF_wu2), k_fixed = 15, SEED = 123)
wu2.comb <- CommunityMining(SELF_wu2)

SELF_wu3 <- self_clustering(inputTags = t(DF_wu3), k_fixed = 15, SEED = 123)
wu3.comb <- CommunityMining(SELF_wu3)

SELF_wu4 <- self_clustering(inputTags = t(DF_wu4), k_fixed = 15, SEED = 123)
wu4.comb <- CommunityMining(SELF_wu4)

SELF_wu5 <- self_clustering(inputTags = t(DF_wu5), k_fixed = 15, SEED = 123)
wu5.comb <- CommunityMining(SELF_wu5)

### We only show codes for subset wu1. Others are the same.

min_num <- 10 
comb_vector_wu1 <- rep(0,length(unique(wu1.comb)))
for(i in 1:length((wu1.comb))){
  comb_vector_wu1[wu1.comb[i]] <- comb_vector_wu1[wu1.comb[i]] + 1
}
Y_min_in_wu1 <- which(comb_vector_wu1 >= min_num)
Y_min_out_wu1 <- which(comb_vector_wu1 < min_num)
Ycomb_wu1 <- wu1.comb
Ycomb_wu1[which(Ycomb_wu1 %in% Y_min_out_wu1)] <- 0
Ycomb_wu1_in <- which(Ycomb_wu1 != 0)
Ycomb_wu1_out <- which(Ycomb_wu1 == 0)

DFpca_wu1 <- stats::prcomp(DF_wu1, rank = 50)$x
set.seed(1)
umap_wu1 <- uwot::umap(DFpca_wu1, n_components = 10)

### 3
wu1_import <- colnames(DF_wu1) %>% toupper
label_wu1 <- NULL
Yclass_wu1 <- NULL
DFclass_wu1 <- NULL
M_wu1 <- NULL
X_wu1 <- NULL
IM_wu1 <- NULL
type_wu1 <- NULL
info_wu1 <- NULL
gain_wu1 <- NULL
Gain_wu1 <- NULL
Ycol_wu1 <- NULL

for(i in 1:length(Y_min_in_wu1)){
  DFclass_wu1 <- NULL
  Yclass_wu1[[i]] <- rep(0, length(Ycomb_wu1))
  Yclass_wu1[[i]][which(Ycomb_wu1 == Y_min_in_wu1[i])] <- 1
  
  Yclass_wu1[[i]] <- Yclass_wu1[[i]] %>% as.numeric
  partial_wu1 <- ceiling(length(Ycomb_wu1) / length(which(Ycomb_wu1 == Y_min_in_wu1[i])) ) 
  DFclass_wu1 <- cbind(Yclass_wu1[[i]], DF_wu1)
  
  if(partial_wu1 > 2){
    for(k in 1:(partial_wu1 - 2)){
      DFclass_wu1 <- rbind(DFclass_wu1, DFclass_wu1[which(Ycomb_wu1 == Y_min_in_wu1[i]),])
    }
  }
  
  M_wu1[[i]] <- xgb.DMatrix(data = as(as.matrix(DFclass_wu1[,-1]),"dgCMatrix"),label = DFclass_wu1[,1])
  X_wu1[[i]] <- xgboost(M_wu1[[i]], 
                        max.depth = 1, 
                        eta = 0.3, 
                        nround = 50, 
                        objective = 'binary:hinge', 
                        eval_metric = "auc")
  IM_wu1[[i]] <- xgb.importance(wu1_import, model = X_wu1[[i]])
  IM_wu1[[i]][,2] <- as.vector(IM_wu1[[i]][,2])

  for(k in 1:nrow(IM_wu1[[i]])){
    this_gene <- IM_wu1[[i]][k,1]
    this_mean <- mean(DF_wu1[which(Ycomb_wu1 == Y_min_in_wu1[i]),
                             which(as.factor(this_gene) == dimnames(DF_wu1)[[2]])])
    other_mean <- mean(DF_wu1[which(Ycomb_wu1 != Y_min_in_wu1[i]),
                              which(as.factor(this_gene) == dimnames(DF_wu1)[[2]])])
    if(other_mean > this_mean){
      IM_wu1[[i]][k,2] <- 0
    }
  }
  
  info_wu1[[i]] <- 0
  for(k in 1:nrow(IM_wu1[[i]])){
    info_wu1[[i]] <- c(info_wu1[[i]], which(as.data.frame(GM_wu)[,1]  == as.data.frame(IM_wu1[[i]])[k,1]))
  }
  info_wu1[[i]] <- info_wu1[[i]][-1]
  
  type_wu1[[i]] <- GM_wu[info_wu1[[i]],2]
  
  gain_wu1 <- NULL
  mid_gene <- NULL
  for(j in 1:length(unique(type_wu1[[i]]))){
    mid_gene <- which(GM_wu[,2] == unique(type_wu1[[i]])[j])
    mid_index <- which(as.data.frame(IM_wu1[[i]])[,1] %in% as.vector(GM_wu[mid_gene,1]))
    gain_wu1[j] <- sum(IM_wu1[[i]][mid_index,2])
    
  }
  Gain_wu1[[i]] <- data.frame(unique(type_wu1[[i]]), gain_wu1)
  Gain_wu1[[i]] <- Gain_wu1[[i]][order(Gain_wu1[[i]][,2],decreasing = T),]
  label_wu1[i] <- Gain_wu1[[i]][1,1]
  rm(gain_wu1)
}


### 4 
train.wu1 <- cbind(Ycomb_wu1, umap_wu1)[which(Ycomb_wu1 !=0), ] %>% as.data.frame()
test.wu1 <- cbind(Ycomb_wu1, umap_wu1)[which(Ycomb_wu1 ==0), ] %>% as.data.frame()
set.seed(1)
knn_wu1 <- kknn(Ycomb_wu1 ~., 
                train = train.wu1,
                test = test.wu1, 
                k = 1,
                kernel = 'rectangular',
)

Ynew_wu1 <- Ycomb_wu1
Ynew_wu1[which(Ycomb_wu1 == 0)] <- knn_wu1$fitted.values
label2_wu1 <- class_wu[match(label_wu1,class_wu[,2]),1]


Yfinal_wu1 <- rep(0,length(Y_wu1))
for(i in 1:length(Y_min_in_wu1)){
  Yfinal_wu1[which(Ynew_wu1 == Y_min_in_wu1[i])] <- label_wu1[i]
}
Yfinal_wu1[which(is.na(Yfinal_wu1))] <- 0
for(i in 1:nrow(class_wu)){
  print(sum(Yfinal_wu1[which(Y_wu1 == class_wu[i,1])] == class_wu[i,1]))
}
for(i in 1:nrow(class_wu)){
  print(sum(Y_wu1[which(Y_wu1 == class_wu[i,1])] == class_wu[i,1]))
}

### results
length(which(Y_wu1.raw == Yfinal_wu1)) 
