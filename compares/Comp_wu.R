# Wu comparisons

Wu_1 <- read.table("GSE103976_MeA_AllCells_DGE.txt", 
                   header = T, na.strings = c("NA"), nrows = 10000)
Wu_2 <- read.table("GSE103976_MeA_AllCells_DGE.txt", 
                   header = F,na.strings = c("NA"), skip = 10000)

dimnames(Wu_2)[[1]] <- Wu_1[,1]
Wu_2 <- Wu_2[,-1]
dimnames(Wu_2)[[2]] <- dimnames(Wu_1)[[2]]

Y_wu <- strsplit(dimnames(Wu_5)[[2]],'_') %>% unlist %>% matrix(ncol = 3, byrow = T)
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

GM_wu <- read.csv("cellmatch_wu.csv")
Wu_3 <- rbind(Wu_1,Wu_2)
genes_wu <- dimnames(Wu_3)[[1]]
genes_wu <- toupper(genes_wu)
index_wu <- which(genes_wu %in% GM_wu[,1])


DFtest_wu <- as.data.frame(Wu_3)
rm(Wu_1,Wu_2)
DFtest_wu <- Seurat::LogNormalize(data = DFtest_wu, 
                                  scale.factor = 10000, 
                                  verbose = F)
DFtest_wu <- as.data.frame(DFtest_wu)
DFtest_wu <- DFtest_wu[which(!duplicated(rownames(DFtest_wu))),]

DFtest_wu <- as.matrix(DFtest_wu)
DFtest_wu <- as(DFtest_wu,'dgCMatrix')
set.seed(1)
cut_wu <- cut(sample(1:length(Y_wu)),breaks = 5,labels = F)

DFtest_wu1 <- DFtest_wu[,which(cut_wu == 1)]
Y_wu1 <- Y_wu[which(cut_wu == 1)]
Y_wu1.raw <- Y_wu.raw[which(cut_wu == 1)]
na_index <- which(apply(DFtest_wu1,1,sum) ==0)
DFtest_wu1 <- DFtest_wu1[-na_index,]
DFtest_wu1 <- DFtest_wu1[!duplicated(rownames(DFtest_wu1)), ]
rownames(DFtest_wu1) <- toupper(rownames(DFtest_wu1))


DFtest_wu2 <- DFtest_wu[,which(cut_wu == 2)]
Y_wu2 <- Y_wu[which(cut_wu == 2)]
Y_wu2.raw <- Y_wu.raw[which(cut_wu == 2)]
na_index <- which(apply(DFtest_wu2,1,sum) == 0)
DFtest_wu2 <- DFtest_wu2[-na_index,]
DFtest_wu2 <- DFtest_wu2[!duplicated(rownames(DFtest_wu2)), ]
rownames(DFtest_wu2) <- toupper(rownames(DFtest_wu2))


DFtest_wu3 <- DFtest_wu[,which(cut_wu == 3)]
Y_wu3 <- Y_wu[which(cut_wu == 3)]
Y_wu3.raw <- Y_wu.raw[which(cut_wu == 3)]
na_index <- which(apply(DFtest_wu3,1,sum) == 0)
DFtest_wu3 <- DFtest_wu3[-na_index,]
DFtest_wu3 <- DFtest_wu3[!duplicated(rownames(DFtest_wu3)), ]
rownames(DFtest_wu3) <- toupper(rownames(DFtest_wu3))


DFtest_wu4 <- DFtest_wu[,which(cut_wu == 4)]
Y_wu4 <- Y_wu[which(cut_wu == 4)]
Y_wu4.raw <- Y_wu.raw[which(cut_wu == 4)]
na_index <- which(apply(DFtest_wu4,1,sum) == 0)
DFtest_wu4 <- DFtest_wu4[-na_index,]
DFtest_wu4 <- DFtest_wu4[!duplicated(rownames(DFtest_wu4)), ]
rownames(DFtest_wu4) <- toupper(rownames(DFtest_wu4))


DFtest_wu5 <- DFtest_wu[,which(cut_wu == 5)]
Y_wu5 <- Y_wu[which(cut_wu == 5)]
Y_wu5.raw <- Y_wu.raw[which(cut_wu == 5)]
na_index <- which(apply(DFtest_wu5,1,sum) == 0)
DFtest_wu5 <- DFtest_wu5[-na_index,]
DFtest_wu5 <- DFtest_wu5[!duplicated(rownames(DFtest_wu5)), ]
rownames(DFtest_wu5) <- toupper(rownames(DFtest_wu5))

# Codes are the same for all subset. Only subset Wu1 is presented here.

### Wu1
markers_wu <- 'markers_wu.txt'
bigGM_wu <- read.csv('bigGM_wu.csv', header = T, check.names=FALSE)
rownames(bigGM_wu) <- bigGM_wu[,1]
bigGM_wu <- bigGM_wu[,-1]
colnames(bigGM_wu) <- unique(GM_wu[,2])
bigGM_wu <- bigGM_wu[which(rownames(bigGM_wu) %in% genes_wu),]


pdata_wu1 <- data.frame(FACS_type = Y_wu1.raw)
rownames(pdata_wu1) <- colnames(DFtest_wu1)
fdata_wu1 <- data.frame(gene_short_name = rownames(DFtest_wu1))
rownames(DFtest_wu1) <- str_to_title(rownames(DFtest_wu1))
rownames(fdata_wu1) <- rownames(DFtest_wu1)
cds_wu1 <- new_cell_data_set(DFtest_wu1,
                             cell_metadata = pdata_wu1,
                             gene_metadata = fdata_wu1)
rownames(DFtest_wu1) <- toupper(rownames(DFtest_wu1))

### garnett
set.seed(1)
wu1_classifier <- train_cell_classifier(cds = cds_wu1,
                                        marker_file = markers_wu,
                                        db = org.Mm.eg.db,
                                        cds_gene_id_type = "SYMBOL",
                                        num_unknown = 50,
                                        marker_file_gene_id_type = "SYMBOL")
set.seed(1)
cds1_wu1 <- classify_cells(cds_wu1, wu1_classifier,
                           db = org.Mm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")

length(which(pData(cds1_wu1)$cluster_ext_type == pData(cds1_wu1)$FACS_type))



### cellassign
# As error in cellassign occurs, we change into raw counts form here.
# DFraw_wu <- as.data.frame(Wu_5)
# DFraw_wu <- DFraw_wu[which(!duplicated(rownames(DFraw_wu))),]
# 
# DFraw_wu <- as.matrix(DFraw_wu)
# DFraw_wu <- as(DFraw_wu,'dgCMatrix')
# set.seed(1)
# cut_wu <- cut(sample(1:length(Y_wu)),breaks = 5,labels = F)
# 
# DFraw_wu1 <- DFraw_wu[,which(cut_wu == 1)]
# Y_wu1 <- Y_wu[which(cut_wu == 1)]
# Y_wu1.raw <- Y_wu.raw[which(cut_wu == 1)]
# na_index <- which(apply(DFraw_wu1,1,sum) ==0)
# DFraw_wu1 <- DFraw_wu1[-na_index,]
# DFraw_wu1 <- DFraw_wu1[!duplicated(rownames(DFraw_wu1)), ]
# rownames(DFraw_wu1) <- toupper(rownames(DFraw_wu1))
sce_wu1 <- SingleCellExperiment(assays = list(counts = DFraw_wu1),
                                colData = Y_wu1.raw,
                                rowData = rownames(DFraw_wu1))
gene_index_wu1 <- which(rownames(DFraw_wu1) %in% rownames(bigGM_wu))
sce_wu1_Size_factors <- edgeR::calcNormFactors(sce_wu1)
sce_wu1 <- sce_wu1[gene_index_wu1,]

Y_wu1_cellassign <- rep(0, length(Y_wu1.raw))
set.seed(1)
cuttest_wu1 <- cut(sample(1:length(Y_wu1.raw)),breaks = 10,labels = F)
for(i in 1:10){
  na_index <- which(apply(DFraw_wu1[gene_index_wu1,which(cuttest_wu1 == i)],1,sum) < 10)
  set.seed(1)
  sce_wu1_fit <- cellassign(exprs_obj = sce_wu1[-na_index,which(cuttest_wu1 == i)], 
                            marker_gene_info = as.matrix(bigGM_wu[which(rownames(bigGM_wu) %in% rownames(sce_wu1[-na_index,])),]), 
                            s = sce_wu1_Size_factors$samples$norm.factors[which(cuttest_wu1 == i)], 
                            shrinkage = TRUE,
                            verbose = TRUE)
  gc()
  Y_wu1_cellassign[which(cuttest_wu1 == i)] <- celltypes(sce_wu1_fit)
}
length(which(Y_wu1.raw == Y_wu1_cellassign))
rm(sce_wu1)


### singleR
celltypes_wu1 <- colnames(bigGM_wu)
cellmaintypes_wu1 <- colnames(bigGM_wu)
wu1_ref <- list(name = 'CellMatch',data = as.matrix(bigGM_wu),
                types = celltypes_wu1,main_types = cellmaintypes_wu1)
wu1_ref$de.genes <- CreateVariableGeneSet(as.matrix(bigGM_wu),celltypes_wu1,200)
wu1_ref$de.genes.main <- CreateVariableGeneSet(as.matrix(bigGM_wu),cellmaintypes_wu1,300)
rm(celltypes_wu1,cellmaintypes_wu1)

wu1_singleR <- CreateSinglerObject(counts = as.matrix(DFtest_wu1),
                                   project.name = 'wu1',
                                   species = 'Mouse',
                                   ref.list = list(wu1_ref))

wu1_singleR <- wu1_singleR$singler
wu1_singleR1 <- wu1_singleR[[1]]
wu1_singleR1 <- wu1_singleR1$SingleR.single.main
Y_wu1_singleR <- NULL
Y_wu1_singleR$cellmatch <- wu1_singleR1$labels1

length(which(Y_wu1_singleR$cellmatch == Y_wu1.raw))
rm(wu1_singleR1)