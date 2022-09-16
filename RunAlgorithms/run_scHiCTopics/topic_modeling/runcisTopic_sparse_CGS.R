#############################################
# Run cisTopic
#############################################

library(cisTopic)
library(data.table)
library(Matrix)

mat_file <- "cs_mat.csv" # cell-LP matrix in sparse matrix format
lp_label_file <- "human_10_1000000_LPnames.txt" # LP annotations (i.e. LP names)
out_dir <- "cistopic_output" # out directory
cell_file <- "summary.txt" #summary file
low_b <- 10 # lowest topic number to try
up_b <- 90 # highest topic number to try
increment <- 5 # increment for topics
resolution <- 1000000 # resolution of contact matrix for name purposes

options(scipen = 10) 
topic_list = seq(low_b, up_b, by = increment)

df = fread(
    mat_file,
    col.names = c("cell.idx", "lp.idx", "count"), header = TRUE)

lp.annotations = read.table(
    lp_label_file,
    col.names = c("lp"),
    colClasses = c("character"))

cell.annotations = read.table(
    cell_file, header = TRUE)

rownames(lp.annotations) = lp.annotations$lp
rownames(cell.annotations) = cell.annotations$name

df$lp.idx = df$lp.idx + 1

# add a dummy cell to ensure that all genes are included in the matrix
# even if a gene isn't expressed in any cell
df = rbind(df, data.frame(
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    lp.idx =  c(1, nrow(lp.annotations)),
    count = c(1, 1)))

temp_idx = match(df$cell.idx, unique(df$cell.idx))

## make sparse matrix
mat = sparseMatrix(i = df$lp.idx + 1, j = temp_idx, x = df$count)
mat = mat[, 1:(ncol(mat)-1)]
mat = mat[1:(nrow(mat)-1), ]

rownames(mat) = lp.annotations$lp
colnames(mat) = cell.annotations$cell

cell_list = cell.annotations$cell_type
idx = order(cell_list)
mat <- mat[,idx]
cell_list <- as.matrix(cell_list[idx])
rownames(cell_list) <- colnames(mat)
colnames(cell_list) = "Cell type"


###
cisTopicObject <- createcisTopicObject(mat, project.name="full_cisTopic", keepCountsMatrix=FALSE)
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = cell_list)
cisTopicObject <- runCGSModels(cisTopicObject, topic=topic_list, seed=999, nCores=length(topic_list))
cisTopicObject = selectModel(cisTopicObject, type='maximum')
saveRDS(cisTopicObject, file="./output/cisTopicObject_CGS.rds")


