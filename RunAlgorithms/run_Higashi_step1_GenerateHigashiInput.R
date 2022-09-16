library("dplyr")
library("Matrix")
library("mgcv")
library("dryhic")
library(data.table)
options(scipen = 7)

path = ""
output = ""
setwd(path)
cell_types = list.files()
cell_type_folders = list.files()
cell_type_folders = paste(path, cell_type_folders, sep = "")
num_types = length(cell_type_folders)

#Getting the path for each cell, and the clustering info.
cell_folders = c()
clusters = c()
cell_names = c()
for (i in 1:num_types){
  path = list.files(cell_type_folders[i])
  cell_names = c(cell_names, path)
  path = paste(cell_type_folders[i], "/", path, sep = "")
  cell_folders = c(cell_folders, path)
  clusters = c(clusters, rep(cell_types[i], length(path)))
}

#This if condition is just to check if you choose Ren's data or Ecker's.
num_chrs = 24
num_cells = length(cell_names)
chr_order = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
              "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
chrs = paste("chr", chr_order, sep = "")

#Loading cells
hic = c()
coverage = rep(0, num_cells)
for (j in 1:num_cells){
  hic_file = fread(cell_folders[j])
  # hic_file$V1 = paste("chr", hic_file$V1, sep = "")
  # hic_file$V3 = paste("chr", hic_file$V3, sep = "")
  # coverage[j] = sum(hic_file$V5)
  hic_file = cbind(hic_file, (j - 1))
  hic_file = cbind(hic_file, cell_names[j])
  hic_file = hic_file[, c(7, 6, 1, 2, 3, 4, 5)]
  hic = rbind(hic, hic_file)
}
# summaries = data.frame(cluster, coverage, libs)
# write.table(summaries, file = "C:/Users/solei/Downloads/Duan2017/summary_duan.txt", row.names = FALSE,
#             col.names = FALSE, quote = FALSE, sep = "\t")
colnames(hic) = c('cell_name','cell_id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'count')
write.table(hic, file = paste0(output, "data.txt"), row.names = FALSE, 
            quote = FALSE, sep = "\t")
