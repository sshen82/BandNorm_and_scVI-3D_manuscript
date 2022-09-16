#Functions to generate evaluation metrics
library(data.table)
library(umap)
library(gridExtra)
library(Rtsne)
library(mclust)
library(cluster)
library(lisi)
library(stringr)
library(harmony)
library(cisTopic)
library(SingleCellExperiment)
library(Seurat)
library(pdfCluster)
require(RANN)
library(kBET)

mlgFindNeighbors = function(latent){
  rownames(latent) = paste0("cell", 1:nrow(latent))
  NN_graph = FindNeighbors(latent, verbose = FALSE, nn.method = "rann")
  G = NN_graph$snn
  return(G)
}

mlgARI = function(embed, summary){
  cluster = summary$cell_type
  clusterN = length(unique(cluster))
  G = mlgFindNeighbors(embed)
  l = 0.1
  h = 2
  iter = 1
  while (l < h) {
    mid = l + (h - l) / 2
    mlgc = try(FindClusters(G, resolution = mid, random.seed = 2020, verbose = FALSE)[[1]])
    if(class(mlgc) == "try-error"){
      h = mid
    }else{
      nml = length(unique(mlgc))
      if (nml == clusterN) break
      if (iter > 40) {
        print("not arrive")
        break
      }
      if (nml < clusterN) l = mid
      else h = mid 
    }
    iter = iter + 1
  }
  mlgARI = adj.rand.index(mlgc, cluster)
  return(mlgARI)
}

evaluation_func = function(embed, summary){
  summary$cell_type = as.factor(summary$cell_type)
  umap_out = umap(embed)
  tsne_out = Rtsne(embed)
  num_centers = length(unique(summary$cell_type))
  kmeansRES = kmeans(embed, centers = num_centers, nstart = 100)$cluster
  kmeansRESumap = kmeans(umap_out$layout, centers = num_centers, nstart = 100)$cluster
  kmeansREStsne = kmeans(tsne_out$Y, centers = num_centers, nstart = 100)$cluster
  return(data.frame(metric = c("Silhouette score + tSNE", "Silhouette score + UMAP",
                               "ARI + Embedding", "ARI + UMAP", "ARI + tSNE", "Batch ARI",
                               "MLG ARI"),
                    score = c(mean(silhouette(as.numeric(summary$cell_type), dist(data.frame(tsne_out$Y)))[, 3]),
                              mean(silhouette(as.numeric(summary$cell_type), dist(data.frame(umap_out$layout)))[, 3]),
                              adjustedRandIndex(kmeansRES, summary$cell_type),
                              adjustedRandIndex(kmeansRESumap, summary$cell_type),
                              adjustedRandIndex(kmeansREStsne, summary$cell_type),
                              ifelse("batch" %in% colnames(summary), adjustedRandIndex(kmeansRES, summary$batch), NA),
                              mlgARI(embed, summary))))
}

evaluationLISI = function(embed, summary){
  if ("batch" %in% colnames(summary)) {
    return(compute_lisi(embed, summary, 'batch')$batch)
  } else {
    return(rep(NA, nrow(summary)))
  }
}

#Function to generate UMAPs
UMAP = function(embed, summary, title_Graph){
  # title_Graph is the title for the graph.
  umap_gg = data.frame(umap(embed)$layout, batch = as.factor(summary$batch), cell = summary$cell_type, 
                       depth = summary$depth,
                       sparsity = summary$sparsity)
  
  ggplot(umap_gg, aes(x = X1, y = X2)) + geom_point(aes(color = cell)) + 
    ylab("UMAP 2") + xlab("UMAP 1") + labs(color = "Cell Type") + theme_bw(base_size = 15) + 
    ggtitle(title_Graph) + scale_color_viridis(option="magma", discrete = TRUE)
}

#Progressive Pooling for BandNorm
library(BandNorm)
bandnormPP = function(path = NULL, hic_df = NULL, save = TRUE, save_path = NULL, res) {
  # Get path and name for all the cells in this path.
  if (!is.null(path)){
    if (!dir.exists(path)){
      stop("path for data doesn't exist!")
    }
    paths = list.files(path, recursive = TRUE, full.names = TRUE)
    names = basename(list.files(path, recursive = TRUE))
    # The input format of the cell should be [chr1, bin1, chr2, bin2, count].
    load_cell = function(i) {
      return(fread(paths[i]) %>% rename(chrom = V1, binA = V2, binB = V4, count = V5) %>%
               mutate(diag = abs(binB - binA), cell = names[i]) %>% select(-V3))
    }
    hic_df = rbindlist(lapply(1:length(paths), load_cell))
  }
  # Calculate band depth and the mean of band depth for bandnorm.
  hic_df = hic_df[diag != 0]
  hic_df$diag = hic_df$diag / res
  hic_df = hic_df[hic_df$diag <= 20, ]
  hic_df$merge = floor((sqrt(8 * (hic_df$diag - 1) + 1) - 1) / 2)
  
  band_info <- hic_df %>% group_by(chrom, merge, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(chrom, merge) %>% summarise(depth = mean(band_depth))
  
  hic_df <- hic_df %>% left_join(alpha_j, by = c("chrom", "merge")) %>% left_join(band_info,
                                                                                  by = c("chrom", "merge", "cell")) %>% mutate(BandNorm = count/band_depth * depth) %>% 
    select(-c(band_depth, depth, count, merge))
  
  if (save) {
    save_path = file.path(save_path, unique(names))
    for (i in 1:length(unique(names))) {
      write.table(hic_df[cell == unique(names)[i], c("chrom", "binA", "chrom",
                                                     "binB", "BandNorm")], file = save_path[i], col.names = FALSE, row.names = FALSE,
                  quote = FALSE, sep = "\t")
    }
  }
  return(hic_df)
}

#Genomic Distance Effect
library(dplyr)
library(Matrix)
library(data.table)
library(umap)
library(pheatmap)
library(ggplot2)
library(GGally)
library(Rtsne)
library(mclust)
library(viridis)
library(dryhic)
bandSummary <- function(band_list, oDepth){
  band_logCPM = list() #
  bandDepth = matrix(0, 4238, 250) # depths of each band at the cell level
  bandMean = matrix(0, 4238, 250)  # means of each band
  bandVar = matrix(0, 4238, 250)  # vars of each band
  cellDepth = matrix(0, 4238, 250)  # cell depth - will repeat
  
  for (j in 1:250){ #go over each band
    band_logCPM[[j]] = log((band_list[[j]]+0.5)/oDepth*1e6) #logCPM
    bandMean[, j] = apply(band_logCPM[[j]], 1, mean) #band mean for each cell
    bandDepth[, j] = apply(band_logCPM[[j]], 1, sum) #band depth for each cell
    bandVar[, j] = apply(band_logCPM[[j]], 1, var) #band depth for each cell
    cellDepth[, j] = oDepth
  }
  return(list(band_list = band_list, band_logCPM = band_logCPM, bandDepth = bandDepth, bandMean = bandMean, bandVar = bandVar, cellDepth  = cellDepth))
}

readInScHic <- function(path){
  setwd(path)
  cell_types = list.files()
  cell_type_folders = list.files()
  cell_type_folders = paste(path, cell_type_folders, sep = "")
  num_types = length(cell_type_folders)
  
  #Getting the path for each cell, and the clustering info.
  cell_folders = c()
  clusters = c()
  for (i in 1:num_types){
    path = list.files(cell_type_folders[i])
    path = paste(cell_type_folders[i], "/", path, sep = "")
    cell_folders = c(cell_folders, path)
    clusters = c(clusters, rep(cell_types[i], length(path)))
  }
  num_cells = 4238
  # Loading cells
  band_list = list()
  temp = 250
  for (i in 1:250){
    band_list[[i]] = matrix(0, nrow = 4238, ncol = temp)
    temp = temp - 1
  }
  sparsity1 = c()
  depth = c()
  for (j in 1:4238){
    hic_file = fread(cell_folders[j]) #1 corresponds to chr 1
    hic_file = hic_file[hic_file$V1 == "chr1", ]
    sc_hic = sparseMatrix(i = hic_file$V2 / 1000000, j = hic_file$V4 / 1000000, x = hic_file$V5, dims = c(250, 250), index1 = FALSE)
    sc_hic = symmetrize_matrix(sc_hic)
    sc_hic = as.matrix(sc_hic)
    sparsity1[j] = sum(sc_hic > 0) / 2
    depth[j] = sum(sc_hic[col(sc_hic) - row(sc_hic) > 1]) #take all the upper triangular entries (excluding the diagonal)
    for (i in 1:250){
      # median_value[j, i] = max(1, median(sc_hic[abs(col(sc_hic) - row(sc_hic)) == (i - 1)]))
      band_list[[i]][j, ] = sc_hic[col(sc_hic) - row(sc_hic) == (i - 1)]
    }
  }
  
  return(list(oDepth = depth, band_list = band_list, sparsity = sparsity1, clusters = clusters))
}
path = "/Ecker2019/human/Counts_1mb_schic/"
## Functions readInScHic and bandSummary are at the end of the file
schic <- readInScHic(path)
bandS <- bandSummary(schic$band_list, schic$oDepth)
cluster = schic$clusters

all.d <- NULL
cell_ind = 1
qs = quantile(schic$oDepth, probs = seq(0, 1, 0.2))
for (cellT in unique(cluster)){
  cell.depth <- schic$oDepth[cellT == cluster]
  numBand = 100 # Get all the bands
  for (cellN in 1:sum(cellT == cluster)){
    if (cell.depth[cellN]>= qs[5]){dType = 5}
    if (cell.depth[cellN]< qs[5] & cell.depth[cellN]>= qs[4]){dType = 4}
    if (cell.depth[cellN]< qs[4] & cell.depth[cellN]>= qs[3]){dType = 3}
    if (cell.depth[cellN]< qs[3] & cell.depth[cellN]>= qs[2]){dType = 2}
    if (cell.depth[cellN]< qs[2] ){dType = 1}
    d = data.frame(dist=1:250, meanCount = bandS$bandMean[cell_ind,])[2:numBand, ] #exclude diagonal
    
    # When using only 50 bands: span = 0.05 = no smoothing; 0.1 is good but still not monotone decreasing, 0.15 seems great
    # With 250 bands, span = 0.12 seems good
    # With 100 bands, span = 0.15 seems good
    fit = loess(meanCount ~ dist, degree=1, span = 0.15, data=d)
    d.fit <- data.frame(d, smooth = fit$fitted, cellDepth = cell.depth[cellN], cellNo = rep(cellN, length(fit$fitted)), depthG = rep(paste("sizeGroup: ", dType, sep = ""), length(fit$fitted)), cellT = rep(cellT,  length(fit$fitted)))
    all.d <- rbind(all.d, d.fit)
    p = d %>% mutate(smooth = fit$fitted) %>%
      ggplot(aes(dist, meanCount)) +
      geom_point(size = 3, alpha = .5, color = "grey") +
      geom_line(aes(dist, smooth), color="red")
    pTitle = paste("Depth = ", cell.depth[cellN], "Cell no: ", cellN, sep = "")
    p1 = p + labs(title = pTitle)
    # print(p1)
    cell_ind = cell_ind + 1
  }
}
all.d$depthG = factor(all.d$depthG, levels = paste("sizeGroup:", 1:5)) # check level ordering

## bands 2-51
ggplot(all.d[all.d$dist<=51, ], aes(dist, meanCount)) +
  geom_point(size = 1, alpha = .5, color = "gray") +
  geom_line(aes(dist, smooth, group = cellNo, color=cellDepth/10000))+
  facet_grid( depthG ~ cellT)+theme(legend.position="top")+scale_color_viridis(option = "D")+
  labs(color=expression(paste("Library size x", 10^5)))+xlab("Band Number")+ylab("Mean contact count (logCPM)")



