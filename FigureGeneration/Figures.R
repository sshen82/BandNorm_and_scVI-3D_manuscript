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





