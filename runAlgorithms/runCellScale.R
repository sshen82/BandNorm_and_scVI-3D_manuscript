#This is the code used to generate cellscale embedding.
library(BandNorm)
input_path = ""
save_path = ""

cellscale = function(path = NULL, hic_df = NULL, save = TRUE, save_path = NULL) {
  if (is.null(path) && is.null(hic_df)){
    stop("Please specify one of path or hic_df.")
  }
  if (!is.null(path) && !is.null(hic_df)){
    warning("Specified both path and hic_df. Will use path file as default for saving memory.")
  }
  if (save) {
    if (is.null(save_path)){
      stop("path for saving the data is NULL!")
    }
    if (!dir.exists(save_path)){
      warning("path for saving the data doesn't exist, will create one according to the directory.")
      dir.create(save_path, recursive = TRUE)
    }
  }
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
  band_info <- hic_df %>% group_by(chrom, cell) %>% summarise(band_depth = sum(count))
  
  hic_df <- hic_df %>% left_join(band_info, by = c("chrom", "cell")) %>% mutate(BandNorm = count/band_depth) %>% 
    select(-c(band_depth, count))
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

output = cellscale(input_path, save = TRUE, save_path = save_path)
