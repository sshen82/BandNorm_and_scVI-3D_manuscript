library(data.table)
library(foreach)
library(doParallel)

generateHiClusterInput = function(input_folder, output_folder, chr, res){
  # both input_folder and output_folder should be full paths containing all data.
  # chr should be a vector of the format of "chr" + i, e.g. chr3.
  # res means resolution.
  files = list.files(input_folder, full.names = TRUE, recursive = TRUE)
  names = list.files(input_folder, full.names = FALSE, recursive = TRUE, include.dirs = FALSE)
  for (c in chr){
    dir.create(paste0(output_folder, "/", c))
  }
  cl <- makeCluster(4)
  registerDoParallel(cl)
  output = foreach(i=1:length(input_folder), .packages=c("dplyr", "data.table", "matrixStats")) %dopar% {
    options(scipen = 200)
    cell = fread(input_folder[i])
    cell = cell[cell$V1 == cell$V3, ]
    cell$V2 = floor(cell$V2 / res)
    cell$V4 = floor(cell$V4 / res)
    temp = c()
    for (j in 1:length(chr)){
      cell2 = cell[cell$V1 == chr[j], ]
      if (nrow(cell2) != 0){
        write.table(cell2[, c(2, 4, 5)], file = paste(output_folder, "/", chr[j], "/",
                                        names[i], "_", chr[j], ".txt", sep = ""), 
                    row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      }
    }
  }
}

runHiCluster = function(names, input_folder, output_folder, chr, res, chr_size){
  # names is a vector of all the names for files
  # both input_folder and output_folder should be full paths containing all data.
  # chr should be a vector of the format of "chr" + i, e.g. chr3.
  # res means resolution.
  # chr_size is a path to the file containing the chromosome sizes
  for (i in 1:length(unique(chrs))){
    dir.create(path = paste(output_folder, "/", unique(chr)[i], sep = ""))
  }
  indirs = paste0(input_folder, "/", chr, "/")
  outdirs = paste0(output_folder, "/", chr, "/")
  cl <- makeCluster(4)
  registerDoParallel(cl)
  foreach(i=1:length(names), .packages="data.table") %dopar% {
    for (j in 1:length(chr)){
      runHiCluster <- sprintf(
        "hicluster impute-cell --indir %s --outdir %s --cell %s --chrom %s --res %s --chrom_file %s",
        indirs[j], outdirs[j], names[i], chrs[j], res, chr_size)
      system(runHiCluster)
    }
  }
  stopCluster(cl)
}
