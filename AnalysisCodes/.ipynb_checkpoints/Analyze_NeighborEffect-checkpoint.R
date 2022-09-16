library(pacman)
p_load(dplyr, Matrix, mgcv, data.table, R.matlab, umap, viridis, ggpubr, gtools, ggplot2, Rtsne, tidyr, stringr, magrittr, cowplot, cluster, pdfCluster, Seurat)

args = commandArgs(trailingOnly=TRUE)
dataName = args[1]
chromSelect = args[2]
 

## 2D binary
getNeighbors <- function(data, topI){
    localNeiN = data.frame(distToTarget = factor(1:10), numNei = (1:10) *4)
    neiI <- data %>% filter(diag!=0) %>% filter(binA >=topI$binA - 5) %>% filter(binA <= topI$binA + 5) %>% filter(binB >= topI$binB - 5) %>% filter(binB <= topI$binB + 5)

    dataFreq = c(abs(neiI$binA - topI$binA) + abs(neiI$binB - topI$binB)) %>% table %>% data.frame
    colnames(dataFreq) <- c("distToTarget", "obsN")
    ## dataFreq$distToTarget <- dataFreq$distToTarget %>% as.numeric
    dataFreqNei <- left_join(localNeiN, dataFreq, by = "distToTarget") 
    dataFreqNei$obsN[which(is.na(dataFreqNei$obsN))] <- 0
    return(dataFreqNei %>% mutate(prob = obsN/numNei))
}


getNeighborsHigh <- function(data, topI){
    localNeiN = data.frame(distToTarget = factor(1:10), numNei = (1:10) *4)
    neiI <- data %>% filter(diag!=0, count > 2) %>% filter(binA >=topI$binA - 5) %>% filter(binA <= topI$binA + 5) %>% filter(binB >= topI$binB - 5) %>% filter(binB <= topI$binB + 5)

    dataFreq = c(abs(neiI$binA - topI$binA) + abs(neiI$binB - topI$binB)) %>% table %>% data.frame
    colnames(dataFreq) <- c("distToTarget", "obsN")
    ## dataFreq$distToTarget <- dataFreq$distToTarget %>% as.numeric
    dataFreqNei <- left_join(localNeiN, dataFreq, by = "distToTarget") 
    dataFreqNei$obsN[which(is.na(dataFreqNei$obsN))] <- 0
    return(dataFreqNei %>% mutate(prob = obsN/numNei))
}




csvPath = paste0("/scHiC/results/", dataName, "/CSV/")
rdsPath = paste0("/scHiC/results/", dataName, "/RDS/")
rdataPath = paste0("//scHiC/results/", dataName, "/RData/")
figPath = paste0("//scHiC/results/", dataName, "/figures/")
cellInfo = readRDS(paste0(rdsPath, "cellInfo.rds"))
cellNameList = readRDS(paste0(rdsPath, "cellNameList.rds"))

cellFilterList = c()
neiI <- c()
neiIhigh <- c()

print(chromSelect)
    
hic_df = readRDS(file = paste0(rdsPath, "hic_df_", chromSelect, ".rds"))
    
## setting 1. top 10 bin-pairs per cell per chromxo
target = hic_df %>% filter(diag != 0, count > 2) %>% group_by(cell, batch, cluster) %>% top_n(10, wt = count) %>% ungroup
print(dim(target))

for(i in 1:nrow(target)){
    if(i %% 5000 == 0){
        print(i)
    }
    topI = target[i, ]
    neiI <- rbind(neiI, getNeighbors(hic_df %>% filter(cell == topI$cell), topI))
    neiIhigh <- rbind(neiIhigh, getNeighborsHigh(hic_df %>% filter(cell == topI$cell), topI))
}

save(neiI, neiIhigh, file = paste0(rdataPath, "neiI_", chromSelect, "_binary.Rdata"))
