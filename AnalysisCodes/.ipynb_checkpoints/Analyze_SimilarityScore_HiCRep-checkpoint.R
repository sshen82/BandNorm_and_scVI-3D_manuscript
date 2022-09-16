library(hicrep)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
cT <- args[1]
method <- args[2]
cT2 <- args[3]
resolution <- 1000000

## True Label
s1Path <- paste0("/Bulk/ContactMatrix/", cT, "/Bulk/")
s2Path <- paste0("/Duan2020/ContactMatrix/", cT, "/", method, "/")
outPath <- paste0("/Duan2020/HiCRep/")

chrList <- paste0("chr", c(1:22, "X"))
reprod <- matrix(0, 1, length(chrList))
for(c in 1:length(chrList)){
    chr <- chrList[c]
    s1M <- fread(paste(s1Path, "/", chr, ".nij.matrix", sep = ""), header = F, fill = T)
    s2M <- fread(paste(s2Path, "/", chr, ".nij.matrix", sep = ""), head = F, fill = T)
    d <- min(c(dim(s1M), dim(s2M)))
    upperB <-resolution * (d - 1)
    print(d)
    s1M <- s1M[1:d, 1:d] %>% as.matrix
    s2M <- s2M[1:d, 1:d] %>% as.matrix

    diag(s1M) <- 0
    diag(s2M) <- 0
    
    ## print(sum(is.na(s1M)))
    ## print(sum(is.na(s2M)))

    ## s1M[is.na(s1M)] <- 0
    ## s2M[is.na(s2M)] <- 0
    ## s1M <- s1M/sum(s1M) * 2000000
    ## s2M <- s2M/sum(s2M) * 2000000

    ## if(method %in% c("3DVI", "higashi", "schicluster")){
    ##     h <- 0
    ## }else{
    ##     h <- 1
    ## }
    h <- 1
    reprod[1, c] <- get.scc(s1M, s2M, resolution, h, lbr = 0, ubr = upperB)$scc
    print(reprod[1, c])
}
write.table(reprod, file = paste(outPath, "/", cT, "_", method, ".scc", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Cluster Label
if(cT != "IMR90"){
    s1Path <- paste0("/Bulk/ContactMatrix/", cT, "/Bulk/")
    s2Path <- paste0("/Duan2020/ContactMatrix/", cT, "_ClusterLabel/", method, "/")
    outPath <- paste0("/Duan2020/HiCRep/")


chrList <- paste0("chr", c(1:22, "X"))
reprod <- matrix(0, 1, length(chrList))
for(c in 1:length(chrList)){
    chr <- chrList[c]
    s1M <- fread(paste(s1Path, "/", chr, ".nij.matrix", sep = ""), header = F, fill = T)
    s2M <- fread(paste(s2Path, "/", chr, ".nij.matrix", sep = ""), head = F, fill = T)
    d <- min(c(dim(s1M), dim(s2M)))
    upperB <-resolution * (d - 1)
    print(d)
    s1M <- s1M[1:d, 1:d] %>% as.matrix
    s2M <- s2M[1:d, 1:d] %>% as.matrix

    diag(s1M) <- 0
    diag(s2M) <- 0
    
    ## print(sum(is.na(s1M)))
    ## print(sum(is.na(s2M)))

    ## s1M[is.na(s1M)] <- 0
    ## s2M[is.na(s2M)] <- 0
    ## s1M <- s1M/sum(s1M) * 2000000
    ## s2M <- s2M/sum(s2M) * 2000000
    ## if(method %in% c("3DVI", "higashi", "schicluster")){
    ##     h <- 0
    ## }else{
    ##     h <- 1
    ## }
    h <- 1
    reprod[1, c] <- get.scc(s1M, s2M, resolution, h, lbr = 0, ubr = upperB)$scc

    print(reprod[1, c])
}
write.table(reprod, file = paste(outPath, "/", cT, "_ClusterLabel_", method, ".scc", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}


## same method between cell line
if(method != "Bulk"){
    s1Path <- paste0("/Duan2020/ContactMatrix/", cT, "/", method, "/")
    s2Path <- paste0("/Duan2020/ContactMatrix/", cT2, "/", method, "/")
}else{
    s1Path <- paste0("/Bulk/ContactMatrix/", cT, "/Bulk/")
    s2Path <- paste0("/Bulk/ContactMatrix/", cT2, "/Bulk/")

}
outPath <- paste0("/Duan2020/HiCRep/")

print(cT)
print(cT2)
print(method)

chrList <- paste0("chr", c(1:22, "X"))
reprod <- matrix(0, 1, length(chrList))
for(c in 1:length(chrList)){
    chr <- chrList[c]
    s1M <- fread(paste(s1Path, "/", chr, ".nij.matrix", sep = ""), header = F, fill = T)
    s2M <- fread(paste(s2Path, "/", chr, ".nij.matrix", sep = ""), head = F, fill = T)
    d <- min(c(dim(s1M), dim(s2M)))
    upperB <-resolution * (d - 1)
    print(d)
    s1M <- s1M[1:d, 1:d] %>% as.matrix
    s2M <- s2M[1:d, 1:d] %>% as.matrix
    
    diag(s1M) <- 0
    diag(s2M) <- 0
    
    h <- 1
    lbr <- 0
    reprod[1, c] <- get.scc(s1M, s2M, resolution, h, lbr = lbr, ubr = upperB)$scc
    print(reprod[1, c])
}
write.table(reprod, file = paste(outPath, "/", cT, "_", cT2, "_", method, ".scc", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

