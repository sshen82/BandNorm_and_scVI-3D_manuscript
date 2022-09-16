library(pacman)
p_load(data.table)
p_load(dplyr)
p_load(HiCcompare)
p_load(ggplot2)
p_load(viridis)
p_load(diffHic)
p_load(edgeR)
p_load(InteractionSet)
p_load(csaw)
p_load(statmod)


diffhicF <- function(data, testN, resolution){

    rep1 <- "cond1"
    rep2 <- "cond2"
    ave.ab <- aveLogCPM(asDGEList(data))
    count.keep <- ave.ab >= aveLogCPM(1, lib.size=mean(data$totals)) ## 1 for 4vs4 5 for 2vs2
    print(aveLogCPM(1, lib.size=mean(data$totals)))
    data <- data[count.keep,]
    dist.keep <-  pairdist(data) > resolution
    data <- data[dist.keep, ]
    data <- normOffsets(data, se.out=TRUE)
    design <- model.matrix(~factor(c(rep(rep1, testN), rep(rep2, testN))))
    colnames(design) <- c("intercept", "replicate")
    y <- asDGEList(data)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    result <- glmQLFTest(fit, coef = 2)
    rowData(data) <- cbind(rowData(data), result$table)
    data <- data %>% interactions %>% data.frame %>% mutate(chr = 1, region1 = (start1 - 1 + end1)/2, region2 = (start2 - 1 + end2)/2, p.value = PValue, p.adj = p.adjust(PValue, method = "BH")) %>% select(chr, region1, region2, logFC, p.value, p.adj)
    
    return( data)
}

runDiffHic <- function(repList, bulkData1, bulkData2, mData1, mData2, resolution){

    repN <- length(repList)
    m1 <- vector("list", length(repList))
    m2 <- vector("list", length(repList))
    ## gSize <- 6231
    ## all.regions <- GRanges(chr, IRanges(as.numeric(0:(gSize - 1)*resolution+1), as.numeric(1:gSize*resolution)))
    to.keep <- 0
    minC <- 1
    resolution <- 1000000

    bulkM <- fread(bulkData1) %>% as.matrix
    methodM <- fread(mData1) %>% as.matrix
    gSize <- min(dim(bulkM), dim(methodM))
    all.regions <- GRanges(chr, IRanges(as.numeric(0:(gSize - 1)*resolution+1), as.numeric(1:gSize*resolution)))
    m1[[1]] <- bulkM[1:gSize, 1:gSize] %>% ContactMatrix(all.regions, all.regions)
    m2[[1]] <- methodM[1:gSize, 1:gSize] %>% ContactMatrix(all.regions, all.regions)
    to.keep <- as.matrix(m1[[1]])>=minC
    to.keep <- to.keep | as.matrix(m2[[1]])>=minC

    bulkM <- fread(bulkData2) %>% as.matrix
    methodM <- fread(mData2) %>% as.matrix
    m1[[2]] <- bulkM[1:gSize, 1:gSize] %>% ContactMatrix(all.regions, all.regions)
    m2[[2]] <- methodM[1:gSize, 1:gSize] %>% ContactMatrix(all.regions, all.regions)
    to.keep <- to.keep | as.matrix(m1[[2]])>=minC
    to.keep <- to.keep | as.matrix(m2[[2]])>=minC
    
    for(rep in 1:repN){
        iset <- deflate(m1[[rep]], extract=to.keep)
        m1[[rep]] <- inflate(iset, all.regions, all.regions)
        iset <- deflate(m2[[rep]], extract=to.keep)
        m2[[rep]] <- inflate(iset, all.regions, all.regions)
    }
    
    diffD <- mergeCMs(m1[[1]], m1[[2]],
                      m2[[1]], m2[[2]])
    res <- diffhicF(diffD, 2, resolution)
    return(res)
    
}

repList <- paste0("rep", 1:2)
repN <- length(repList)
resolution <- 1000000
chromList <- paste0("chr", c(1:22, "X"))
cellList <- c("GM12878", "H1ESC", "HAP1", "HFF", "IMR90")
methodList <- c("BandNorm", "3DVI", "Higashi", "scHiCluster")
summaryD <- c()
for(cT in cellList){
    print(cT)
    bulkPath <- paste0("/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Bulk/ContactMatrix/", cT)
    methodPath <- paste0("/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Duan2020/ContactMatrix/", cT)

    for(chr in chromList){
        print(chr)
        bulkData1 <- paste0(bulkPath, "/rep1/", chr, ".nij.matrix")
        bulkData2 <- paste0(bulkPath, "/rep2/", chr, ".nij.matrix")
    
        for(method in methodList){
            print(method)
            mData1 <- paste0(methodPath, "/", method, "_rep1/", chr, ".nij.matrix")
            mData2 <- paste0(methodPath, "/", method, "_rep2/", chr, ".nij.matrix")
        
            res <- runDiffHic(repList, bulkData1, bulkData2, mData1, mData2, resolution)

            summaryD <- data.frame(signN = c(res %>% filter(p.adj <= 0.001) %>% nrow(), res %>% filter(p.adj <= 0.01) %>% nrow(), res %>% filter(p.adj <= 0.05) %>% nrow(), res %>% filter(p.adj <= 0.1) %>% nrow()), padj = c(0.001, 0.01, 0.05, 0.1), method1 = "Bulk", method2 = method, chrom = chr, cellType = cT) %>% rbind(summaryD, .)

        }
        
    }
    

}
save(summaryD, file = "/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Duan2020/RData/diffHic_summaryD.RData")

summaryD$method2 <- factor(summaryD$method2, levels = methodList)
pdf("/p/keles/schic/volumeA/Figures/diffHic_BulkVSmethod.pdf", width = 20, height = 13)
summaryD %>% filter(padj >= 0.01) %>% ggplot(aes(x = method2, y = signN, fill = method2)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    facet_grid(padj~cellType, scale = "free") +
    theme_bw(base_size = 25) +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    xlab("") +
    ylab("# of Differentially Interacting Locus-pair") +
    ggpubr::rremove("legend") +
    ggpubr::rotate_x_text(angle = 20)
dev.off()

for(c1 in 1:(length(cellList) - 1)){
    for(c2 in (c1+1):length(cellList)){
        cT1 <- cellList[c1]
        cT2 <- cellList[c2]
        print(cT1)
        print(cT2)
        
        bulkPath1 <- paste0("/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Bulk/ContactMatrix/", cT1)
        bulkPath2 <- paste0("/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Bulk/ContactMatrix/", cT2)
        methodPath1 <- paste0("/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Duan2020/ContactMatrix/", cT1)
        methodPath2 <- paste0("/p/keles/schic/volumeC/YeZheng_workingDirectory/results/Duan2020/ContactMatrix/", cT2)

    for(chr in chromList){
        print(chr)
        bulkData1 <- paste0(bulkPath1, "/rep1/", chr, ".nij.matrix")
        bulkData2 <- paste0(bulkPath1, "/rep2/", chr, ".nij.matrix")

        bulkData3 <- paste0(bulkPath2, "/rep1/", chr, ".nij.matrix")
        bulkData4 <- paste0(bulkPath2, "/rep2/", chr, ".nij.matrix")
        
        resBulk <- runDiffHic(repList, bulkData1, bulkData2, bulkData3, bulkData4, resolution)
        resBulkTS <- resBulk %>% filter(p.adj <= 0.1)
        resBulkSign <- resBulk %>% filter(p.adj <= 0.05)
       
        for(method in methodList){
            print(method)
            mData1 <- paste0(methodPath1, "/", method, "_rep1/", chr, ".nij.matrix")
            mData2 <- paste0(methodPath1, "/", method, "_rep2/", chr, ".nij.matrix")

            mData3 <- paste0(methodPath2, "/", method, "_rep1/", chr, ".nij.matrix")
            mData4 <- paste0(methodPath2, "/", method, "_rep2/", chr, ".nij.matrix")

            resMethod <- runDiffHic(repList, mData1, mData2, mData3, mData4, resolution)
            
            summaryD <- data.frame(signN = c(res %>% filter(p.adj <= 0.001) %>% nrow(), res %>% filter(p.adj <= 0.01) %>% nrow(), res %>% filter(p.adj <= 0.05) %>% nrow(), res %>% filter(p.adj <= 0.1) %>% nrow()), padj = c(0.001, 0.01, 0.05, 0.1), method1 = "Bulk", method2 = method, chrom = chr, cellType = cT)

        }
        
    }
       
    
}




            
