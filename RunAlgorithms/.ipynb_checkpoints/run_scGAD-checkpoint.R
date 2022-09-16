library(BandNorm)

inPath <- "/data/Lee2019/Counts_100kb/"
outPath <- "/results/Lee2019_scGAD/RDS/"

gad_score = scGAD(path = inPath, genes = hg19Annotations, depthNorm = TRUE, cores = 3, threads = 8)
head(gad_score)

saveRDS(gad_score, file = paste0(outPath, "gad_score_100kb.rds")) 
