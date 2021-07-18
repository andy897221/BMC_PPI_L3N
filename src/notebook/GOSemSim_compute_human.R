# install libraries instruction
# if(!requireNamespace('BiocManager', quietly=TRUE))
#   install.packages('BiocManager')
# BiocManager::install('GOSemSim')
# BiocManager::install('GO.db')
# BiocManager::install('org.Hs.eg.db')


library(RJSONIO)
library("rjson")
library(GOSemSim)
library('org.Hs.eg.db')
hs <- org.Hs.eg.db
hsGO <- godata('org.Hs.eg.db', ont='MF', keytype='SYMBOL')


dir = "./linkPred_out_human_reduced" # ./linkPred_out_reduced but only human datasets
# dir = "./tmp"
for (jsonF in list.files(path=dir)) {
    # parse all ppi that needs to compute GOSemSim
    if (grepl('topScore', jsonF, fixed=TRUE) == TRUE) next

    allPPIs <- list()
    prefix <- strsplit(jsonF, "_topPPI")[[1]][1]
    data <- fromJSON(file = paste(dir,"/",jsonF, sep=""))
    for (i in 1:length(data)) {
        for (j in 1:length(data[[i]]))
            allPPIs[[length(allPPIs)+1]] <- data[[i]][[j]]
    }

    allPPIs <- unique(allPPIs)
    geneScores <- list()
    curLen <- 0
    for (i in 1:length(allPPIs)) {
        ppi <- allPPIs[[i]]
        res <- geneSim(ppi[[1]], ppi[[2]], semData=hsGO, measure='Wang', combine='BMA')
        if (is.na(res)) {
            geneScores[i] <- list(NULL) # allow NA score
        } else {
            geneScores[[i]] <- res$geneSim
        }
    }

    writeF <- toJSON(allPPIs)
    write(writeF, paste("./GOSemSim_out/", prefix, "_PPI.json"))

    writeF <- toJSON(geneScores)
    write(writeF, paste("./GOSemSim_out/", prefix, "_GOSemSim.json"))
}