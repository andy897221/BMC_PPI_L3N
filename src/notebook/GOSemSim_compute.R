# install libraries instruction
# if(!requireNamespace('BiocManager', quietly=TRUE))
#   install.packages('BiocManager')
# BiocManager::install('GOSemSim')
# BiocManager::install('GO.db')
# BiocManager::install('org.Sc.sgd.db')


library(RJSONIO)
library("rjson")
library(GOSemSim)
library('org.Sc.sgd.db')
hs <- org.Sc.sgd.db
hsGO <- godata('org.Sc.sgd.db', ont='MF', keytype='GENENAME')

dir = "./linkPred_out_reduced"
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
        # print(geneScores[[i]])
    }

    writeF <- toJSON(allPPIs)
    write(writeF, paste("./GOSemSim_out/", prefix, "_PPI.json"))

    writeF <- toJSON(geneScores)
    write(writeF, paste("./GOSemSim_out/", prefix, "_GOSemSim.json"))
}