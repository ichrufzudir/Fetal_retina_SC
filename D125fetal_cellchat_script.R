library(CellChat)
library(patchwork)
library(ggalluvial)
options(stringsAsFactors = FALSE)

D125CPfetalScellchat <- data.subset
D125CPfetalScellchat$labels <- D125CPfetalScellchat@active.ident
DimPlot(D125CPfetalScellchat, label = TRUE)

cellchat <- createCellChat(object = D125CPfetalScellchat, group.by = "labels")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
showDatabaseCategory(CellChatDB)
cellchat <- projectData(cellchat, PPI.human)

cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = TRUE) #check these options

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht1
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht2
ht1+ht2

netAnalysis_contribution(cellchat, signaling = "MK")
netVisual_aggregate(cellchat, signaling = c('MK'), layout = "chord", show.legend = TRUE)
pairLR.NT <- extractEnrichedLR(cellchat, signaling = 'MK', geneLR.return = FALSE)

LR.show <- pairLR.NT[2,]
netVisual_individual(cellchat, signaling = 'MK', pairLR.use = LR.show,
                     layout = "chord", show.legend = TRUE) #generate the chord plot

netAnalysis_contribution(cellchat, signaling = "ACTIVIN")
netVisual_aggregate(cellchat, signaling = c('MK'), layout = "chord", show.legend = TRUE)
pairLR.NT <- extractEnrichedLR(cellchat, signaling = 'MK', geneLR.return = FALSE)
LR.show <- pairLR.NT[c(1, 3, 5),]
