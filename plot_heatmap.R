library(pheatmap)

mapping <- read.table("./DR_OM_SE_CP_AS_MA_PA_MappingTable.csv", sep="\t", header=T, row.names=1)
mapping[1:6, 1:6]
mapping[row(mapping) == col(mapping)] <- 1
mapping[1:6, 1:6]

pdf("./heatmap.pdf", width=50, height=50)
pheatmap(mapping, color=colorRampPalette(c("lightblue2", "white", "#FC7956"))(100), fontsize=30)
dev.off()
