
CG <<- read.csv(file = "ClusterGroups_6-21-2017_1.csv", head = FALSE, sep = ",")
CG <<- na.omit(CG)

CG[which(CG$V2 %in% CG[which(CG$V8==1),"V2"]),"V8"] <- ( 1 + CG[which(CG$V2 %in% CG[which(CG$V8==1),"V2"]),"V8"] )

for (element in unique(CG$V1))
{
  CG[which(CG$V1 == element),]
}