library(ggplot2)
library(Rmisc)
ClusterGroups <- read.csv(file="ClusterGroups_8-23.csv",head=FALSE,sep=",")
ClusterGroups <- na.omit(ClusterGroups)
ClusterGroups<-ClusterGroups[which(ClusterGroups$V1!="Fip"),]
ClusterGroups<-ClusterGroups[which(ClusterGroups$V1!="mpR"),]

x=c()
for (id in unique(ClusterGroups$V1))
{
  x<-c(x,(nrow(ClusterGroups[which(ClusterGroups$V1==id),])))
}
maxx <- max(x)

for (backbone in unique(ClusterGroups$V1))
{
  CG2 <- ClusterGroups[which(ClusterGroups$V1==backbone),]
  cols <- rainbow(length(unique(CG2$V2)))
  CG2$V6 <- c(1:nrow(CG2))
  
  index = 1
  for (i in unique(CG2$V2))
  {
    CG2$V6[which(CG2$V2==i)] <- cols[index]
    index = index + 1
  }
  # rC <- nrow(CG2)
  # if (nrow(CG2) != maxx)
  # {
  #   hax<-data.frame(matrix("empty", nrow=(maxx-nrow(CG2)), ncol=6))
  #   names(hax) <- names(CG2)
  #   hax$V1 <- "empty"
  #   hax$V2 <- "empty"
  #   hax$V3 <- as.character(c(1:nrow(hax)))
  #   hax$V4 <- 1
  #   hax$V6 <- "#FFFFFF"
  #   CG2<-rbind(CG2,hax)
  # }
  CG2$V5[which(CG2$V5=="*")] <- "100.0"
  CG2$V5<-as.double(as.character(CG2$V5))
  CG2$V5[which(is.na(CG2$V5))] <- 0.0
  CG2$V4 <- as.double(CG2$V4)
  CG2$V5 <- CG2$V5/100.00
  #remove(hax)
  
  CG3 <- CG2[order(CG2$V4),] # CRITICAL: Allow the reordering of colors!
  #CG3 <- CG2 # BY AA CLUSTER INSTEAD!
  CG3$V7<-c(1:nrow(CG3))
  
  references <- tolower(c("R7K","R46","pRA3","pKJK5","RK2","RP4","pNDM-1_Dok01"))
  cids <- c()
  for (cRow in c(1:nrow(CG3)))
  {
    plasmid <- strsplit(as.character(CG3$V3[cRow]),"_")[[1]]
    if(length(plasmid)>2)
    {
      plasmid <- plasmid[2:length(plasmid)]
      plasmid <- paste(plasmid,collapse="_")
    }
    else
    {
      plasmid <- plasmid[2]
    }
    if (plasmid %in% references)
    {
      cids <- c(cids, CG3$V2[cRow])
    }
  }
  cids2 <- unique(cids)
  CG3$V8 <- FALSE
  CG3[which(CG3$V2 %in% cids),"V8"] <- TRUE
  
  if (is.null(CG3))
    break
  
  ggplot(CG3,aes(x=c(1:nrow(CG3)),y=CG3$V4))+ # Notice X position is determined by V7
    geom_bar(stat="identity",
             fill=CG3$V6,
             alpha=CG3$V5
    )+
    scale_y_continuous()+
    scale_x_discrete()+ # Remove extra space
    geom_text(size=1.00,y=max(CG3$V4)/2.0,
              fontface=ifelse(CG3$V8,"bold","plain"),
              label=CG3$V3,
              color="black") +
    geom_text(size=1.5,y=-5,
               aes(label=rev(CG3$V7)),
               color="black") +
    coord_flip() +
    ggtitle(backbone)+
    labs(y="Sequence Length",x="")+
    theme(
      axis.text.y=element_text(size=1),
      axis.title.y=element_text(size=2))
  ggsave(filename = paste("BarGraphs/",backbone,".png",sep=''), 
         height = 25+as.integer(25*(11*nrow(CG3))/150.0), 
         width=2.2*25, 
         units="mm", 
         limitsize = FALSE)}
 
# 
# multiplot(plotlist=myPlots[1:length(myPlots)],cols=length(myPlots))
# 
# dev.copy2pdf(file="PLOT.pdf",width=200,height=200,out.type = "pdf")
