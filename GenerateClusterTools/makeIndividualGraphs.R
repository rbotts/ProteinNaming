library(ggplot2)
library(Rmisc)
ClusterGroups <- read.csv(file="ClusterGroups.csv",head=FALSE,sep=",")
ClusterGroups <- na.omit(ClusterGroups)
ClusterGroups<-ClusterGroups[which(ClusterGroups$V1!="Fip"),]
ClusterGroups<-ClusterGroups[which(ClusterGroups$V1!="mpR"),]

x=c()
for (id in unique(ClusterGroups$V1))
{
  x<-c(x,(nrow(ClusterGroups[which(ClusterGroups$V1==id),])))
}
maxx <- max(x)

#==========================================================
# I don't know why, but this for loop had to be moved outside
# of the function. R had problems with the cids. It might have
# something to do with multithreading. -ZL (5/18/2017)
#==========================================================
ClusterGroups$V8 <- FALSE
ClusterGroups$V9 <- FALSE
references <- tolower(c("R7K","R46","pRA3","pKJK5","RK2","RP4","pNDM-1_Dok01"))
cids <- c()
for (cRow in c(1:nrow(ClusterGroups)))
{
  plasmid <- strsplit(as.character(ClusterGroups$V3[cRow]),"_")[[1]]
  if(length(plasmid)>2)
  {
    plasmid <- plasmid[3:length(plasmid)]
    plasmid <- paste(plasmid,collapse="_")
  }
  else
  {
    plasmid <- plasmid[2]
  }
  if (tolower(plasmid) %in% references)
  {
    cids <- c(cids, ClusterGroups$V2[cRow])
    ClusterGroups$V9[cRow] <- TRUE
  }
}

ClusterGroups[which(ClusterGroups$V2 %in% cids),"V8"] <- TRUE


makePlot <- function (backbone)
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
  CG3_AA <- CG2 # BY AA CLUSTER INSTEAD!
  CG3$V7<-c(1:nrow(CG3))
  CG3_AA$V7<-c(1:nrow(CG3_AA))
  
  ggplot(CG3,aes(x=c(1:nrow(CG3)),y=CG3$V4))+ # Notice X position is determined by V7
    geom_bar(stat="identity",
             fill=CG3$V6,
             alpha=CG3$V5
    )+
    scale_y_continuous()+
    scale_x_discrete()+ # Remove extra space
    geom_label(size=1.00,y=max(CG3$V4)/2.0,
              fontface=ifelse(CG3$V8,"bold","plain"),
              label=CG3$V3,
              color=ifelse(CG3$V9,"red","black"),
              label.padding = unit(0.05,"lines"),
              label.r = unit(0.05,"lines"),
              label.size = .01,
              alpha = .5) +
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
         width=3.0*25, 
         units="mm", 
         limitsize = FALSE)

  ggplot(CG3_AA,aes(x=c(1:nrow(CG3_AA)),y=CG3_AA$V4))+ # Notice X position is determined by V7
    geom_bar(stat="identity",
             fill=CG3_AA$V6,
             alpha=CG3_AA$V5
    )+
    scale_y_continuous()+
    scale_x_discrete()+ # Remove extra space
    geom_label(size=1.00,y=max(CG3_AA$V4)/2.0,
              fontface=ifelse(CG3_AA$V8,"bold","plain"),
              label=CG3_AA$V3,
              color=ifelse(CG3$V9,"red","black"),
              label.padding = unit(0.05,"lines"),
              label.r = unit(0.05,"lines"),
              label.size = .01,
              alpha = .5) +
    geom_text(size=1.5,y=-5,
              aes(label=rev(CG3_AA$V7)),
              color="black") +
    coord_flip() +
    ggtitle(backbone)+
    labs(y="Sequence Length",x="")+
    theme(
      axis.text.y=element_text(size=1),
      axis.title.y=element_text(size=2))
  ggsave(filename = paste("BarGraphs/",backbone,"_AA.png",sep=''), 
         height = 25+as.integer(25*(11*nrow(CG3))/150.0), 
         width=3.0*25, 
         units="mm", 
         limitsize = FALSE)
  }
 
sapply(unique(ClusterGroups$V1),makePlot)
# 
# multiplot(plotlist=myPlots[1:length(myPlots)],cols=length(myPlots))
# 
# dev.copy2pdf(file="PLOT.pdf",width=200,height=200,out.type = "pdf")
