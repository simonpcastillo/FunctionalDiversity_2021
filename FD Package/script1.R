  #0. Preamble
  # Load required package
  
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, ggrepel, gridExtra, data.table, geometry, plyr, picante, stringr, RColorBrewer,
               qpcR, rowr, reshape2,SYNCSA)

  
  #Load required fucntions
  source("hullarea.R")
  source("FDind.R")
  
  #Create directories to save PC plots
  dir.create("PC_plots")
  dir.create("FD_plots")
  dir.create("Reports")
  
  
  #Load data
  data0 <- read.csv("EcoGenera.csv", header=T, row.names = 1, sep=",", check.names = F)
  #colnames(data0) <- paste0("time", 0:(ncol(data0)-1))  
  timeSums <- colSums(data0)
  for (i in 1:nrow(data0)) {
    for(j in 1:ncol(data0)){
      data0[i,j] <- data0[i,j]/timeSums[j]
    }
  }
  data1<-melt(data0)
  data1<-cbind(rep(rownames(data0, ncol(data0))),data1)
  data1<- cbind(data1$`rep(rownames(data0, ncol(data0)))`,
                as.data.frame(str_split_fixed(data1$`rep(rownames(data0, ncol(data0)))`, n=8, pattern = "")),
                data1$variable,
                data1$value)
  names(data1) <- c("ecocode",letters[1:8],"time", "abundance")
  #Setting up traits data & mergind data for analyses
  traits<-as.data.frame(matrix(0, ncol = 48, nrow = nrow(data1)))
  colnames(traits)<-paste0(rep(letters[1:8], each=6),1:6)
  data1<-cbind(No=array(NA, nrow(data1)),data1,traits)
  
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,3])
    data1[j,(b+12)]<-1
  }
  
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,4])
    data1[j,(b+18)]<-1
  }
  
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,5])
    data1[j,(b+24)]<-1
  }  
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,6])
    data1[j,(b+30)]<-1
  }
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,7])
    data1[j,(b+36)]<-1
  }
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,8])
    data1[j,(b+42)]<-1
  }
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,9])
    data1[j,(b+48)]<-1
  }
  for (j in 1:nrow(data1)) {
    b<-as.numeric(data1[j,10])
    data1[j,(b+54)]<-1
  }
  
  #1. Analyses
  times<-as.character(unique(data1$time))
  summary1<-data.frame()
  for (i in 1:length(unique(data1$time))) {
    stime<-times[i]
    print(stime)
    
    dataAb<-t(data1$abundance[which(data1$time == stime)])
    colnames(dataAb)<- data1$ecocode[which(data1$time == stime)]
    
    dataTraits<-data1[which(data1$time == stime),14:ncol(data1)]
    rownames(dataTraits) <- as.character(data1[which(data1$time == stime),]$ecocode)
    rd<-rao.diversity(dataAb,traits = dataTraits)
    summary1<-rbind(summary1, cbind(Time=stime, Simpson=rd$Simpson, FunRao=rd$FunRao, FunRedundancy=rd$FunRedundancy))
  }
  
  ## Save report summary 1: Functional metrics (Simpson, Rao's quadratic Entropy and Funcional redundancy)
  write.csv(summary1, "Reports/Summary1_RaoFD.csv")
  print("<Summary1_RaoFD.csv> has been done!")
  
  # Summary 2
  data2 <- data1
  times<-as.character(unique(data2$time))
  summary2<-data.frame()
  
  for (i in 1:length(times)) {
    weight <- as.data.frame(t(subset(data2, time== times[i], select = "abundance")))
    coord  <- subset(data2, time== times[i], select = 3:5)
    nnn    <- subset(data2, time== times[i], select = "ecocode")
    colnames(weight)<-rownames(coord)<-nnn$ecocode
    rownames(weight)<- times[i]
    names <- names(weight[which(weight == 0)])
    
    coord<-coord[-which(rownames(coord) %in% names),]
    weight <- weight[-which(colnames(weight) %in% names)]
    for (m in 1:ncol(coord)) {
      coord[,m] = as.numeric(coord[,m])
      
    }
    
    
    GT_pca<- prcomp(coord, retx=TRUE, center=TRUE, scale =TRUE)
    pc1<-round(as.data.frame(summary(GT_pca)[6])[2,1],3)
    pc2<-round(as.data.frame(summary(GT_pca)[6])[2,2],3)
    pc3<-round(as.data.frame(summary(GT_pca)[6])[2,3],3)
    ch1<-paste0("pc1_", times[i])
    ch2<-paste0("pc2_", times[i])
    ch3<-paste0("pc3_", times[i])
    assign(ch1, pc1)
    assign(ch2,pc2)
    assign(ch3, pc3)
    
    tGT_pca <- GT_pca$x
    GT_dat <- data.frame(tGT_pca)
    GPC1<-GT_dat$PC1
    GPC2<-GT_dat$PC2
    GPC3<-GT_dat$PC3
    df0 <-rbind(GPC1,GPC2)
    df1=rbind(GPC1,GPC3) 
    df10<-t(df0)
    df11=t(df1)
    con.hull.pos1 <- chull(df10)
    con.hull.pos2 <- chull(df11) 
    con.hull1 <- rbind(df10[con.hull.pos1,],df10[con.hull.pos1[1],])
    con.hull2 <- rbind(df11[con.hull.pos2,],df11[con.hull.pos2[1],])
    tr <- coord
    ab <- as.matrix(weight)
    FRic<- FDind(traits = tr, abundances = ab)$FRic
    FSpe<- FDind(traits = tr, abundances = ab)$FSpe
    FDiv<- FDind(traits = tr, abundances = ab)$FDiv
    FEve<- FDind(traits = tr, abundances = ab)$FEve
    chull_PC12 <- convhulln(cbind(GPC1, GPC2), "FA")$area#hullarea(GPC1,GPC3) 
    chull_PC13 <- convhulln(cbind(GPC1, GPC3), "FA")$area#hullarea(GPC1,GPC3) 
    PCS<-data.frame(GT_dat)
    ab<-t(weight)
    ab<-data.frame(ab)
    PCS<-cbind.fill(PCS,ab)
    PCS<-data.frame(PCS)
    rownames(PCS)<-rownames(ab)
    con.hull1<-data.frame(con.hull1)
    con.hull2<-data.frame(con.hull2)
    datamultidim1<-c(time=times[i],FEve=FEve, FDiv=FDiv, FSpe=FSpe,FRic = FRic, chull_PC12= chull_PC12, chull_PC13=chull_PC13)
    summary2.1<-as.data.frame(t(datamultidim1))
    colnames(summary2.1)<-c( "time","FEve","FDiv","FSpe","FRic","chull_PC12", "chull_PC13")
    summary2 <- rbind(summary2, summary2.1)
    npcs<-paste0("PCS_", times[i])
    assign(npcs, PCS)
    ch1<-paste0("con.hull1_", times[i])
    ch2<-paste0("con.hull2_", times[i])
    assign(ch1, con.hull1)
    assign(ch2, con.hull2)
    print(times[i])
  }
  
  
  # Save a report named "Summary2_FDmetrics.csv" with FEve, FSpe, and FDis
  write.csv(summary2, "Reports/Summary2_FDmetrics.csv")
  print("<Summary2_FDmetrics.csv> has been created!")
  
  # This chunk prepares summary data normalisation (range 0-1) for plots
      norm_summary1 <- summary1
      for (j in 2:ncol(norm_summary1)) {
        norm_summary1[,j] <- as.numeric(as.character(norm_summary1[,j]))
      }
      norm_summary2 <- summary2
      for (j in 2:ncol(norm_summary2)) {
        norm_summary2[,j] <- as.numeric(as.character(norm_summary2[,j]))
      }
    
      norm_summary1 <- melt(norm_summary1, id=c("Time"))
      
      norm_summary1$norm_value <-NA
      for(i in 1:nrow(norm_summary1)){
        norm_summary1$norm_value[i] <- (norm_summary1$value[i]-min(norm_summary1[which(norm_summary1$variable == norm_summary1$variable[i]),]$value))/(max(norm_summary1[which(norm_summary1$variable == norm_summary1$variable[i]),]$value)-min(norm_summary1[which(norm_summary1$variable == norm_summary1$variable[i]),]$value))
        }
      norm_summary2 <- melt(norm_summary2, id=c("time"))
      norm_summary2$norm_value <-NA
      for(i in 1:nrow(norm_summary2)){
        norm_summary2$norm_value[i] <- (norm_summary2$value[i]-min(norm_summary2[which(norm_summary2$variable == norm_summary2$variable[i]),]$value))/(max(norm_summary2[which(norm_summary2$variable == norm_summary2$variable[i]),]$value)-min(norm_summary2[which(norm_summary2$variable == norm_summary2$variable[i]),]$value))
      }

#3. Plots
    #FD plots
    a <- ggplot(data = norm_summary1, aes(x=Time, y = value, colour = variable, group= variable)) +
      geom_point() +
      geom_line()+
      labs(x ="Time", y = "Functional diversity measure", color="Measure") +
      theme_bw()
    
    b <- ggplot(data = norm_summary1, aes(x=Time, y = norm_value, colour = variable, group= variable)) +
      geom_point() +
      geom_line()+
      labs(x ="Time", y = "Functional diversity normalized measure", color="Measure") +
      theme_bw()
    
    c <- ggplot(data = norm_summary2[which(norm_summary2$variable %in% c("FEve", "FDiv", "FSpe", "FRic")),], 
           aes(x=time, y = value, colour = variable, group= variable)) +
      geom_point() +
      geom_line()+
      labs(x ="Time", y = "Functional diversity measure", color="Measure") +
      theme_bw()
      
    d <- ggplot(data = norm_summary2[which(norm_summary2$variable %in% c("FEve", "FDiv", "FSpe", "FRic")),], 
           aes(x=time, y = norm_value, colour = variable, group= variable)) +
      geom_point() +
      geom_line()+
      labs(x ="Time", y = "Functional diversity normalized measure", color="Measure") +
      theme_bw()
    
    
    e <- ggplot(data = norm_summary2[which(norm_summary2$variable %in% c("FRic")),], 
                aes(x=time, y = norm_value, colour = variable, group= variable)) +
      geom_point() +
      geom_line()+
      labs(x ="Time", y = "Functional Richness (normalized measure)") +
      theme_bw()
    
    f <- ggplot(data = norm_summary2[which(norm_summary2$variable %in% c("FRic")),], 
                aes(x=time, y = value, colour = variable, group= variable)) +
      geom_point() +
      geom_line()+
      labs(x ="Time", y = "Functional Richness") +
      theme_bw()
    
    tiff(paste0("FD_plots/1.1",".tiff"), res = 300, height = 1200, width = 2400)
    print(a)
    dev.off()  
   
    tiff(paste0("FD_plots/plot1.2", ".tiff"), res = 300, height = 1200, width = 2400)
    print(b)
    dev.off()  
   
    tiff(paste0("FD_plots/2.1",".tiff"), res = 300, height = 1200, width = 2400)
    print(c)
    dev.off()  
   
    tiff(paste0("FD_plots/plot2.2", ".tiff"), res = 300, height = 1200, width = 2400)
    print(d)
    dev.off() 
    
    tiff(paste0("FD_plots/plot3.1_FRic", ".tiff"), res = 300, height = 1200, width = 2400)
    print(e)
    dev.off() 
    
    tiff(paste0("FD_plots/plot3.2_FRic", ".tiff"), res = 300, height = 1200, width = 2400)
    print(f)
    dev.off()

   
    #PC plots
    
    for (i in 1:length(times)) {
      plotit<-data.frame(id=data1$ecocode[which(data1$time == times[i])], abplot= data1$abundance[which(data1$time == times[i])])
      ab<-plotit$abplot[-which(plotit$abplot == 0)]
      colPoints <-rainbow(n = length(ab))
      PCSname <- paste0("PCS_",times[i])
      dataPCS<- get(PCSname)
      dataconhull1 <-  get(paste0("con.hull1_",times[i]))
      dataconhull2 <-  get(paste0("con.hull2_",times[i]))
      pc1<-get(paste0("pc1_",times[i]))
      pc2<-get(paste0("pc2_",times[i]))
      pc3<-get(paste0("pc3_",times[i]))
      
      p<-ggplot(data=dataPCS, aes(x=PC1, y= PC2 ))  + 
        geom_polygon(data = dataconhull1, aes( x= GPC1, y=GPC2), alpha = 0.6, fill=rgb(.7,.8,.9,0.2)) +  #convex hull
        geom_point(shape=21, size=ab*100, fill= colPoints, colour="black", alpha=.75) + #points
        geom_point(data=dataPCS, aes(x=mean(PC1), y=mean(PC2)), shape=3, size=4, colour = rgb(1, 0, 0)) + #Mean point (red + at the middle)
       # geom_text_repel(aes(label = rownames(dataPCS)), size = 4, col = colPoints, nudge_y = -0.3, nudge_x = -.1, segment.colour = "black") +  #labels
        scale_colour_manual("black") +
        labs( x=paste("PC1 ",pc1*100, "%", sep =""), y=paste("PC2 ",pc2*100, "%", sep =""))  +  #labs name
        annotate("text", x=-2, y=2.8, label=paste0("FRic=", round(as.numeric(as.character(summary2$FRic[i])), 2)), size = 4, family="Times",fontface="bold", colour="black") +
        annotate("text", x=-2, y=2.5, label=paste0("FRed=", round(as.numeric(as.character(summary1$FunRedundancy[i])),2)), size = 4, family="Times",fontface="bold", colour="black") +
        xlim(-3,3)+
        ylim(-3,3)+
        theme_bw() +  
        theme(axis.title.x=element_text(size=15), axis.title.y = element_text(size=15),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      q <-ggplot(data=dataPCS, aes(x=PC1, y= PC3 ))  + 
        geom_polygon(data = dataconhull2, aes( x= GPC1, y=GPC3), alpha = 0.6, fill=rgb(.7,.8,.9,0.2)) +  #convex hull
        geom_point(shape=21, size=ab*100, fill= colPoints, colour="black", alpha=.75) + #points
        geom_point(data=dataPCS, aes(x=mean(PC1), y=mean(PC3)), shape=3, size=4, colour = rgb(1, 0, 0)) + #Mean point (red + at the middle)
        #geom_text_repel(aes(label = rownames(dataPCS)), size = 4, col = colPoints, nudge_y = -0.3, nudge_x = -.1, segment.colour = "black") +  #labels
        scale_colour_manual("black") +
        labs( x=paste("PC1 ",pc1*100, "%", sep =""), y=paste("PC3 ",pc3*100, "%", sep =""))  +  #labs name
        annotate("text", x=-2, y=2.8, label=paste0("FRic=", round(as.numeric(as.character(summary2$FRic[i])), 2)), size = 4, family="Times",fontface="bold", colour="black") +
        annotate("text", x=-2, y=2.5, label=paste0("FRed=", round(as.numeric(as.character(summary1$FunRedundancy[i])),2)), size = 4, family="Times",fontface="bold", colour="black") +
        xlim(-3,3)+
        ylim(-3,3)+
        theme_bw() +  
        theme(axis.title.x=element_text(size=15), axis.title.y = element_text(size=15),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      
      tiff(paste0("PC_plots/Plot_PC12-PC13_", times[i], ".tiff"), res = 300, height = 1200, width = 2400)
      print(grid.arrange(p, q, ncol=2, nrow=1))
      dev.off()  
      
    }   
    
    
    ### Overlayed PC plots
    for (i in 1:(length(times)-1)) {
      plotit<-data.frame(id=data1$ecocode[which(data1$time == times[i])], abplot= data1$abundance[which(data1$time == times[i])])
      ab<-plotit$abplot[-which(plotit$abplot == 0)]
      colPoints <-rainbow(n = length(ab))
      PCSname <- paste0("PCS_",times[i])
      dataPCS<- get(PCSname)
      PCSname2 <- paste0("PCS_",times[i+1])
      dataPCS2<- get(PCSname2)
      dataconhull1 <-  get(paste0("con.hull1_",times[i]))
      dataconhull1.2 <-  get(paste0("con.hull1_",times[i+1]))
      dataconhull2 <-  get(paste0("con.hull2_",times[i]))
      dataconhull2.2 <-  get(paste0("con.hull2_",times[i+1]))
      pc1<-get(paste0("pc1_",times[i]))
      pc2<-get(paste0("pc2_",times[i]))
      pc3<-get(paste0("pc3_",times[i]))
      pc1.2<-get(paste0("pc1_",times[i+1]))
      pc2.2<-get(paste0("pc2_",times[i+1]))
      pc3.2<-get(paste0("pc3_",times[i+1]))
      
      p<- ggplot(data=dataPCS, aes(x=PC1, y= PC2 ))  + 
        geom_polygon(data = dataconhull1, aes( x= GPC1, y=GPC2), alpha = 0.6, fill=rgb(0,0,0)) +  #convex hull
        geom_polygon(data = dataconhull1.2, aes( x= GPC1, y=GPC2), alpha = 0.6, fill=rgb(1,0,0)) +  #convex hull
        #geom_point(shape=21, size=ab*100, fill= colPoints, colour="black", alpha=.75) + #points
        #geom_point(data=dataPCS, aes(x=mean(PC1), y=mean(PC2)), shape=3, size=4, colour = rgb(1, 0, 0)) + #Mean point (red + at the middle)
        #geom_text_repel(aes(label = rownames(dataPCS)), size = 4, col = colPoints, nudge_y = -0.3, nudge_x = -.1, segment.colour = "black") +  #labels
        scale_colour_manual("black") +
        labs( x="PC1", y="PC2")  +  #labs name
        #annotate("text", x=-2, y=2.8, label=i) + 
        # annotate("text", x=-2, y=2.8, label=paste0("FRic=", round(as.numeric(as.character(summary2$FRic[i])), 2)), size = 4, family="Times",fontface="bold", colour="black") +
        #annotate("text", x=-2, y=2.5, label=paste0("FRed=", round(as.numeric(as.character(summary1$FunRedundancy[i])),2)), size = 4, family="Times",fontface="bold", colour="black") +
        xlim(-3,3)+
        ylim(-3,3)+
        theme_bw() +  
        theme(axis.title.x=element_text(size=15), axis.title.y = element_text(size=15),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      
      q <-ggplot(data=dataPCS, aes(x=PC1, y= PC3 ))  + 
        geom_polygon(data = dataconhull2, aes( x= GPC1, y=GPC3), alpha = 0.6, fill=rgb(0,0,0)) +  #convex hull
        geom_polygon(data = dataconhull2.2, aes( x= GPC1, y=GPC3), alpha = 0.6, fill=rgb(1,0,0)) +  #convex hull
        #geom_point(shape=21, size=ab*100, fill= colPoints, colour="black", alpha=.75) + #points
        #geom_point(data=dataPCS, aes(x=mean(PC1), y=mean(PC3)), shape=3, size=4, colour = rgb(1, 0, 0)) + #Mean point (red + at the middle)
        #geom_text_repel(aes(label = rownames(dataPCS)), size = 4, col = colPoints, nudge_y = -0.3, nudge_x = -.1, segment.colour = "black") +  #labels
        scale_colour_manual("black") +
        labs( x="PC1", y="PC2")  +  #labs name
        #labs( x=paste("PC1 ",pc1*100, "%", sep =""), y=paste("PC3 ",pc3*100, "%", sep =""))  +  #labs name
        #annotate("text", x=-2, y=2.8, label=paste0("FRic=", round(as.numeric(as.character(summary2$FRic[i])), 2)), size = 4, family="Times",fontface="bold", colour="black") +
        #annotate("text", x=-2, y=2.5, label=paste0("FRed=", round(as.numeric(as.character(summary1$FunRedundancy[i])),2)), size = 4, family="Times",fontface="bold", colour="black") +
        xlim(-3,3)+
        ylim(-3,3)+
        theme_bw() +  
        theme(axis.title.x=element_text(size=15), axis.title.y = element_text(size=15),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      
      tiff(paste0("PC_plots/Plot_PC12-PC13_Over", times[i],"-",times[i+1], ".tiff"), res = 300, height = 2400, width = 1200)
      print(grid.arrange(p, q, ncol=1, nrow=2))
      dev.off()  
      
    } 
  
    