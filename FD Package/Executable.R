##############################################################
#This code runs for an input matrix saved as csv called 'EcoGenera.csv'
# if you want to replace that file, you MUST change it in this script (line 32) and in the Script1.R (line 28)
#
# INPUT REQUIREMENTS: Matrix S x T; with S being species and T being communities or time bins.
# 
#
# This script returns metrics (as csv stored in folders) and plots related to functional diversity calculated from
# Villeger et al  2009, Ecological Applications (Function FDind.R)
# Rao's related measures (Simpson, FunRao and FunRedundancy) were obatined using SYNCSA::rao.diversity; see referneces therein
# 
# Note that this plot also runs simmulation based on a simple reordering of the initial matrix (see nullcom_FD.R)
# to define the number of null communities to run set the value of 'nullrepl' (line 27; this script)

# By:  Simon P. Castillo spcastil[at].uc.cl

#0. Preamble

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, ggrepel, gridExtra, data.table, geometry, plyr, picante, stringr, RColorBrewer,
               qpcR, rowr, reshape2,SYNCSA)
source("hullarea.R")
source("FDind.R")
    #Load required scripts
    source("script1.R") #Here we call all the required packages
    source("nullcom_FD.R") # this is the routine that creates null communities
    
    #Create directories to save PC plots
    dir.create("null_FD_plots")
    dir.create("null_Reports")
    
    #To define number of null replicates 6 empty dataframe
    nullrepl<-100
    summary1null<-array(dim = 4)
    summary2null<-array(dim = 7)

    #Load original data
    data0 <- read.csv("EcoGenera.csv", header=T, row.names = 1, sep=",", check.names = F)
    #colnames(data0) <- paste0("time", 0:(ncol(data0)-1)) 
    
  for (l in 1:nullrepl) {
    data_n <- nullcomm_FD(data0,type="all")
    timeSums <- colSums(data_n)
    for (i in 1:nrow(data_n)) {
      for(j in 1:ncol(data_n)){
        data_n[i,j] <- data_n[i,j]/timeSums[j]
      }
    }
    data1<-melt(data_n)
    data1<-cbind(rep(rownames(data_n, ncol(data_n))),data1)
    data1<- cbind(data1$`rep(rownames(data_n, ncol(data_n)))`,
                  as.data.frame(str_split_fixed(data1$`rep(rownames(data_n, ncol(data_n)))`, n=8, pattern = "")),
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
    summaryn1<-data.frame()
    for (i in 1:length(unique(data1$time))) {
      stime<-times[i]
      print(stime)
      
      dataAb<-t(data1$abundance[which(data1$time == stime)])
      colnames(dataAb)<- data1$ecocode[which(data1$time == stime)]
      
      dataTraits<-data1[which(data1$time == stime),14:ncol(data1)]
      rownames(dataTraits) <- as.character(data1[which(data1$time == stime),]$ecocode)
      rd<-rao.diversity(dataAb,traits = dataTraits)
      summaryn1<-rbind(summaryn1, cbind(Time=stime, Simpson=rd$Simpson, FunRao=rd$FunRao, FunRedundancy=rd$FunRedundancy))
    }
    
    ## Save report summary 1: Functional metrics (Simpson, Rao's quadratic Entropy and Funcional redundancy)
    print("<Summary1_RaoFD.csv_rep> has been done!")
    
    # Summary 2
    data2 <- data1
    times<-as.character(unique(data2$time))
    summaryn2<-data.frame()
    
    for (i in 1:length(times)) {
      
      weight <- as.data.frame(t(subset(data2, time== times[i], select = "abundance")))
      coord  <- subset(data2, time== times[i], select = 3:5)
      #biovol <- as.data.frame(t(subset(datafunc, time== times[i], select = "av.biovol")))
      nnn    <- subset(data2, time== times[i], select = "ecocode")
      colnames(weight)<-rownames(coord)<-nnn$ecocode
      rownames(weight)<- times[i]
      names <- names(weight[which(weight == 0)])
      
      coord<-coord[-which(rownames(coord) %in% names),]
      weight <- weight[-which(colnames(weight) %in% names)]
      for (m in 1:ncol(coord)) {
        coord[,m] = as.numeric(coord[,m])
        
      }
      
      
      GT_pca<- prcomp(coord, retx=TRUE, center=TRUE, scale =FALSE)
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
      summaryn2 <- rbind(summaryn2, summary2.1)
      npcs<-paste0("PCS_", times[i])
      assign(npcs, PCS)
      ch1<-paste0("con.hull1_", times[i])
      ch2<-paste0("con.hull2_", times[i])
      assign(ch1, con.hull1)
      assign(ch2, con.hull2)
    }
    # Save a report named "Summary2_FDmetrics.csv" with FEve, FSpe, FOri and FDis
    print("<Summary2_FDmetrics.csv_rep_> has been created!")
    
    summary1$nullcomm<-rep(l, times= nrow(summaryn1))
    summary2$nullcomm <- rep(l, times=nrow(summaryn2))
    
    summary1null<-rbind(summary1null,summaryn1)
    summary2null<-rbind(summary2null,summaryn2)
    print(l)
  }
    summary1null<-summary1null[-c(1),]
    summary2null<-summary2null[-c(1),]
    
    for (i in 2:(ncol(summary1null)-1)) {
      summary1null[,i] <- as.numeric(levels(summary1null[,i]))[summary1null[,i]]
    }
    for (i in 2:(ncol(summary2null)-1)) {
      summary2null[,i] <- as.numeric(levels(summary2null[,i]))[summary2null[,i]]
    }
    for (j in 2:ncol(summary1null)) {
      summary1null[,j] <- as.numeric(as.character(summary1null[,j]))
    }
    for (j in 2:ncol(summary2null)) {
      summary2null[,j] <- as.numeric(as.character(summary2null[,j]))
    }
    
    
    for (n in 1:length(times)) {
      summary1.2<-ddply(subset(summary1null, Time == times[n]), c("Time"), summarise, 
                        mSimpson = mean(Simpson, na.rm = T),
                        mFunRao = mean(FunRao, na.rm = T),
                        mFunRedun = mean(FunRedundancy, na.rm = T),
                        vSimpson = var(Simpson, na.rm = T),
                        vFunRao = var(FunRao, na.rm = T),
                        vFunRedun = var(FunRedundancy, na.rm = T),
                        Simpson.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(Simpson)/ length(Simpson)), 
                        Simpson.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(Simpson)/ length(Simpson)),
                        FunRao.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(FunRao)/ length(FunRao)), 
                        FunRao.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(FunRao)/ length(FunRao)),
                        FunRed.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(FunRedundancy)/ length(FunRedundancy)), 
                        FunRed.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(FunRedundancy)/ length(FunRedundancy))
                        )
      nm <- paste0("summary1_", times[n])
      assign(nm, summary1.2)
    }
    
    summary1.2null<-rbind(lapply(ls(pattern = "summary1_"), get)) 
    summary1.2null <- do.call("rbind", summary1.2null)
    
   
    
    
    for (n in 1:length(times)) {
      summary2.2<-ddply(subset(summary2null, time == times[n]), c("time"), summarise, 
                        mFEve = mean(FEve, na.rm = T),
                        mFDiv = mean(FDiv, na.rm = T),
                        mFSpe = mean(FSpe, na.rm = T),
                        mFRic = mean(FRic, na.rm = T),
                        mchull_PC12 = mean(chull_PC12, na.rm = T),
                        mchull_PC13 = mean(chull_PC13, na.rm = T),
                        vFEve = var(FEve, na.rm = T),
                        vFDiv = var(FDiv, na.rm = T),
                        vFSpe = var(FSpe, na.rm = T),
                        vFRic = var(FRic, na.rm = T),
                        vchull_PC12 = var(chull_PC12, na.rm = T),
                        vchull_PC13 = var(chull_PC13, na.rm = T),
                        FEve.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(FEve)/ length(FEve)), 
                        FEve.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(FEve)/ length(FEve)),
                        FDiv.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(FDiv)/ length(FDiv)), 
                        FDiv.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(FDiv)/ length(FDiv)),
                        FSpe.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(FSpe)/ length(FSpe)), 
                        FSpe.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(FSpe)/ length(FSpe)),
                        FRic.ic95= qt(1-((1-0.95)/2), (nullrepl- 1))* sqrt(var(FRic)/ length(FRic)), 
                        FRic.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))* sqrt(var(FRic)/ length(FRic)),
                        chull_PC12.ic95= qt(1-((1-0.95)/2), (nullrepl - 1))*sqrt(var(chull_PC12)/ length(chull_PC12)), 
                        chull_PC12.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))*sqrt(var(chull_PC12)/ length(chull_PC12)),
                        chull_PC13.ic95= qt(1-((1-0.95)/2), (nullrepl - 1))*sqrt(var(chull_PC13)/ length(chull_PC13)), 
                        chull_PC13.ic5= qt(1-((1-0.05)/2), (nullrepl - 1))*sqrt(var(chull_PC13)/ length(chull_PC13)))
      nm <- paste0("summary2_", times[n])
      assign(nm, summary2.2)
    }
    
    
    summary2.2null<-rbind(lapply(ls(pattern = "summary2_"), get)) 
    summary2.2null <- do.call("rbind", summary2.2null)
    
  write.csv(summary1.2null, "null_Reports/null_Summary1_RaoFD.csv")
  write.csv(summary2.2null, "null_Reports/null_Summary2__FDmetrics.csv")  
  
  
  #Plots
  for (j in 2:ncol(summary1)) {
    summary1[,j] <- as.numeric(as.character(summary1[,j]))
  }
  for (j in 2:ncol(summary2)) {
    summary2[,j] <- as.numeric(as.character(summary2[,j]))
  }
  
  n1<-ggplot(summary1, aes(Time, Simpson, group=1))+
    geom_point()+
    geom_line() +
    geom_line(data=summary1.2null, aes(x=Time, y=(mSimpson), group=1), colour="red") +
      geom_line(data=summary1.2null, aes(x=Time, y=(mSimpson + Simpson.ic95), group=1), colour="red", lty=2) +
    geom_line(data=summary1.2null, aes(x=Time, y=(mSimpson - Simpson.ic95), group=1), colour="red", lty=2) +
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))
  
    
  n2<-ggplot(summary1, aes(Time, FunRao, group=1))+
    geom_point()+
    geom_line() +
    geom_line(data=summary1.2null, aes(x=Time, y=(mFunRao), group=1), colour="red") +
    geom_line(data=summary1.2null, aes(x=Time, y=(mFunRao + FunRao.ic95), group=1), colour="red", lty=2) +
    geom_line(data=summary1.2null, aes(x=Time, y=(mFunRao - FunRao.ic95), group=1), colour="red", lty=2) +
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))
  
  n3<-ggplot(summary1, aes(Time, FunRedundancy, group=1))+
    geom_point()+
    geom_line() +
    geom_line(data=summary1.2null, aes(x=Time, y=(mFunRedun), group=1), colour="red") +
    geom_line(data=summary1.2null, aes(x=Time, y=(mFunRedun + FunRed.ic95), group=1), colour="red", lty=2) +
    geom_line(data=summary1.2null, aes(x=Time, y=(mFunRedun - FunRed.ic95), group=1), colour="red", lty=2) +
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))  
  
  n4<-ggplot(summary2, aes(time, FRic, group=1))+
    geom_point()+
    geom_line() +
    geom_line(data=summary2.2null, aes(x=time, y=(mFRic), group=1), colour="red") +
    geom_line(data=summary2.2null, aes(x=time, y=(mFRic + FRic.ic95), group=1), colour="red", lty=2) +
    geom_line(data=summary2.2null, aes(x=time, y=(mFRic - FRic.ic95), group=1), colour="red", lty=2) +
    labs(y="Functional richness") +
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))
  
  
  
  tiff(paste0("null_FD_plots/plotSimpson", ".tiff"), res = 300, height = 1200, width = 2400)
  print(n1)
  dev.off()   
  
  tiff(paste0("null_FD_plots/plotFunRao", ".tiff"), res = 300, height = 1200, width = 2400)
  print(n2)
  dev.off()   
  
  tiff(paste0("null_FD_plots/plotFunRedundancy", ".tiff"), res = 300, height = 1200, width = 2400)
  print(n3)
  dev.off() 
  
  tiff(paste0("null_FD_plots/plotFRic12", ".tiff"), res = 300, height = 1200, width = 2400)
  print(n4)
  dev.off() 
  