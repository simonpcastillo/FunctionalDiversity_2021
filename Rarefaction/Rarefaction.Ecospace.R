# DATA UPLOAD ###########################
# if you updating the file make sure to label the columns in the csv file as follows (all caps):
# TIME  COLL  CLIM  ENV  CODE (other columns don't matter)
ecoall<- read.csv("Ecodata.csv")		 
ecoall<- ecoall[-which(ecoall$TIME %in% c(587.5)),] # removes the oldest timebin(s)
gencodes<- as.factor(round(tapply(ecoall$CODE,ecoall$GENUS,mean)))
ecomajors<- c(115,261,312,315,324,325,331,351,361,531)	# identify important ecocategory
ecoall$CODE<- as.factor(ecoall$CODE)			# changes the ecocode to factor

# SET NUMBER ITERATIONS ###################
times<- 10000	# 100 takes about 3 minutes on my laptop

# Create subsets of data based on grouping variable of your choice ###########################
ecotemp<- subset(ecoall,ecoall$CLIM %in% c("temperate","Temperate")) 	# temperate only
ecotrop<- subset(ecoall,ecoall$CLIM %in% c("tropical","Tropical")) 	# tropical only
ecoclast<- subset(ecoall,ecoall$ENV %in% "siliciclastic") 			# clastic only
ecocarb<- subset(ecoall,ecoall$ENV %in% "carbonate") 				# carbonate only

eco<- ecoall	# choose dataset to use (set to "all data" right now)

# TIME COVERAGE ###########################
tcov<- table(eco$TIME)				 # time coverage by occurrences
mylength<- function(x) y<- length(table(x))
ocov<- tapply(eco$COLL,eco$TIME,mylength)		# time coverage by number of collections
ss<-min(ocov)-1						# subsampling level (minimum coverage - 1)

# plots for data coverage ########################
# for different data "ylim", "at" and "labels" for y-axis may need to be adjusted
# occurrences
plot(-as.numeric(names(tcov)),tcov,ylim=c(0,50000),ylab="",xlab="",axes=F)
 rect(-550,-20,10,50000,col="gray",border=NA)
 abline(v=c(-250,-65),col="white",lwd=2)
 points(-as.numeric(names(tcov)),tcov,type="o",col="black",pch=16,cex=0.9)
 axis(1,tck=-0.02,padj=-0.8, cex.axis=0.8,at=seq(-600,0,100),labels=c("600","500","400","300","200","100","0"))
 axis(2,las=1,tck=-0.02,hadj=0.8, cex.axis=0.8,at=seq(0,50000,10000),labels=seq(0,50000,10000))
 mtext("age [MA]", side=1,line=1.5)
 mtext("number of occurrences", side=2,line=3)

# collections
plot(-as.numeric(names(ocov)),ocov,ylim=c(0,5000),ylab="",xlab="",axes=F)
 rect(-550,-10,10,5000,col="gray",border=NA)
 abline(v=c(-250,-65),col="white",lwd=2)
 points(-as.numeric(names(ocov)),ocov,type="o",col="black",pch=16,cex=0.9)
 axis(1,tck=-0.02,padj=-0.8, cex.axis=0.8,at=seq(-600,0,100),labels=c("600","500","400","300","200","100","0"))
 axis(2,las=1,tck=-0.02,hadj=0.8, cex.axis=0.8,at=seq(0,5000,1000),labels=seq(0,5000,1000))
 mtext("age [MA]", side=1,line=1.5)
 mtext("number of collections", side=2,line=3)

# SUBSAMPLING BY COLLECTIONS ###########################
tbins<- as.numeric(names(tcov))		# find all time bins

outdiv<- NULL
oute<- NULL
outeg<- NULL
for (i in tbins)
 {
 bin<- eco[which(eco$TIME==i),]		# select one time bin
 cid<- as.numeric(names(table(bin$COLL)))	# find all collections present in this bin
 div<- NULL
 ediv<- NULL
 egen<- NULL
 eocc2<- NULL
    for (n in 1:times)
      {    
      one<- bin[which(bin$COLL %in% sample(cid,ss)),]
      gen<- table(droplevels(one$GENUS))
	eg<- table(gencodes[which(names(gencodes) %in% names(gen))])
      ecot1<- table(droplevels(as.factor(one$CODE)))
	eocc1<- rbind(table(one$CODE))
	eocc2<- rbind(eocc2,eocc1)
      div<- rbind(div,length(gen))
      ediv<- rbind(ediv,length(ecot1))
      egen<- rbind(egen,eg)
      }
     divm<- c(i,median(div),mean(div),median(ediv),mean(ediv))
     eocc3<- apply(eocc2,2,mean)
     egen2<- apply(egen,2,mean)
     outdiv<- rbind(outdiv,divm)
     oute<- rbind(oute,eocc3)
     outeg<- rbind(outeg,egen2)
 }

egen_all<- cbind(tbins,outeg)		# number of genera per ecocategory by time bin
rownames(egen_all)<- 1:nrow(egen_all)
eocc_all<- cbind(tbins,oute)		# number of occurrences per ecocategory by time bin
rownames(eocc_all)<- 1:nrow(eocc_all)
colnames(outdiv)<- c("TIME","MED_DIV","MEAN_DIV","MED_EDIV","MEAN_EDIV")
rownames(outdiv)<- 1:nrow(outdiv)
outdiv<- as.data.frame(outdiv)	# diversity and ecodiversity estimates by time bin

# plot diversity of genera #####################
# for different data "ylim", "at" and "labels" for y-axis may need to be adjusted
plot(-outdiv$TIME,outdiv$MED_DIV,ylim=c(0,600),ylab="",xlab="",axes=F)
 rect(-550,-10,10,600,col="gray",border=NA)
 abline(v=c(-250,-65),col="white",lwd=2)
 points(-outdiv$TIME,outdiv$MED_DIV,type="o",col="black",pch=16,cex=0.9)
 axis(1,tck=-0.02,padj=-0.8, cex.axis=0.8,at=seq(-600,0,100),labels=c("600","500","400","300","200","100","0"))
 axis(2,las=1,tck=-0.02,hadj=0.8, cex.axis=0.8,at=seq(0,600,100),labels=seq(0,600,100))
 mtext("age [MA]", side=1,line=1.5)
 mtext("number of genera", side=2,line=2.2)

# plot diversity of ecotypes ####################
# for different data "ylim", "at" and "labels" for y-axis may need to be adjusted
plot(-outdiv$TIME,outdiv$MEAN_EDIV,ylim=c(0,35),ylab="",xlab="",axes=F)
 rect(-550,0,0,35,col="gray",border=NA)
 abline(v=c(-250,-65),col="white",lwd=2)
 points(-outdiv$TIME,outdiv$MEAN_EDIV,type="o",col="black",pch=16,cex=0.9)
 axis(1,tck=-0.02,padj=-0.8, cex.axis=0.8,at=seq(-600,0,100),labels=c("600","500","400","300","200","100","0"))
 axis(2,las=1,tck=-0.02,hadj=0.8, cex.axis=0.8,at=seq(0,35,5),labels=seq(0,35,5))
 mtext("age [MA]", side=1,line=1.5)
 mtext("number of ecotypes", side=2,line=2)

# find mean diversity by time bin partitioned into major ecotypes
major<- which(colnames(egen_all)%in% as.character(ecomajors))
important<- egen_all[,major]
others<- apply(egen_all[,-major][,-1],1,sum)
ecogps<- cbind(tbins,important,others)	# diversity partitioned by major ecotypes
ecogps<- rbind(c(tbins[1],ecogps[1,-1]*0),ecogps)
ecogps<- rbind(ecogps,c(tbins[length(tbins)],ecogps[1,-1]*0))
colnames(ecogps)[1]<- "tbins"

# plot diversity of genera #####################
# for different data "ylim", "at" and "labels" for y-axis may need to be adjusted
plot(-ecogps[,1],ecogps[,2],ylim=c(0,600),ylab="",xlab="",axes=F,cex=0)
polygon(-ecogps[,1],apply(ecogps[,2:12],1,sum),col="black",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:11],1,sum),col="gray",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:10],1,sum),col="red3",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:9],1,sum),col="red1",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:8],1,sum),col="orange",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:7],1,sum),col="brown",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:6],1,sum),col="yellow",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:5],1,sum),col="green1",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:4],1,sum),col="green3",border=NA)
polygon(-ecogps[,1],apply(ecogps[,2:3],1,sum),col="blue2",border=NA)
polygon(-ecogps[,1],ecogps[,2],col="blue4",border=NA)
 points(-outdiv$TIME,outdiv$MEAN_DIV,type="o",col="black",pch=16,cex=0.9)
 axis(1,tck=-0.02,padj=-0.8, cex.axis=0.8,at=seq(-600,0,100),labels=c("600","500","400","300","200","100","0"))
 axis(2,las=1,tck=-0.02,hadj=0.8, cex.axis=0.8,at=seq(0,600,100),labels=seq(0,600,100))
 mtext("age [MA]", side=1,line=1.5)
 mtext("number of genera", side=2,line=2.2)
legend(-550,600,rev(colnames(ecogps)[-1]),cex=0.7,lwd=4,
 col=rev(c("blue4","blue2","green3","green1","yellow","brown","orange","red1","red3","gray","black")))

# SAVE OUTPUTS TO EXTERNAL FILES ################
write.csv(outdiv,"Diversity.csv")
write.csv(eocc_all,"EcoOccurrences.csv")
write.csv(egen_all,"EcoGenera.csv")
write.csv(ecogps,"EcoGroup.csv")
