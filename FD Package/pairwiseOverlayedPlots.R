



dataBray<-as.matrix(betapabun$beta.bray)
plotEqual<- matrix(ncol = 2, nrow = (ncol(dataBray)-1))
colnames(plotEqual)<-c("timeContrast", "Bray-Curtis")

for (i in 1:(ncol(dataBray)-1)) {
  plotEqual[i,1] <- paste0(colnames(dataBray)[i], "-",rownames(dataBray)[i+1])
  plotEqual[i,2] <- dataBray[i,(i+1)]
  
}

plotEqual<- as.data.frame(plotEqual)
plotEqual$`Bray-Curtis` <- as.numeric(as.character(plotEqual$`Bray-Curtis`))
plot(plotEqual$timeContrast, plotEqual$`Bray-Curtis`, "l")


write.csv(plotEqual, "plotEqual.csv")
