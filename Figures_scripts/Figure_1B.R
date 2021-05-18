mydata <- data.frame()
colnames(mydata) <- c("time", "memory")
rownames(mydata) <- c("ROTS", "DiffBind","MAnorm2", "diffReps","PePr", "THOR")

mydata["ROTS","time"] <- 31
mydata["ROTS","memory"] <- 1.23
mydata["DiffBind","time"] <- 100
mydata["DiffBind","memory"] <- 42.35
mydata["MAnorm2","time"] <- 507
mydata["MAnorm2","memory"] <- 15.35
mydata["diffReps","time"] <- 694
mydata["diffReps","memory"] <- 2.11
mydata["PePr","time"] <- 347
mydata["PePr","memory"] <- 2.35
mydata["THOR","time"] <- 164
mydata["THOR","memory"] <- 10.05

mydata2 <- mydata
mydata2$memory<- mydata2$memory+3
mydata2["DiffBind","memory"]<-mydata2["DiffBind","memory"]-3
mydata2["DiffBind","time"]<-mydata2["DiffBind","time"]+110
mydata2["diffReps","time"]<-mydata2["diffReps","time"]-70
mydata2["ROTS","time"]<-mydata2["ROTS","time"]+50

png(filename = "time_vs_memory.png",height = 600,width = 600)
par(mar=c(4,4,2,2))
plot(mydata$time,mydata$memory,
     pch=19, col=c("black","blue","violet","green","red","orange"),
     cex=2, cex.lab=2, xlab="",
     ylab="", main="",cex.main=2,axes=FALSE, frame=TRUE)
text(mydata2$time,mydata2$memory,c("ROTS","DiffBind","MAnorm2","diffReps","PePr","THOR"),cex=2)
axis(1,seq(100,700,100),cex.axis=2)
axis(2,seq(0,40,10),cex.axis=2)
#mtext(side = 1,text = "time in minutes", cex=2,line =0.1)
#mtext(side = 2,text = "memory in Gigabites", cex=2,line = 0.1)
dev.off()
