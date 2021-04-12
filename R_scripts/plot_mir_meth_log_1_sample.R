args <- commandArgs(trailingOnly=T)
### 1) miR_name
### 2) sample_output_dir
### 3) annotation_folder
### 4) sample_short_name
### 5) color
### 6) log/linear
#args




sample_1_bed <- read.table(paste(args[2],"/data_for_R/",args[1],".bed",sep=""))

sequenza <- readChar(paste(args[3],"/anno_for_R/",args[1],".fa",sep=""), file.info(paste(args[3],"/anno_for_R/",args[1],".fa",sep=""))$size-1)
s_e_q <- strsplit(sequenza[1], "")
pos <- read.table(paste(args[3],"/anno_for_R/",args[1],".pos",sep=""))
sample1_cov <- read.table(paste(args[2],"/data_for_R/",args[1],".cov",sep=""))
max_sample1 <- max(sample1_cov[,3])
if (max_sample1 > 29) {
cov1 <- sample1_cov[,3]/max(sample1_cov[,3])
sig <- sample_1_bed[sample_1_bed[,5] < 0.01,]

m <- rbind(c(0, 1,0.45, 1), c(0, 1, 0,0.45 ))
par(mar=c(1,1,1,1))

if(min(sample_1_bed[,5]) < 0.01) {status <- "methylated/"}
else {status <- "not_methylated/"}
png(file=paste(args[2],"/images/",status,args[4],"_",args[1],".png",sep=""),width=1500,height=500)

silent  <- split.screen(m)
screen(1)
par(cex=2,mar=c(0,4,3,1))

if(args[6]=="log") {
	plot(sample_1_bed[,2],sample_1_bed[,7]+0.001,ylim=c(0.001,4),type="h",lwd=10,lend=1,xlim=c(0,length(s_e_q[[1]])+2),xlab="Position in pre-miRNA",main =args[1],col=args[5],ylab="", xaxt="n",yaxt="n", log="y")
	rect(-.5,0.00072,length(s_e_q[[1]])+2, 0.01, density = NULL, border = NA,col = "white")
	axis(2,at=c(0.01,0.1,1), labels=c("0.01","0.1","1.0"), las=2, cex.axis=.8)
	text((length(s_e_q[[1]]))/2,2.1,"non-conversion frequency")
	for (i in seq(1:(length(s_e_q[[1]])))) {text(i-1,0.005,s_e_q[[1]][i],cex=1-(length(s_e_q[[1]])-90)/100)}
	segments(pos[1,1]-.5,.002,pos[1,2]+.5,.002, col="black",lwd=8,lend=1)
	if (length(pos[,1]) > 1 ) {
		segments(pos[2,1]-.5,.002,pos[2,2]+.5,.002, col ="black",lwd=8,lend=1)
	}

	for (i in sig[,2]) {text(i,1.25,"*",cex=1.5)}



} else {plot(sample_1_bed[,2],sample_1_bed[,7],ylim=c(-0.5,1.35),type="h",lwd=10,lend=1,xlim=c(0,length(s_e_q[[1]])+2),xlab="Position in pre-miRNA",main =args[1],col=args[5],ylab="", xaxt="n",yaxt="n")
	#rect(-.5,0.00072,length(s_e_q[[1]])+2, 0.01, density = NULL, border = NA,col = "white")
	axis(2,at=c(0,0.5,1), labels=c("0","0.5","1"), las=2, cex.axis=.8)
	text((length(s_e_q[[1]]))/2,1.25,"non-conversion frequency")
	for (i in seq(1:(length(s_e_q[[1]])))) {text(i-1,-0.2,s_e_q[[1]][i],cex=1-(length(s_e_q[[1]])-90)/100)}
	segments(pos[1,1]-.5,-.4,pos[1,2]+.5,-.4, col="black",lwd=8,lend=1)
	if (length(pos[,1]) > 1 ) {
		segments(pos[2,1]-.5,-.4,pos[2,2]+.5,-.4, col ="black",lwd=8,lend=1)
	}

	for (i in sig[,2]) {text(i,1.1,"*",cex=1.5)}



}	





#text(30,1.5,paste("relative_coverage; max reads in input= ",max_input,"; max reads in AGO1=", max_ago,sep=""))
screen(2)
par(cex=2,mar=c(5,4,0,1))
plot(sample1_cov[,2],(cov1) ,ylim=c(0,max(cov1)*1.5),type="h", lwd=33-(length(s_e_q[[1]])/45),col=args[5],lend=1, xlim =c(0,length(s_e_q[[1]])+2), yaxt="n", xaxt="n",ylab="")
axis(2,at=1, labels=max(sample1_cov[,3]), las=2, cex.axis=.8)
text((length(s_e_q[[1]]))/2,1.1,"Coverage (raw reads)")
par(xpd=TRUE)
legend(-3.5,-1.3,c(args[4],"mature miRNA (mirbase annotation)"), col=c(args[5],"black"), fill=c(args[5],"black"),border=c(args[5],"black"),ncol=3)
text((length(s_e_q[[1]]))/2,1.1,"Coverage (raw reads)")


done <- close.screen(all=TRUE)
fine <- dev.off()
}



