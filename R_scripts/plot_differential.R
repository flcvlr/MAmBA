args <- commandArgs(trailingOnly=T)
####
# usage: plot_differential.R species_anno folder_output_of_sample_1 folder_output_of_sample_2 short_name_sample_1 short_name_sample_2 species output_dir
# 1) folder containing output of main pipeline for sample 1 
# 2) folder containing output of main pipeline for sample 2 
# 3) short name to use in legends for sample 1
# 4) short name to use in legends for sample 2
# 5) annotation_directory (e.g hsa)
# 6) output_dir
meth_data <- read.table(paste(args[6],"/meth_data_to_test",sep=""));
meth_data[6] <- meth_data[6]-meth_data[5]
meth_data[8] <- meth_data[8]-meth_data[7]

ft_FUN = function(x) {fisher.test(matrix(as.integer(x[1:4]), ncol = 2, nrow = 2, byrow = T))$p.value} ;
adj <- p.adjust(apply(meth_data[,5:8],1,ft_FUN),method="bonferroni") ;
all_data <- cbind(meth_data,adj); 
to_plot <- all_data[all_data$adj < 2,]; 

colnames(all_data) <- c("pre-miRNA","start","end","pre-miRNA",paste(args[3],"_5(h)mC",sep=""),paste(args[3],"_C",sep=""),paste(args[4],"_5(h)mC",sep=""),paste(args[4],"_C",sep=""),"adj.P-value")
write.table(all_data[,2:9], file=paste(args[6],"/methylated_miRNAs.txt",sep=""), quote=F,sep ="\t", col.names=T,row.names=F)
for (mirna in unique(to_plot[,1])) {
	sequence <- readChar(paste(args[5],"/anno_for_R/",mirna,".fa",sep=""), file.info(paste(args[5],"/anno_for_R/",mirna,".fa",sep=""))$size-1)
	s_e_q <- strsplit(sequence[1], "")
	pos <- read.table(paste(args[5],"/anno_for_R/",mirna,".pos",sep=""))
	
	sample_1_bed <- read.table(paste(args[1],"/data_for_R/",mirna,".bed",sep=""))
	sample1_cov <- read.table(paste(args[1],"/data_for_R/",mirna,".cov",sep=""))
	max_sample1 <- max(sample1_cov[,3])
	cov1 <- sample1_cov[,3]/max(sample1_cov[,3])

	sample_2_bed <- read.table(paste(args[2],"/data_for_R/",mirna,".bed",sep=""))
	sample2_cov <- read.table(paste(args[2],"/data_for_R/",mirna,".cov",sep=""))
	max_sample2 <- max(sample2_cov[,3])
	cov2 <- sample2_cov[,3]/max(sample2_cov[,3])
	m <- rbind(c(0, 1,0.60, 1), c(0, 1, 0.42,0.60 ),c(0,1,0,0.42))
	par(mar=c(1,1,1,1))
	png(file=paste(args[6],"/",mirna,".png",sep=""),width=1500,height=600)
	silent  <- split.screen(m)

	screen(1)
		par(cex=2,mar=c(0,4,3,1))
		plot(sample_1_bed[,2]-.25,sample_1_bed[,7]+0.001,ylim=c(0.0005,2),type="h",lwd=10,lend=1,xlim=c(0,length(s_e_q[[1]])+2),xlab="Position in pre-miRNA",main =mirna,col="red",ylab="", xaxt="n",yaxt="n", log="y")
		lines(sample_2_bed[,2]+.25,sample_2_bed[,7]+0.001,type="h",col="darkgrey",lwd=10,lend=1)
		rect(0,0.0002,length(s_e_q[[1]])+2, 0.005, density = NULL, border = NA,col = "white")
		axis(2,at=c(0.01,0.1,1), labels=c("0.01","0.1","1.0"), las=2, cex.axis=.8)
		text((length(s_e_q[[1]]))/2,1.1,"Frequency of 5mC")
		for (i in seq(1:(length(s_e_q[[1]])))) {text(i-1,0.0018,s_e_q[[1]][i],cex=1-(length(s_e_q[[1]])-90)/100)}
		segments(pos[1,1]-.5,.0007,pos[1,2]+.5,.0007, col="black",lwd=8,lend=1)
		#for (i in pos[1,1]:pos[1,2]) {text(i,.0015,"-",cex=2)}
		if (length(pos[,1]) > 1 ) {
			segments(pos[2,1]-.5,.0007,pos[2,2]+.5,.0007, col ="black",lwd=8,lend=1)
			#for (i in pos[2,1]:pos[2,2]) {text(i,0.0015,"-",cex=2)}
		}
		for (signif in to_plot[to_plot[,1]==mirna & (to_plot$adj <0.1),2]) { text(signif,1,"*",cex=2)}


	screen(2)
		par(cex=2,mar=c(0,4,0,1))
		plot(sample1_cov[,2],(cov1) ,ylim=c(0,max(cov1)*1.5),type="h", lwd=33-(length(s_e_q[[1]])/45),col="red",lend=1, xlim =c(0,length(s_e_q[[1]])+2), yaxt="n", xaxt="n",ylab="")
		axis(2,at=1, labels=max(sample1_cov[,3]), las=2, cex.axis=.8)
		text((length(s_e_q[[1]]))/2,1.3,"Coverage (raw reads)")
	screen(3)
		par(cex=2,mar=c(5,4,0,1))
		plot(sample2_cov[,2],(cov2) ,ylim=c(0,max(cov2)*1.5),type="h", lwd=33-(length(s_e_q[[1]])/45),col="darkgrey",lend=1, xlim =c(0,length(s_e_q[[1]])+2), yaxt="n", xaxt="n",ylab="")
		axis(2,at=1, labels=max(sample2_cov[,3]), las=2, cex.axis=.8)
		text((length(s_e_q[[1]]))/2,1.3,"Coverage (raw reads)")

	par(xpd=TRUE)
	legend(-3.5,-1.3,c(args[3],args[4],"mature miRNA (mirbase annotation)"), col=c("red","darkgrey","black"), fill=c("red","darkgrey","black"),border=c("red","darkgrey","black"),ncol=3)
	done <- close.screen(all=TRUE)
	fine <- dev.off()
}





