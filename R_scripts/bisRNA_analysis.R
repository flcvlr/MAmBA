args <- commandArgs(trailingOnly=T)
### 1) sample_output_dir

library("BisRNA");

aa <- read.table(paste(args[1],"/bisRNA_input",sep=""))
colnames(aa) <- c("RNA","Cpos","coverage","ncratio")
lambda_1 <- RNAmeth.poisson.par(aa)$estimate
output <- RNAmeth.poisson.test(aa,lambda_1)

final_out <- data.frame(output$nonconv.ratio, output$pv.adj, row.names=output$RNA.pos)
signif <- subset(final_out, (final_out$output.pv.adj < 0.01 & final_out$output.nonconv.ratio > 0.2))
write.table(final_out,sep = "\t",quote=F, col.names = F,file=paste(args[1],"/bisRNA_output",sep =""))












