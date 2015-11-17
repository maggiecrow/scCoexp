sunset_plot_and_topology_analysis_pairwise <- function (data, pheno, file_ext) {

data=as.matrix(data)
temp=data
temp[data==0]=NA
nd_exp <-matrix(0,ncol=1,nrow=8)
assort_mat <-matrix(0,ncol=1,nrow=8)
source("assort.r")
library(RColorBrewer)
col1=colorRampPalette(brewer.pal(11,"RdYlBu"))(100)
col2=rev(col1)

for (i in 1:length(levels(as.factor(pheno$Batch)))){
		
		#count the number of zero values for each gene and filter out genes that are always 0
		zeros=apply(data[,which(pheno$Batch==i)],1, function(x) sum(x==0))
		filt=zeros!=length(which(pheno$Batch==i))
		ord=order(zeros[filt])

		#generate the coexpression network using Spearman correlation (NB: slow)
		coexp=abs(cor(t(temp[filt,which(pheno$Batch==i)]),use="pairwise.complete.obs",method="s"))

		#convert correlation coefficients to ranks and standardize
		rank.coexp=coexp
		rank.coexp[]=rank(coexp,ties.method="average",na.last="keep")
		rank.coexp=rank.coexp/length(which(!is.na(rank.coexp)))

		#order the network by the number of zeros each gene takes
		rank.coexp=rank.coexp[ord,ord]

		#calculate assortativity
		assort_mat[i,1]=assortativity(rank.coexp, pairwise=T)

		#calculate the correlation between node degree and median expression
		nd=apply(rank.coexp,1,function(x) median(x,na.rm=T))
    		med_exp=apply(data[filt,which(pheno$Batch==i)],1,function(x) median(x))
		nd_exp[i,1]=cor(nd,med_exp,method="spearman")
		
		#sunset plot
		file_plot=paste(file_ext,i,"sunset.png",sep=".")
		breaks=seq(min(rank.coexp,na.rm=T),max(rank.coexp,na.rm=T),length=101)
		x=rep(c(rep(0,each=9),1),length.out=dim(rank.coexp)[1])
		png(file_plot)
		image(rank.coexp[which(x==1),which(x==1)],col=col2,breaks=breaks,useRaster=T)
		dev.off()
	}

write.table(nd_exp,file=paste(file_ext,"nd_vs_exp.txt",sep="."),sep="\t")
write.table(assort_mat,file=paste(file_ext,"assortativity.txt",sep="."),sep="\t")

}