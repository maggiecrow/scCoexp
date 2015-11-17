plot_figure1 <- function (data, pheno, file_base) {

data=as.matrix(data)
tpm=data*0
for(i in 1:dim(data)[2]){
tpm[,i]=(data[,i]/colSums(data)[i])*10^6
}

pca=prcomp(t(log2(tpm+1)))
z1 = matrix(0,ncol=3,nrow=8)

for(i in 1:8){
z1[i,1]=mean(pca$x[,1][which(pheno$Batch==i)])
z1[i,2]=mean(pca$x[,2][which(pheno$Batch==i)])
}
z1[,3]=10

file1=paste(file_base,"1A.png",sep=".")
file2=paste(file_base,"1B.png",sep=".")
file3=paste(file_base,"1B_inset.png",sep=".")

png(file=file1)
plot(pca$x[,1],pca$x[,2], col=c(rep("blue",each=63),rep("red",each=63)), pch=20,ylab=paste("PC2 (",round((((pca$sdev^2/sum(pca$sdev^2))[2])*100),digits=2),"%)"),xlab=paste("PC1 (",round((((pca$sdev^2/sum(pca$sdev^2))[1])*100),digits=2),"%)"), xlim=c(-200,200),ylim=c(-150,150),axes=F) 

col1=rgb(0,0,1,alpha=0.3)
col2=rgb(1,0,0,alpha=0.3)

symbols(x=z1[1:4,1], y=z1[1:4,2], circles=z1[1:4,3], inches=1/3, ann=F, bg=col1, fg=NULL,xlim=c(-200,200),ylim=c(-150,150),axes=F,add=T)
symbols(x=z1[5:8,1], y=z1[5:8,2], circles=z1[1:4,3], inches=1/3, ann=F, bg=col2, fg=NULL,xlim=c(-200,200),ylim=c(-150,150),axes=F,add=T)

axis(1)
axis(2)
dev.off()

library_complexity=apply(data,2,function(x) sum (x!=0))

png(file=file2)
plot(pca$x[,1],library_complexity, xlab="PC1 loading", ylab="Number of expressed genes",pch=20,axes=F)
text(-150,2000, labels=paste("r = ", round(cor(pca$x[,1],library_complexity),2),sep=""))
axis(1)
axis(2)
dev.off()

png(file=file3)
boxplot(library_complexity~pheno$Batch, boxwex=0.5, whisklty=1, xlab="Batch", ylab="Number of expressed genes", axes=F)
axis(1)
axis(2)
dev.off()
}
