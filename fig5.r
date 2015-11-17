### Figure 5C and 5D
make_networks_and_plot_ROCs <- function(data, pheno, DE_genes, filename) {  
load("run_GBA.Rdata")
library(gplots)
library(RColorBrewer)
data=as.matrix(data)
DE_genes=as.matrix(DE_genes)
batches=levels(as.factor(pheno$Batch))
combinations=combn(batches,2)	
GBA_mat_naive<-matrix(0,ncol=8,nrow=8)
GBA_mat_robust<-matrix(0,ncol=8,nrow=8)

#make networks for each combination and save ROCs 
for (i in 1:dim(combinations)[2]){
		batch_ind1=which(pheno$Batch==combinations[1,i])
		batch_ind2=which(pheno$Batch==combinations[2,i])
		temp=data[,c(batch_ind1,batch_ind2)]
		coexp=cor(t(temp),method="spearman")
		rank.coexp<-coexp
		rank.coexp[]<-rank(coexp,na.last="keep",ties.method="average")
		rank.coexp=rank.coexp/length(which(!is.na(rank.coexp)))
		nas=which(is.na(rank.coexp),arr.ind=T)
		rank.coexp[nas]=median(rank.coexp,na.rm=T)
		temp2=neighbor.voting.CV(DE_genes,rank.coexp,3)
		GBA_mat_naive[batch_ind1,batch_ind2]=temp2[1,1]
		GBA_mat_robust[batch_ind1,batch_ind2]=temp2[2,1]
	}

#make networks for each batch and save ROCs
for (i in 1:8){
		temp=data[,which(pheno$Batch==i)]
		coexp=cor(t(temp),method="spearman")
		rank.coexp<-coexp
		rank.coexp[]<-rank(coexp,na.last="keep",ties.method="average")
		rank.coexp=rank.coexp/length(which(!is.na(rank.coexp)))
		nas=which(is.na(rank.coexp),arr.ind=T)
		rank.coexp[nas]=median(rank.coexp,na.rm=T)
		temp2=neighbor.voting.CV(DE_genes,rank.coexp,3)
		GBA_mat_naive[i,i]=temp2[1,1]
		GBA_mat_robust[i,i]=temp2[2,1]
	}

#make ROC matrices symmetrical
ind=lower.tri(GBA_mat_naive,diag=F)
GBA_mat_naive[ind]=t(GBA_mat_naive)[ind]
GBA_mat_robust[ind]=t(GBA_mat_robust)[ind]

#plot
col1=colorRampPalette(brewer.pal(11,"RdYlBu"))(100)
col2=rev(col1)
breaks=seq(0.4,0.8,length=101)

file1=paste(filename,"naiveDE_heatmap.png",sep=".")
file2=paste(filename,"robustDE_heatmap.png",sep=".")

png(file=file1)
heatmap.2(GBA_mat_naive,trace="none",density.info="none",keysize=1,dendrogram="none",Rowv=FALSE,Colv=FALSE,col=col2, breaks=breaks)
dev.off()

png(file=file2)
heatmap.2(GBA_mat_robust,trace="none",density.info="none",keysize=1,dendrogram="none",Rowv=FALSE,Colv=FALSE,col=col2, breaks=breaks)
dev.off()

}


######### HELPER FUNCTIONS -- batch_naive_DE, batch_robust_DE to get gene sets for network analysis; run_GBA can be found in fig4.r  ###########

batch_naive_DE <- function(data,pheno){

data<-as.matrix(data)
rank.dat=apply(data,2,function(x) rank(x,ties.method="random"))
pval<-matrix(0,ncol=1,nrow=dim(data)[1])
rownames(pval)=rownames(data)

for (i in 1:dim(pval)[2]){
	celltype_ind=which(pheno$Cell==j)
	other_ind=which(pheno$Cell!=j)
		for (j in 1:dim(data)[1]){
			pval[j,i]=wilcox.test(rank.dat[j,celltype_ind],rank.dat[j,other_ind])$p.value
		}
	}

qval=as.matrix(p.adjust(pval[,1],method="BH"))

fc=matrix(0,ncol=length(levels(as.factor(pheno$Celltype))),nrow=dim(data)[1])
for (j in 1:dim(fc)[2]){
	celltype_ind=which(pheno$Cell==j)
	other_ind=which(pheno$Cell!=j)
	for (i in 1:dim(data)[1]){
	fc[i,j]=mean(data[i,celltype_ind])/mean(data[i,other_ind])
	}
}

DE_matrix=cbind(qval,fc[,1])
DE_genes=rownames(data)[DE_mat[,1]<=0.05&DE_mat[,2]>=4]
return(DE_genes)

}

batch_robust_DE <- function(data,pheno){
## requires fold_change, metap and thresholding functions (below)

data<-as.matrix(data)
pval<-matrix(0,ncol=4,nrow=dim(data)[1])
rownames(pval)=rownames(data)
  
for (i in 1:dim(pval)[2]){
    batch_ind1=which(pheno$Batch==i)
    batch_ind2=which(pheno$Batch==i+4)
    for (j in 1:dim(data)[1]){
      pval[j,i]=wilcox.test(data[j,batch_ind1],data[j,batch_ind2])$p.value
    }
  }

qvals=pval*0
for(i in 1:dim(pval)[2]){
qvals[,i]=p.adjust(pval[,i],method="BH")
}

pval.t<-thresholding(pval,qvals)
meta.q<-metaP_pairwise(pval.t)
fc<-fold_change(data,pheno)
DE_matrix=cbind(meta.q,fc[,1])
DE_genes=rownames(data)[DE_mat[,1]<=0.05&DE_mat[,2]>=4]
return(DE_genes)
}

### run functions and make DE gene list for GBA

DE_naive=batch_naive_DE(data,pheno)
DE_robust=batch_robust_DE(data,pheno)
DE_genes=matrix(0,ncol=2,nrow=dim(data)[1])
rownames(DE_genes)=rownames(data)
m<-match(rownames(DE_genes),DE_naive)
DE_genes[!is.na(m),1]=1
m<-match(rownames(DE_genes),DE_robust)
DE_genes[!is.na(m),2]=1
colnames(DE_genes)=c("Naive","Robust")

############# DE HELPER FUNCTIONS ##############

fold_change <- function (data,pheno) {
fc=matrix(0,ncol=length(levels(as.factor(pheno$Celltype))),nrow=dim(data)[1])
for (j in 1:dim(fc)[2]){
	celltype_ind=which(pheno$Cell==j)
	other_ind=which(pheno$Cell!=j)
	for (i in 1:dim(data)[1]){
	fc[i,j]=mean(data[i,celltype_ind])/mean(data[i,other_ind])
	}
}
return(fc)
}

thresholding <- function (pval, qval) {
  for(i in 1:dim(qval)[2]){
    qval_thresh=rownames(qval[qval[,i]<=0.05,])
    m<-match(rownames(pval),qval_thresh)
    x=pval[!is.na(m),]
    y=which.max(x[,i])
    pval[pval[,i]<=x[y,i],i]=x[y,i]
  }
  return(pval)
}


##meta-analysis code from the MADAM library on CRAN (deprecated; citation - Kugler KG, Mueller LA, Graber A. MADAM - An open source meta-analysis toolbox for R and Bioconductor. Source Code for Biology and Medicine. 2010;5:3. doi:10.1186/1751-0473-5-3.) 
fisher.method <- function(pvals, method=c("fisher"), p.corr=c("bonferroni","BH","none"), zero.sub=0.00001, na.rm=FALSE, mc.cores=NULL){
  stopifnot(method %in% c("fisher"))
  stopifnot(p.corr %in% c("none","bonferroni","BH"))
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  if(is.null(dim(pvals)))
    stop("pvals must have a dim attribute")
  p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
  ##substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  if(is.null(mc.cores)){
    fisher.sums <- data.frame(do.call(rbind, apply(pvals, 1, fisher.sum, zero.sub=zero.sub, na.rm=na.rm)))
  } else {
    fisher.sums <- multicore::mclapply(1:nrow(pvals), function(i){
      fisher.sum(pvals[i,], zero.sub=zero.sub, na.rm=na.rm)
    }, mc.cores=mc.cores)
    fisher.sums <- data.frame(do.call(rbind, fisher.sums))
  }
    
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S, df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
                              bonferroni = p.adjust(fisher.sums$p.value, "bonferroni"),
                              BH = p.adjust(fisher.sums$p.value, "BH"),
                              none = fisher.sums$p.value)
  return(fisher.sums)
}

### get meta-analytic BH adjusted p-values
metaP_pairwise<-function(pvals) {
  metaP=matrix(0,ncol=1,nrow=dim(pvals)[1])
  rownames(metaP)=rownames(pvals)
  
  t_pval=t(pvals)	
  meta_p<-matrix(0,ncol=dim(t_pval)[2],nrow=1)
  colnames(meta_p)=colnames(t_pval)
  for (j in 1:dim(t_pval)[2]){
    y<-as.matrix(t_pval[,j])
    meta_p[,j]=fisher.method(t(y),zero.sub=1e-08)$p.value
  }
  metaP[,1]<-t(meta_p)

  qval<-metaP*0
  rownames(qval)=rownames(metaP)
  
  for (i in 1:dim(metaP)[2]){
    qval[,i]<-p.adjust(metaP[,i],method="BH")
  }
return(qval)
}
