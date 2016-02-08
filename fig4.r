## Figure 4A and 4B - network aggregation method 
sample_aggs_get_ROC <- function (go, translate, file_base, num_reps, file_ext){

#make a table to store AUROC scores
AUROC=matrix(0,nrow=num_reps,ncol=7)

#randomly sample networks and aggregate
for (n in 2:7){
sample_table=matrix(0,ncol=num_reps,nrow=n)

	for (i in 1:num_reps){
	sample_n=sample(c(1:8),n)

	agg<-matrix(0,ncol=17457,nrow=17457)

		for (j in 1:length(sample_n)){
		filename=paste(file_base,sample_n[j],"rank_coexp.Rdata",sep=".")
		#load a randomly chosen network (rank.coexp)
		load(filename)
		#add to aggregate
		agg=agg+rank.coexp
		}

	#divide the aggregate by the number of networks added
	agg=agg/n

	#change rownames from gene symbols to Entrez IDs
	m<-match(rownames(rank.coexp),translate[,1])
	f.a=!is.na(m)
	f.b=m[f.a]
	agg=agg[f.a,f.a]
	rownames(agg)=translate[f.b,2]
	colnames(agg)=translate[f.b,2]

	#run the neighbor voting algorithm 
	GBA<-run_GBA(agg,go)
	GBA_samples<-list(GBA,list(sample_n))
	file2=paste(file_ext,n,"net",i,"GBA.Rdata",sep=".")
	
	#save the results
	save(GBA_samples,file=file2)
	sample_table[,i]=sample_n

	#store the AUROC for each aggregate network in the AUROC table
	AUROC[i,n]=GBA[[3]]
	}
#save the batch IDs for each aggregate network
write.table(sample_table,file=paste(file_ext,n,"samples.txt",sep="."),sep="\t")
}
write.table(AUROC,file=paste(file_ext,"AUROCs.txt",sep="."),sep="\t")
}

## Figure 4C - performance vs expression level plot 
nfold_ROC <- function (data, pheno, geneset, file_base, file_ext) {
data=as.matrix(data) 
geneset=as.matrix(geneset)

for (i in 1:length(levels(as.factor(pheno$Batch)))) {
	filename=paste(file_base,i,"rank_coexp.Rdata",sep=".")
	load(filename)
	
	#get median expression data for the batch
	temp=data[,which(pheno$Batch==i)]
	medians=apply(temp,1,function (x) median(x))
	medians=as.matrix(medians)
	rownames(medians)=rownames(temp)

	#run n-fold neighbor voting
	rocs=neighbor.voting.CV.nfold(as.matrix(geneset[,1]),rank.coexp,sum(geneset))
	rocs=as.matrix(rocs)
	rownames(rocs)=rownames(geneset)[which(geneset==1)]

	#create a table with genes as rows and medians and rocs as columns
	m<-match(rownames(medians),rownames(rocs))
	f.a=!is.na(m)
	f.b=m[f.a]
	temp2=cbind(medians[f.a,],rocs[f.b,])

	#plot rocs vs log2 median expression 
	png(file=paste(file_ext,i,"png",sep="."))
	plot(log2(temp2[,1]+1),temp2[,2],xlab="log2 median count + 1", ylab="ROC")
	mtext(paste("Spearman=",round(cor(temp2[,2],temp2[,1],method="s"),2),", Pearson=",round(cor(temp2[,2],temp2[,1]),2),sep=""),adj=1)
	dev.off()
	file_out=paste(file_ext,i,"txt",sep=".")

	#save roc and median expression table
	write.table(temp2,file=file_out,sep="\t")
	}
}

# Figure 4D - thresholding networks and assessing performance 
threshold_GBA <- function (data, pheno, geneset, file_base, file_ext) {
data=as.matrix(data) 
geneset=as.matrix(geneset)
GBA=matrix(0,ncol=3,nrow=length(levels(as.factor(pheno$Batch))))

for (i in 1:length(levels(as.factor(pheno$Batch)))) {
	filename=paste(file_base,i,"rank_coexp.Rdata",sep=".")
	load(filename)

	#calculate gene expression median 
	temp=data[,which(pheno$Batch==i)]
	medians=apply(temp,1,function (x) median(x))
	medians=as.matrix(medians)
	rownames(medians)=rownames(temp)

	#filter networks and geneset to contain only genes with expression >16 counts
	filt=medians>=16
	rank.coexp=rank.coexp[filt,filt]
	geneset2=as.matrix(geneset[filt,])

	#run the neighbor voting algorithm and store the AUROC 
	GBA[i,]=neighbor.voting.CV(geneset2,rank.coexp,3)
	}
write.table(GBA,file=file_ext,sep="\t")
}


##### HELPER FUNCTIONS -- run_GBA.Rdata (includes filter_network_cols, neighbor.voting.CV and predictions) and neighbor.voting.CV.nfold #######

run_GBA <- function (network, labels) {

m = match(rownames(network), rownames(labels))
f = !is.na(m)
g = m[f]

network.sub  = network[f,f]
genes.labels = filter_network_cols(labels[g,],20,1000)

roc.sub   = neighbor.voting.CV(genes.labels, network.sub, 3)
genes = predictions(genes.labels, network.sub)
auroc = mean(roc.sub[,1], na.rm=T)

results = list(roc.sub, genes, auroc)
return(results)

}

filter_network_cols <- function(network, min, max)
{
	# Filter terms between min and max
	colsums <- colSums(network)
	col.filter <- which( (colsums > min & colsums < max) )

	network <- network[,col.filter]
	
        return(network)
}

neighbor.voting.CV <- function(genes.labels,network,nFold){

        # Filter for common genes between network and labels
        ord = order(rownames(network))
        network = network[ord,ord]

        ord = order(rownames(genes.labels))
        genes.labels = genes.labels[ord,]

        match.lab <- match(rownames(genes.labels),rownames(network))
        filt.lab  <- !is.na(match.lab)
        filt.net  <- match.lab[filt.lab]
        network   <- network[filt.net,filt.net]
        genes.labels  <- genes.labels[filt.lab,]

        # genes.label : needs to be in 1s and 0s
        l <- dim(genes.labels)[2]
        g <- dim(genes.labels)[1]
        ab <- which(genes.labels != 0, arr.ind=T)
        n  <- length(ab[,1])

        #print("Make genes label CV matrix")
        test.genes.labels = matrix(genes.labels, nrow=g, ncol=nFold*l)

        # For each fold in each GO group, remove 1/nth of the values of the genes.label
        for (j in 1:l){
                d <- which(ab[,2]==j)  # Which indices the genes are in this particular GO group
                t <- length(d)         # Total number of genes in the GO group
                r <- sample(1:t, replace=F)
                f <- t/nFold
                for (i in 1:nFold){
                        e = c( (1:f) + f*(i-1) )
                        e = sort(r[e])
                        c = j + l*(i-1)        # GO group to look at (ie column)
                        test.genes.labels[ab[d],c][e] <- 0
                }
        }

        #print("Get sums - mat. mul.")
        #sumin    = ( t(network) %*% test.genes.labels)
        sumin    = ( (network) %*% test.genes.labels)

        #print("Get sums - calc sumall")
        sumall   = matrix(apply(network,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])

        #print("Get sums - calc predicts")
        predicts = sumin/sumall

        #predicts0 = predicts
        #predicts = get_subnet_rows(predicts, rownames(genes.labels))

        #print("Hide training data")
        nans = which(test.genes.labels == 1, arr.ind=T)
        predicts[nans] <- NA

        #print("Rank test data")
        predicts = apply(abs(predicts), 2, rank,na.last="keep",ties.method="average")
        filter = matrix(genes.labels, nrow=g, ncol=nFold*l)
        negatives = which(filter == 0, arr.ind=T)
        positives = which(filter == 1, arr.ind=T)
        predicts[negatives] <- 0

        #print("Calculate ROC - np")
        np = colSums(filter) - colSums(test.genes.labels) # Postives

        #print("Calculate ROC - nn")
        nn = dim(test.genes.labels)[1] - colSums(filter)      # Negatives

        #print("Calculate ROC - p")
        p =  apply(predicts,2,sum,na.rm=T)

        #print("Calculate ROC - rocN")
        rocN = (p/np - (np+1)/2)/nn
        rocN = matrix(rocN, ncol=nFold, nrow=l)
        rocN = rowMeans(rocN)

        #print("Calculate node degree")
        node_degree <- rowSums(network)
        colsums <- colSums(genes.labels)

        #print("Calculate node degree - sum across gene labels")
        node_degree <- matrix(node_degree)
        temp <- t(node_degree) %*% genes.labels

        #print("Calculate node degree - average")
        average_node_degree <- t(temp)/colsums

        #print("Calculate node degree roc - rank node degree")
        ranks <- apply(abs(node_degree), 2, rank,na.last="keep",ties.method="average")
        ranks <- matrix(ranks, nrow=length(ranks), ncol=dim(genes.labels)[2])

        #print("Calculate node degree roc - remove negatives")
        negatives = which(genes.labels == 0, arr.ind=T)
        ranks[negatives] = 0

        #print("Calculate node degree roc - np")
        np = colSums(genes.labels)

        #print("Calculate node degree roc - nn")
        nn = dim(genes.labels)[1] - np

        #print("Calculate node degree roc - p")
        p =  apply(ranks,2,sum,na.rm=T)

        #print("Calculate node degree roc - roc")
        roc  <-  (p/np - (np+1)/2)/nn

        scores = cbind(rocN,matrix(average_node_degree)[,1],roc)

}

predictions <- function(genes.labels, network)
{

        # Filter for common genes between network and labels
	ord = order(rownames(network))
	network = network[ord,ord]

        ord = order(rownames(genes.labels))
        genes.labels = genes.labels[ord,]

	match.lab <- match(rownames(genes.labels),rownames(network))
       	filt.lab  <- !is.na(match.lab)
	filt.net  <- match.lab[filt.lab]
        network   <- network[filt.net,filt.net]
        genes.labels  <- genes.labels[filt.lab,]

        # genes.label : needs to be in 1s and 0s
        l <- dim(genes.labels)[2]
        g <- dim(genes.labels)[1]
        ab <- which(genes.labels != 0, arr.ind=T)
        n  <- length(ab[,1])

        # no self predictions
        diag(network) <- 0

        # Get sums - mat. mul.
        sumin    = ( t(network) %*% genes.labels)

        # Sum of all edge in network
        sumall   = matrix(apply(network,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])

        # Predictions
        predicts = sumin/sumall
        return(predicts)
}

neighbor.voting.CV.nfold<-function(genes.labels,network,nFold){

        # Filter for common genes between network and labels
        ord = order(rownames(network))
        network = network[ord,ord]

        ord = order(rownames(genes.labels))
        genes.labels = as.matrix(genes.labels[ord,])

        match.lab <- match(rownames(genes.labels),rownames(network))
        filt.lab  <- !is.na(match.lab)
        filt.net  <- match.lab[filt.lab]
        network   <- network[filt.net,filt.net]
        genes.labels  <- as.matrix(genes.labels[filt.lab,])

        # genes.label : needs to be in 1s and 0s
        l <- dim(genes.labels)[2]
        g <- dim(genes.labels)[1]
        ab <- which(genes.labels != 0, arr.ind=T)
        n  <- length(ab[,1])

        #print("Make genes label CV matrix")
        test.genes.labels = matrix(genes.labels, nrow=g, ncol=nFold*l)

	for (i in 1:nFold){
		d<-which(ab[,2]==1)
		test.genes.labels[ab[d],i][d][i]<-0
		}
		
        #print("Get sums - mat. mul.")
        sumin    = ( (network) %*% test.genes.labels)

        #print("Get sums - calc sumall")
        sumall   = matrix(apply(network,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])

        #print("Get sums - calc predicts")
        predicts = sumin/sumall

        #print("Hide training data")
        nans = which(test.genes.labels == 1, arr.ind=T)
        predicts[nans] <- NA

        #print("Rank test data")
        predicts = apply(-abs(predicts), 2, rank,na.last="keep",ties.method="average")
        filter = matrix(genes.labels, nrow=g, ncol=nFold*l)
        negatives = which(filter == 0, arr.ind=T)
        positives = which(filter == 1, arr.ind=T)
        predicts[negatives] <- 0

        #print("Calculate ROC - np")
        np = colSums(filter) - colSums(test.genes.labels) # Postives

        #print("Calculate ROC - nn")
        nn = dim(test.genes.labels)[1] - colSums(filter)      # Negatives

        #print("Calculate ROC - p")
        p =  apply(predicts,2,sum,na.rm=T)

        #print("Calculate ROC - rocN")
        rocN = 1-(p/np - (np+1)/2)/nn
        scores = rocN
}
