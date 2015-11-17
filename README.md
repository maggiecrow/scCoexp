##Experimental design and analysis for single cell RNA-seq co-expression
Megan Crow, Anirban Paul, Z. Josh Huang, Jesse Gillis

####Abstract
Single cell RNA-sequencing (scRNA-seq) is becoming increasingly widespread, yielding vast amounts of molecular data that will need to be interpreted. Co-expression analysis is a good technique to do this, but it is unclear how technical variation associated with scRNA-seq will influence results. To address this, we sequenced two well-characterized adult mouse inhibitory interneuron cell-types in an experimental design which allows us to detect the influence of technical properties on co-expression networks. We find that normalization induces strong co-variation among genes in scRNA-seq co-expression networks, which leads to altered network topology and reduced functional connectivity, quantified through cross-validation to learn Gene Ontology functions. Network aggregation across library batches increases performance modestly for most functions (mean AUROC across individual networks=0.538 +/- 0.0047, mean aggregate AUROC=0.556 +/- 0.0073, p<0.05, Wilcoxon rank sum test, n=108 GO groups), but yields much higher performance for synaptic genes (AUROC=0.793). However, we discovered that performance is highly dependent on gene expression level and controlling for expression reduces learnability for synaptic genes to baseline levels (mean AUROC=0.53 +/- 0.02). Finally, we demonstrate how the lessons learned from co-expression can improve the interpretation of differential expression by discriminating between robust and fragile expression profile co-variation.  Our results provide strong evidence that scRNA-seq co-expression networks are highly sensitive to normalization, and that expression level is a primary feature driving gene-gene connectivity in functional networks. We encourage researchers to explicitly control for technical variation in their experimental design, and to use meta-analysis in combination with expression-level matching and thresholding to assess the robustness of co-expression results.

####This site contains scripts to reproduce analyses in this manuscript
* Expression data and metadata available upon request.
* The generic GO slim was downloaded from the Gene Ontology Consortium and filtered for mammalian functions (http://geneontology.org/page/go-slim-and-subset-guide, access date August 2015). 
* The synaptic geneset was downloaded from the Genes to Cognition database (https://www.genes2cognition.org/db/GeneList/L00000001, access date February 2015).
* assort.R and run_GBA.Rdata were provided by Sara Ballouz

######under construction
