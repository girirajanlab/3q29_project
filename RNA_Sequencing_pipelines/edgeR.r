#This R script implements an edgeR pipeline (v.3.20.1), using the generalized linear model option,
#to identify genes differentially expressed with knockdown of 3q29 homologs in Drosophila.

#Quantified counts were generated from the Tophat2-HTSeq pipeline, and technical replicates were
#merged to produce a final set of quantified counts for each biological replicate.
#These data are available on NCBI GEO, accession number GSE128094.

#Load edgeR library
library("edgeR")

#Load raw counts file
biol_reps_data<-read.table("biological_replicates.txt", sep="\t", stringsAsFactors=FALSE)

#Move row names and column names out of data matrix
rownames(biol_reps_data)<-biol_reps_data[,1]
colnames(biol_reps_data)<-biol_reps_data[1,]
biol_reps_data<-biol_reps_data[,-1]
biol_reps_data<-biol_reps_data[-1,]

#Convert to matrix and convert characters to numeric type
biol_reps_matrix<-as.matrix(biol_reps_data)
class(biol_reps_matrix)<-"numeric"

#Load data matrix into DGEList format to begin edgeR pipepline
groups<-c(rep("NCBP2_FBXO45",3),rep("NCBP2_DLG1",3),rep("NCBP2",3),rep("PAK2",3),rep("FBXO45",3),rep("DLG1",3),rep("KK_Control",3),rep("GD_Control",3))
y<-DGEList(biol_reps_matrix, group=groups)

#Filter variants: Recommended defaults here are 5 counts in smallest library, translated to counts/million (here ~27 million = 0.2 cpm), 
#and in at least 3 samples (or expressed in all members of one group).
keep<-rowSums(cpm(y)>0.2)>=3
y<-y[keep, , keep.lib.sizes=FALSE] #Recalculate library sizes afterwards

#Calculate normalization factor
y<-calcNormFactors(y)

#Create MDS plot to display cliustering of samples relative to each other. Save as PDF.
plotMDS(y)
title(main="Multi-dimensional scaling plot of knockdown model sample clustering")

#Create the design matrix for all samples
design<-model.matrix(~ 0 + groups)
#Estimate dispersion of variances along model
y <- estimateDisp(y, design, robust=TRUE)
#Model the quasi-lielihood dispersions of variation within model
fit <- glmQLFit(y, design, robust=TRUE)

#Export normalized expression counts for downstream analysis
effectiveLibSize<-y$samples$lib.size*y$samples$norm.factors
normCounts<-log2(t(t(y$counts+0.5) / (effectiveLibSize+0.5)))
write.csv(normCounts,file="edgeR_norm_log_counts.csv")

#Perform differential comparison tests. Controls are wild-type flies with same VDRC genetic background (GD or KK)
NCBP2<-glmQLFTest(fit, contrast=c(0,0,0,-1,1,0,0,0))
PAK2<-glmQLFTest(fit, contrast=c(0,0,0,-1,0,0,0,1))
FBXO45<-glmQLFTest(fit, contrast=c(0,1,-1,0,0,0,0,0))
DLG1<-glmQLFTest(fit, contrast=c(1,0,-1,0,0,0,0,0))
NCBP2_FBXO45<-glmQLFTest(fit, contrast=c(0,0,0,-1,0,0,1,0))
NCBP2_DLG1<-glmQLFTest(fit, contrast=c(0,0,0,-1,0,1,0,0))

#Export results for each differential comparison to csv files
write.csv(topTags(NCBP2, n=nrow(NCBP2$table))$table, file="NCBP2_edgeR.csv")
write.csv(topTags(PAK2, n=nrow(PAK2$table))$table, file="PAK2_edgeR.csv")
write.csv(topTags(FBXO45, n=nrow(FBXO45$table))$table, file="FBXO45_edgeR.csv")
write.csv(topTags(DLG1, n=nrow(DLG1$table))$table, file="DLG1_edgeR.csv")
write.csv(topTags(NCBP2_DLG1, n=nrow(NCBP2_DLG1$table))$table, file="NCBP2_DLG1_edgeR.csv")
write.csv(topTags(NCBP2_FBXO45, n=nrow(NCBP2_FBXO45$table))$table, file="NCBP2_FBXO45_edgeR.csv")
