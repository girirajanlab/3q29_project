#This R script first finds the overlap of genes that are annotated for apoptosis function (apoptotic process, GO:0006915 or children terms) 
#among three sets of candidate neurodevelopmental (NDD) genes: schizophrenia (derived from Purcell et al, Nature 2014); autism (derived from SFARI Gene database);
#and intellectual disablity (derived from DD G2P database). It then performs 100,000 simulations of randomly selecting the same number of genes in each 
#NDD geneset from the set of all RefSeq genes, and determining the number of genes with apoptosis in each random set . 

#A list of all apoptosis genes is provided in the file human_apoptotic_process_raw.txt, a list of all RefSeq genes is provided 
#in the file refseq_genes.txt, and lists of the candidate genes are provided in the files schiz_genes.txt, ASD_genes.txt, and refseq_genes.txt. 
#See README.MD for more information on these gene sets.

#Open gene lists
apoptosis<-read.table("human_apoptotic_process_raw.txt",sep="\t")
schiz<-read.table("schiz_genes.txt")
asd<-read.table("ASD_genes.txt",sep="\t")
id.dd<-read.table("ID-DD_genes.txt",sep=",")
refseq<-read.table("refseq_genes.txt")

#Find overlap between NDD genes and apoptosis genes
length(intersect(apoptosis$V1,asd$V1))
length(intersect(apoptosis$V1,schiz$V1))
length(intersect(apoptosis$V1,id.dd$V1))

#Output overlapping genes to text files
asd_apoptosis_genes<-intersect(apoptosis$V1,asd$V1)
write.table(asd_apoptosis_genes,file="autism_apoptosis_genes_overlap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

schiz_apoptosis_genes<-intersect(apoptosis$V1,schiz$V1)
write.table(schzi_apoptosis_genes,file="schizophrenia_apoptosis_genes_overlap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

id.dd_apoptosis_genes<-intersect(apoptosis$V1,id.dd$V1)
write.table(id.dd_apoptosis_genes,file="ID-DD_apoptosis_genes_overlap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Pick random genes (based on all NDD genes) from all RefSeq genes
#Find number of genes that are in the apoptosis list
#Repeat simulation 100,000 times, and save in data frame
n2_asd<-length(asd$V1)
sim2_asd<-replicate(100000, {
	sim=sample(refseq$V1,n2_asd,replace=FALSE) 
	freq=sum(sim %in% apoptosis$V1) 
	})

n2_schiz<-length(schiz$V1)
sim2_schiz<-replicate(100000, {
	sim=sample(refseq$V1,n2_schiz,replace=FALSE) 
	freq=sum(sim %in% apoptosis$V1) 
	})

n2_id.dd<-length(id.dd$V1)
sim2_id.dd<-replicate(100000, {
	sim=sample(refseq$V1,n2_id.dd,replace=FALSE) 
	freq=sum(sim %in% apoptosis$V1) 
	})

#Get summary statistics, and output simulation results to text files
fivenum(sim2_asd)
mean(sim2_asd<length(intersect(apoptosis$V1,asd$V1))) #Find percentile score
write.table(sim2_asd,file="asd_in_apoptosis_100k_sim.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

fivenum(sim2_schiz)
mean(sim2_schiz<length(intersect(apoptosis$V1,schiz$V1))) #Find percentile score
write.table(sim2_schiz,file="schiz_in_apoptosis_100k_sim.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

fivenum(sim2_id.dd)
mean(sim2_id.dd<length(intersect(apoptosis$V1,id.dd$V1))) #Find percentile score
write.table(sim2_id.dd,file="ID-DD_in_apoptosis_100k_sim.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


#Generate histogram of data, with line added for the number of genes in the observed apoptosis/NDD datasets
hist(sim2_asd,main="Simulated autism genes in apoptosis genes",xlab="Overlap")
abline(v=length(intersect(apoptosis$V1,asd$V1)),col="red",lwd=3,lty=2)

hist(sim2_id.dd,main="Simulated ID/DD genes in apoptosis genes",xlab="Overlap")
abline(v=length(intersect(apoptosis$V1,id.dd$V1)),col="red",lwd=3,lty=2)

hist(sim2_schiz,main="Simulated schizophrenia genes in apoptosis genes",xlab="Overlap")
abline(v=length(intersect(apoptosis$V1,schiz$V1)),col="red",lwd=3,lty=2)