#This R script randomly selects 2,546 RefSeq genes (matching the number of schizophrenia candidate genes from
#Purcell et al, Nature 2014) and determines the number of genes in each random set that are annotated for 
#apoptosis genes (apoptotic process, GO:0006915 or children terms). A list of all apoptosis genes is provided in
#the file human_apoptotic_process_raw.txt, and a list of all RefSeq genes is provided in the file refseq_genes.txt.

apoptosis<-read.table("human_apoptotic_process_raw.txt",sep="\t")
refseq<-read.table("refseq_genes.txt")

#Pick random genes (based on all schizophrenia genes) from all RefSeq genes
#Find number of genes that are in the apoptosis list
#Repeat simulation 100,000 times and save in data frame
n2<-2546
sim2<-replicate(100000, {
	sim=sample(refseq$V1,n2,replace=FALSE) 
	freq=sum(sim %in% apoptosis$V1) 
	})

fivenum(sim2)
mean(sim2<268) #Find percentile score
write.table(sim2,file="schizophrenia_in_apoptosis_100k_sim.txt",quote=FALSE,row.names=FALSE,col.names=FALSE) #Output simulation results

#Generate histogram of data, with line added for the 262 genes in the observed apoptosis/schizophrneia dataset
hist(sim2,main="Simaulated schizophrenia genes in apoptosis genes",xlab="Overlap")
abline(v=262,col="red",lwd=3,lty=2)