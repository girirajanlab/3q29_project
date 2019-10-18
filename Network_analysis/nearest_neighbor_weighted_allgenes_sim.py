#This script calculates the connectivity between groups of random genes in the genome. It picks 50 random start genes from a set of RefSeq genes in the network, 
#and then detrmines the connectivity of each of the target genes to 100 sets of 20 randomly-selected target genes. The output file lists each of the 50 start genes 
#and the average shortest distance to the sets of target genes. The script takes a replicate number as input, to allow for multiple replciates to be run simultaneously in the same directory:
# python nearest_neighbor_weighted_allgenes_sim.py 1

#Import modules
import csv
import sys
import networkx
import itertools
import numpy as np
import statistics

#Open network file (fromat: Entrez_ID_1, Entrez_ID_2, Weighted_connectivity) and list of genes in network (Entrez IDs in column 1, gene names in column 3)
network_file=open("brain.degnorm-ge2.prob-gept02.dat",'rU')
network_lines=network_file.readlines()
gene_file=open("brain.genes.protein-coding.txt",'rU')
gene_lines=gene_file.readlines()

#Populate list of protein-coding genes in the network, and dictionary for Entrez-gene symbol conversion
coding_genes_list=[]
entrez_dict={}
for gene_line in gene_lines:
	gene_line=gene_line.strip().split("\t")
	gene_ID=gene_line[0]
	gene_symbol=gene_line[2]
	coding_genes_list.append(gene_ID)
	entrez_dict[gene_ID]=gene_symbol
coding_genes=set(coding_genes_list)

#Populate network data for each gene with NetworkX
network_graph=networkx.Graph()
for network_line in network_lines:
	network_line=network_line.strip().replace('\x00','').split("\t")
	node1=network_line[0]
	node2=network_line[1]
	weight=1/float(network_line[2]) #Take inverse of weight for edge lengths

	#Add edge to NetworkX network if both genes are in protein-coding list	
	if node1 in coding_genes:
		if node2 in coding_genes:
			network_graph.add_edge(node1,node2,weight=weight)

#Initialize output file, using replicate number from system input
rep_number=sys.argv[1]
outwriter = csv.writer(file("weighted_connectivity_sim_20_rep"+rep_number+".txt",'w'),delimiter='\t')
header=["Target gene","Average distance","Notes"]
outwriter.writerow(header)

i=0 #Initialize start gene counter 

#Randomly select 50 starting genes from the list of genes in the network
gene_start_targets=np.random.choice(list(coding_genes),size=50,replace=False)
for gene_start_target in gene_start_targets:
	i=i+1 #Increment start gene counter
	j=0 #Initialize target gene counter

	if network_graph.has_node(gene_start_target)==False: #Skip start gene if not in network
		outrow=[gene_start_target,'',"not_in_network"]
		outwriter.writerow(outrow)
		continue		

	#Randomly select 2000 target genes in network (20 other genes in CNV region * 100 replicates)
	coding_genes_sim=coding_genes
	coding_genes_sim.remove(gene_start_target) #Prevent gene from targeting itself
	gene_targets=np.random.choice(list(coding_genes_sim),size=2000,replace=True)
	shortest_path_lengths=[] #Initialize list of shortest paths from start gene to target gene

	#Iterate through target genes
	for gene_end_target in gene_targets:
		#Increment target gene counter and print counters
		j=j+1 
		print str(i)+"\t"+str(j)

		#Use NetworkX to generate the shortest paths between start and target genes, using Dijkstra's algorithm based on the weighted edges.
		#Append the length of the shortest path (based on edge weights) to master list.
		try:
			shortest_path_lengths.append(networkx.dijkstra_path_length(network_graph,source=gene_start_target,target=gene_end_target))
		except(networkx.exception.NetworkXNoPath): #Keep field empty if no path exists between the two genes.
			continue

	#Calculate average shortest path lengths
    try:
        avg_shortest_path_length=statistics.mean(shortest_path_lengths)
    except(statistics.StatisticsError):
        outrow=[gene_start_target_new,"No_connected_genes"]
        outwriter.writerow(outrow)
        continue

	#Convert start gene Entrez ID to gene symbol
	gene_start_target_new=entrez_dict[gene_start_target]

	#Output to file
	outrow=[gene_start_target_new,avg_shortest_path_length]
	outwriter.writerow(outrow)