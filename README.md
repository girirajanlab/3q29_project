# Description

This GitHub repository contains the scripts and pipelines used to generate bioinformatic data related to the analysis of *Drosophila* homologs of 3q29 genes.

There are three directories in this repository: 

1. Pipelines for identifying differentially-expressed fly genes from RNA-Sequencing data for homologs of 3q29 genes. 
  * This directory contains a batch script for aligning and quantifying read counts using [TopHat2 v.2.1.1](https://ccb.jhu.edu/software/tophat/index.shtml) and [HTSeq v.0.6.1](https://htseq.readthedocs.io/en/release_0.11.1/), and an R pipeline for identifying differentially-expressed genes using [edgeR v.3.20.1](https://bioconductor.org/packages/release/bioc/html/edgeR.html).  
  * The raw RNA-Seq reads and quantified read counts for biological replicates are available at [NCBI GEO accession number GSE128094](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE128094).

2. Simulation of apoptosis gene enrichment among candidate neurodevelopmental genes.
 * This directory contains an R script which simulates enrichment of apoptosis genes in schizophrenia gene sets, and raw text files containing apoptosis, RefSeq, and candidate NDD gene sets used in the analysis.
 * The candidate gene sets are derived from [Purcell et al, *Nature* 2014](https://www.ncbi.nlm.nih.gov/pubmed/24463508) (schizophrenia), [SFARI Gene](https://gene.sfari.org/) (autism), and [Developmental Delay G2P database](https://www.ebi.ac.uk/gene2phenotype) (ID/DD).

3. Network analysis of CNV genes and simulated random gene sets.
 * This directory contains two Python scripts for analyzing the connectivity of genes within a human brain-specific interaction network. 
 * One script (`nearest_neighbor_weighted_allgenes.py`) takes an individual gene as standard input (gene name and Entrez ID), and outputs the shortest distance to every other gene in the network, as well as the genes located in the shortest path between the two target genes. For example, to find the connectivity of NCBP2 (Entrez ID 22916) to all other genes in the genome, one would run:
 `python nearest_neighbor_weighted_allgenes.py NCBP2_22916`
 * The other script (`nearest_neighbor_weighted_allgenes_sim.py`) calculates the connectivity between groups of random genes in the genome. It picks 50 random start genes from a set of RefSeq genes in the network, and then detrmines the connectivity of each of the target genes to 100 sets of 20 randomly-selected target genes. The output file lists each of the 50 start genes and the average shortest distance to the sets of target genes. The script takes a replicate number as input, to allow for multiple replciates to be run simultaneously in the same directory:
 `python nearest_neighbor_weighted_allgenes_sim.py 1`
 * The network file used in this analysis (`brain.degnorm-ge2.prob-gept02.dat`) is described in [Greene *et al, Nat. Genet.* 2015](https://www.ncbi.nlm.nih.gov/pubmed/25915600) and [Krishnan *et al, Nat. Neurosci.* 2016](https://www.ncbi.nlm.nih.gov/pubmed/27479844); we generated a sub-network that only contained edges with weights >2.0 (the top 0.5% of interactions in the network).

Bash pipeline scripts can be run in any Unix environment. R scripts can be run using any R version (scripts were generated using R v.3.4.2.). Python scripts for network analysis can be run in Python2 (scripts were generated using Python v.2.7.16) and require the [NetworkX package v.2.4](https://networkx.github.io/). 

# Citation
[Singh MD, Jensen M, Lasser M, Huber E, Yusuff T, Pizzo L, Lifschutz B, Desai I, Kubina A, Yennawar S, Kim S, Iyer J, Rincon-Limas DE, Lowery LA, Girirajan S. *NCBP2* modulates neurodevelopmental defects of the 3q29 deletion in *Drosophila* and *X. laevis* models. **BioRxiv** 614750; 16 Oct. 2019.](https://www.biorxiv.org/content/10.1101/614750v2)

# Copyright/License
The code in this repository is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <https://www.gnu.org/licenses/>.

# Contact
For questions or comments, please contact Matthew Jensen (mpj5142@psu.edu) or Santhosh Girirajan (sxg47@psu.edu).