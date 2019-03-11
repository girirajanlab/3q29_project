#Description

This GitHub repository contains the scripts and pipelines used to generate bioinformatic data related to the analysis of 3q29 gene homlogs in *Drosophila* models.

There are two directories in this repository: 

1. Pipelines for identifying differentially-expressed genes from RNA-Sequencing data of 3q29 homologs. 
..* This directory contains a batch script for aligning and quantifying read counts using [TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml) and [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/), and an R pipeline for identifying differentially-expressed genes using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
..* The raw RNA-Seq reads and quantified read counts for biological replicates are available at [NCBI GEO accession number GSE128094](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE128094).

2. Simulation of apoptosis gene enrichment among schizophrenia genes.
..* This directory contains an R script which simulates enrichment of apoptosis genes in schizophrenia gene sets, and two raw text files containing apoptosis and RefSeq gene sets used in the analysis.

# Citation
Singh MD, Jensen M, Lasser M, Huber E, Yusuff T, Pizzo L, Lifschutz B, Desai I, Kubina A, Yennawar S, Kim S, Iyer J, Rincon-Limas DE, Lowery LA, Girirajan S. *NCBP2*-mediated interactions contribute to neurodevelopmental defects of the schizophrenia-associated 3q29 deletion. 

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