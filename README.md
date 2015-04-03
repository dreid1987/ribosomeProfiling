One-step ribosome profiling analysis
----------------------

The objective of this software is to create a single tool that can execute analysis of ribosome profiling data with a single step, without the need to seek out the several programs that are typically required for this analysis.

INPUT:

1. fastq files from deep sequencing analysis

2. A gtf (gene annotation) file that includes start and stop codons

3. A set of fa (fasta format) genome sequence data

OUTPUT:

1. totalExpression.csv: a tab-separated spreadsheet that indicates the level of expression for each gene (expressed as RPKM). This may be used for downstream statistical analysis.

2. individGenes: graphs of the distribution of ribosomes on each mRNA

3. riboPosition: metagene plots that indicate the positions of ribosomes relative to the start codon across the transcriptome


QUICK START
----------------------
1. Create a work folder

	Make a new folder where you will perform ribosome profiling analysis. Extract all of the files included in this package into this new folder

2. Downloading reference sequences

	This analysis requires a reference genome to map reads to. These can be downloaded in whole from Illumina iGenomes: http://ccb.jhu.edu/software/tophat/igenomes.shtml

	After downloading a reference genome, extract the entire compressed package into a folder called 'database' in your work folder.

	Most of the Refseq iGenomes include all of the required components to successfully complete this analysis. These are:
		1. a gtf file that indicates exon positions, gene names, start codons, and stop codons
		2. a set of fa files (fasta format) that indicate the sequence of each chromosome, all in one file.

3. Extract deep sequencing reads into the 'reads' folder

	Most deep sequencing reads are downloaded as compressed .fastq files. Use any decompression utility (eg 7zip) to extract these files into a folder called 'reads' that is in your work folder

	Each of the different .fastq files will be analyzed and quantified separately
	
	You may use both ribosome profiling and mRNA-seq data as input. Each will be quantified identically
	
4. Execute analysis

	Open Terminal and navigate to your work folder. 
	
	Type in 'python profiling.py'
	
	After the program starts, indicate the name of your reference genome and the primers that you used in library preparation. The analysis will then execute without the need for further input
	
	This analysis may take several hours to execute, depending on your reference genome and the size and number of your deep sequencing libraries


OUTLINE OF APPROACH
----------------------
These are the steps that this program automatically executes:

1. Installation of all required software
		Tests whether the required software is installed, and if not, then installs it.
2. Adapter removal from deep sequencing reads
		Uses Cutadapt
							https://cutadapt.readthedocs.org/en/stable/
							http://journal.embnet.org/index.php/embnetjournal/article/view/200/479
3. Preparation of reference transcriptome
		Uses Tophat and Cufflinks to identify isoform abundances, then selects the most abundant isoform for analysis using chooseTranscripts.py
							http://ccb.jhu.edu/software/tophat/index.shtml
							http://cole-trapnell-lab.github.io/cufflinks/
4. Mapping of deep sequencing reads
		Uses Bowtie with mostly default settings, allowing for no mismatches
							http://bowtie-bio.sourceforge.net/
5. Assessment of ribosome profiling mapped reads and output generation
		Uses diff.py

CUSTOMIZATION
----------------------
There are several variables in the analysis that the user may wish to change. For each of these, open 'diff.py' in a text editor and change the appropriate number at the top of the file.

readCutoff: the number of reads that must map to a gene in order for it to be considered for analysis (default=3)

individReadCutoff: the number of reads that must map to a gene for it to be plotted in the individGenes output (default=2500)

posxmin: the minimum position on the metagene plot (in the ribosomePos folder) relative to the start codon (default=-20)

posxmax: the maximum position on the metagene plot (in the ribosomePos folder) relative to the start codon (default=150)

ribo5Add: number of nucleotides to add to the start of the read when defining the position of the ribosome (default=14)

TROUBLESHOOTING
----------------------
Some of the Illumina iGenomes lack all of the annotation required - for example, Saccharomyces cerevisiae. In this case, the user must generate a gtf file with a format resembling the following:

1	unknown	exon	196788861	196789032	.	+	.	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	CDS	196788975	196789032	.	+	0	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	start_codon	196788975	196788977	.	+	.	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	CDS	196794607	196794801	.	+	2	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	exon	196794607	196794801	.	+	.	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	CDS	196800927	196801126	.	+	2	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	exon	196800927	196801319	.	+	.	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";
1	unknown	stop_codon	196801127	196801129	.	+	.	gene_id "CFHR1"; gene_name "CFHR1"; p_id "P234"; transcript_id "NM_002113.2"; tss_id "TSS23824";

The user will then enter the location of this gtf file after running profile.py
