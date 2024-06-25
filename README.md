# Transportome Decision
**Overview:** This pipeline takes GBLAST results of a genome and determines all transport systems found in the genome.

## Getting Started
Ensure the following dependencies are installed to run ![GBLAST](https://github.com/SaierLaboratory/TCDBtools/blob/master/manuals/BioV_manual.pdf). In addition, the pfam has to be installed as well.

Download chebi ontologies to get substrate annotations for each transport systems. Download chebi.obo.gz: https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/?C=D;O=A

File named tcdb.faa should be present as a FASTA file containing the sequences of each TCDB system.

Additionally, genomic information such as the feature table as well gbff file containing orientation of the genome must be present


## Run Analysis
Run `python microbiome_main.py -g GBLAST FOLDER -f PATH TO FEATURE TABLE -b PATH TO GBFF FILE`

Some helpful flags:

- -q: query coverage threshold (0-100) (Default: 80)
- -s: subject coverage threshold (0-100) (Default: 80)
- -r: automatic coverage for rejection (Default: 20)
- -o: Minimum TMS overlap to be considered a good hit (Default: 50)
- -m: Membrane Protein Threshold
- -t: Minimum TMS regions to be classified a membrane protein (Default: 3)

Results of analysis will be found in analysis folder of the GBLAST results.

Output:

- Green.tsv: Transport systems that are found with high confidence.
- Yellow.tsv: Transport systems that require further analysis.
- Red.tsv: Transport systems that were not found in the genome.

## Create Master Table
Collect all the greens files from your analysis using this command:

`mkdir greens`

`cp */analysis/Green.tsv greens/`




