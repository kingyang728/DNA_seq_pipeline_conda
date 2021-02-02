# Control Bam file folder

This is a instruction for control bam folder. Control bam file which required for CNV calling should be downloaded in this folder, 'genesHG38.bed' is the gene interval bed file for genome hg38.


## Authors

* Jingyu Yang


## Usage


#### Prepare for the control bam file.

If there is no corresponding control bam, you can download the public control cram file [NA12878](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/exome_alignment/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram) 
And run the 'CramToBam.sh' script to generate the control bam file.
Start terminal at this directory and run the command below:
	
	sh CramToBam.sh


