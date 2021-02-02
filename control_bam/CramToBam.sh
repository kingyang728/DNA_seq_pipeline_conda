cd ..
samtools view -@ 2 -b  -T ./Ref_data/GRCh38.d1.vd1.fa -o ./control_bam/NA12878.exome.bam ./control_bam/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram
samtools index ./control_bam/NA12878.exome.bam
