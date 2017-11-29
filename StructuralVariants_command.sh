
## DevRO command line options

Step1. run variantcallers for calling SVs 
 perl VariantCaller_dup.pl genomeFile.txt prefix_Output PATHBAMFiles.config
 perl VariantCaller_inv.pl genomeFile.txt prefix_Output PATHBAMFiles.config
 perl VC_dels_ins_refdel.pl regionFile.txt OutputFile.name
 
Output variantcallers, more details in DevRO/README_VariantCaller.txt
	1.bin
	2.coordinates 1k window scan
	3.size if raw SV
	4.coordinates of mates window
	5.median_mate_postion
	6.coordinates_raw_SV
	7.median position of singletons on forward strand
	8.median position of singletons on forward strand
	9.Hash seperated inverted reads counts in each population on forward strand
	10.Hash seperated inverted mate counts in each population on forward strand
	11.Hash seperated inverted reads counts in each population on reverse strand
	12.Hash seperated inverted mate counts in each population on reverse strand 
	13-16 Softclips same as 9-12
	17-20 Singletons same as 9-12
	21 Duplication reads counts in each population hash seperated
	22 Duplication mate counts in each population hash seperated
	23 Total reads counts in each population hash seperated
	24 Total mate counts in each population hash seperated

Step2. run parser to parse the variantcaller results and compare the domestic to wild populations
 perl Parser_2.2.pl Output_VariantCaller.txt OutputFile.name
 perl Parser_Reference_dels_1.pl Output_VC_dels_ins_refdel.txt OutputFile.name

Step3. Annotate the CNVs using normalized depth of coverage in 1kbp windows with M-value and p-value

### Breakdancer command line options
#step1. Generating the config files using breakdancer perl script:
~/breakdancer_perl/bam2cfg.pl -h -g Sample1.RG.merged.remdup.bam Sample2.RG.merged.remdup.bam Sample3.RG.merged.remdup.bam Sample4.RG.merged.remdup.bam >sampleAll.cfg

#step2. Running breakdancer cpp main script: Jobs_Breakdancer.txt
breakdancer-max -o chr1 -q 10 -d Fastq_SVReadsFile_1 -g GBrowse_OUTFile_1.bed -l -h sampleAll.cfg >chr1_SV_allSample.ctx
... on all chromosomes

#Output file format:
#Chr1	Pos1	Orientation1	Chr2	Pos2	Orientation2	Type	Size	Score	num_Reads	num_Reads_lib	Allele_frequency	ind1.bam	ind2.bam ind3.bam ind4.bam

### SVDetect command line options

perl BAM_preprocessingPairs.pl Sample.bam
SVDetect linking filtering -conf config_Sample.conf
SVDetect links2circos -conf config_Sample.conf
SVDetect links2bed -conf config_Sample.conf
SVDetect links2compare -conf config_Sample.conf
SVDetect links2SV -conf config_Sample.conf
SVDetect cnv -conf config_Sample.conf
 ## SVDetect is run on domestic and wild populations separately
### Parsing the outputs.
1. For Breakdancer and SVDetect, take only observations above score 80
2. Only extract the Structural variants with significant difference between domestics or wild populations.