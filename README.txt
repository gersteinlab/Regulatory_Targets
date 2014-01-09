compile: make clean 
	 make
requirement: gsl, boost regex and alglib,please alter the 'makefile' if you didn't install the boost librar in the default directory.

Usage: ./DRM2Target <options>
	Version: 1.07
	-b [string:required]	Bed file for DRM regions. test file: test.bed
	-g [string:required]	GTF file for gencode_v7 tss based annotation ../input_data/opt_g.gencode.v7.tss.gff
	-r [string:required]	RPM file in gff format ../input_data/opt_r.gencode.v7.tss.rpm.tsv
	-s [string:required]	Correlation significance test file is: opt_s.bed.gz ../input_data/opt_s.bed.gz
	-d [string:required]	Directory for compressed reads alignment file for methylation, H3K4me1, and H3K27ac ../bgzips/
	-l [string:required]	Meta information file for aligned data file ../input_data/opt_l.drm.meta
	-m [int:10000]      	Minimum distance between DRM and its closest gene
	-M [int:1000000]    	Maximum distance between DRM and genes,
	-p [int:0.05]       	pvalue cutoff, should be not greater than 0.05, and if '-j' is set, it is used as the adjusted pvalue cutoff for any of four adjustment methods: Bonferroni, Holm, BH and BY. Tips: if you want to output all results, please set '-p' greater than 0.05
	-j                  	Whether p-value is adjusted based on pre-calculated data. Default: not adjusted.
	-J                  	To enable '-J' will disable '-j'.  Whether p-value is adjusted within user's data. Default: not adjusted.
	-N [int: 3000]      	The maximum number of DRM regions can be suggested to run with this program. This tool only works for a small number of regions in a batch
	-c [string:1,1,1]   	Correlation methods used for methylation, H3K4me1 and H3K27ac, string delimeted by ',' with 0 and 1.[0, Pearson; 1,Spearman]
	-t [string:0,1,1]   	Tailtype for correlation significance test for methylation, H3K4me1 and H3K27ac. [0,left-tail; 1,right-tail; 2,both-tails]. When using '-j' option, use default only. Other combinations are not available for adjustment with pre-calculated data. If you want to get the raw p-value or adjustment within user's data, please don't enable '-j' option! 
	-o [string:required]	output file
	-h                  	print this help


An example:./DRM2Target -b test.bed ../input_data/opt_g.gencode.v7.tss.gff -r ../data/opt_r.gencode.v7.tss.rpm.tsv -s ../input_data/opt_s.bed.gz -d ../bgzips/ -l ../input_data/opt_l.drm.meta -o test.out
