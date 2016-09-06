path_w4CSeq <- "/var/www/html/" # path to the software
path_bwa <- "/var/www/html/w4cseq/bin/bwa-0.7.12/" # path to bwa
path_samtools <- "/var/www/html/w4cseq/bin/samtools-1.2/" # path to samtools
path_bedtools <- "/var/www/html/w4cseq/bin/bedtools2-2.25.0/bin/" # path to bedtools
path_RCircos <- "/var/www/html/w4cseq/bin/R-3.1.2/library" # path to RCircos R package
path_quantsmooth <- "/var/www/html/w4cseq/bin/R-3.1.2/library" # path to quantsmooth R package
#######################################################################

proc <- 1 # number of threads
exp_name <- "4C_exp1"
file_in <- "/var/www/html/w4cseq/work/860/query1.fq" # path/name to input fastq file
unzip <- "no" #whether your data is uncompressed, note that you cannot judge by the suffix like .fq or .fq.gz, since it's not reliable.
		#Nevertheless, generally, it should be "yes" when your file's suffix is .fq or .fastq
						 # and "no" when your file's suffix is .fq.gz or .fastq.gz

build <- "mm10" # reference genome name can be chosen from mm9/mm10/hg18/hg19
primer_frag <- "ATCTGCTATTGAGGAAGCTT" # primer sequence in bait region
enzyme <- "AAGCTT" # recoginition site of enzyme
bait_ch <- "chr17"
bait_st <- "35497036"
bait_en <- "35504156"

size_inter <- 500
size_intra <- 100
window_intra <- 3000
FDR <- 0.05


