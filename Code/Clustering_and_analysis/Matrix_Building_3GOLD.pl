#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Matrix_Building_3GOLD;

my $TAG_LIST = ''; 
my $SUBTAG_LIST = ''; 
my $ERROR_THRESHOLD = ''; 
my $WEIGHT_THRESHOLD = ''; 
my $TOLERATED_FRAME_SHIFT = ''; 
my $INSERTION_WEIGHT = ''; 
my $DELETION_WEIGHT = ''; 
my $SUBSTITUTION_WEIGHT = '';
my $SINGLETONS_OUT = ''; 
my $AVE_DIST_OUT = ''; 
my $CLUSTERS_OUT = '';
my $SIZE_OUT = ''; 
my $QUIET = ''; 
my $THREADS = ''; 
my $usage = "\n\n$0 [options] \n
Required Options:
	-i_full		A text file of barcodes or UMIs in linearized fasta format. The full dataset. 
	-i_sub		A text file of barcodes or UMIS in linearized fasta format. The subset to run in parallel. 
	-error		Total number of errors allowed between related sequences 
	-weight		Cumulative weight of errors allowed between related sequences 
	-fs     	The number of prefix or suffix bases allowed to not match between two related sequences
	-ins		The weight assigned to an insertion error
	-del		The weight assigned to a deletion error
	-sub		The weight assigned to a substitution error
	-s_out		The output file containing singletons
	-dist_out	The output file containing average distance between the centroid and clustered sequences
	-clust_out  The output file containing clusters
	-size_out	The output file containing cluster size look up information 
Optional Options:
	-quiet		Silences progress reporting by assigning 'on' to the argument.
	-threads	The number of threads to use. Default is the number of cores available on your machine. 
				An input of '0' will result in default parameters. 
    -help		Show this message.\n";
GetOptions(
	'i_full=s' => \$TAG_LIST, 
	'i_sub=s' => \$SUBTAG_LIST, 
	'error=s' => \$ERROR_THRESHOLD, 
	'weight=s' => \$WEIGHT_THRESHOLD, 
	'fs=s' => \$TOLERATED_FRAME_SHIFT, 
	'ins=s' => \$INSERTION_WEIGHT, 
	'del=s' => \$DELETION_WEIGHT, 
	'sub=s' => \$SUBSTITUTION_WEIGHT,
	's_out=s' => \$SINGLETONS_OUT, 
	'dist_out=s' => \$AVE_DIST_OUT,
    'clust_out=s' => \$CLUSTERS_OUT, 
	'size_out=s' => \$SIZE_OUT, 
	'quiet=s' => \$QUIET,
	'threads=s' => \$THREADS,
	 help => sub{pod2usage($usage); },
) or pod2usage(2);
unless ($TAG_LIST){ die "\nProvide a text file of tags, barcodes, or umis. One per line in either fasta or list format. $usage"; } 
unless ($SUBTAG_LIST){ die "\nProvide a text file of tags, barcodes or umis. One per line in either fasta or list format. Must be a subdataset for parallelized computing. $usage";}
unless ($ERROR_THRESHOLD){ die "\nProvide the total number of errors allowed between related sequences. $usage"; } 
unless ($WEIGHT_THRESHOLD){ die "\nProvide the cumulative weight of errors allowed between related sequences. $usage"; } 
unless ($TOLERATED_FRAME_SHIFT){ die "\nProvide the number of prefix or suffix bases allowed to not match between two related sequences. $usage"; } 
unless ($INSERTION_WEIGHT){ $INSERTION_WEIGHT = 0;} 
unless ($DELETION_WEIGHT){ die "\nProvide the weight assigned to a deletion error. $usage"; } 
unless ($SUBSTITUTION_WEIGHT){ $SUBSTITUTION_WEIGHT = 0;} 
unless ($SINGLETONS_OUT){ die "\nProvide the output name for the file containing singletons. $usage"; } 
unless ($AVE_DIST_OUT){ die "\nProvide the output name for the file containing the average distance between the centroid and clustered sequences. $usage"; }
unless ($CLUSTERS_OUT){ die "\nProvide the output name for the file containing the clustered sequences. $usage"; }
unless ($SIZE_OUT){ die "\nProvide the output name for the file containing the cluster size lookup hash. $usage"; }

my $ThreeGold = Matrix_Building_3GOLD->new($ERROR_THRESHOLD, $WEIGHT_THRESHOLD, $TOLERATED_FRAME_SHIFT, $INSERTION_WEIGHT, $DELETION_WEIGHT, $SUBSTITUTION_WEIGHT, $SINGLETONS_OUT, $AVE_DIST_OUT, $CLUSTERS_OUT, $SIZE_OUT, $QUIET, $THREADS);

$ThreeGold->loadDataset($TAG_LIST, $SUBTAG_LIST);

$ThreeGold->createSequenceClusters();