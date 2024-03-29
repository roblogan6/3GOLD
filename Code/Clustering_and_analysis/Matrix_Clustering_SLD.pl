#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Matrix_Clustering_SLD;

my $TAG_LIST = '';
my $SINGLETONS_IN = ''; 
my $AVE_DIST_IN = ''; 
my $CLUSTERS_IN = ''; 
my $SIZE_IN = ''; 
my $ERROR_THRESHOLD = '';
my $TAG_OUTPUT_SEQ = ''; 
my $TAG_OUTPUT_ID = ''; 
my $QUIET = ''; 
my $THREADS = ''; 
my $usage = "\n\n$0 [options] \n
Required Options:
    -i          A text file of barcodes or UMIs in linearized fasta format.
	-s  		The text file of singletons 
    -dist       The text file of distances 
    -clust      The text file of clusters
    -size       The text file of sizes for lookup 
	-error		Total number of errors allowed between related sequences
	-o_seq		The output file containing clusters of barcodes or UMIs
	-o_id		The output file containing clusters of fasta headers
Optional Options:
	-quiet		Silences progress reporting by assigning 'on' to the argument.
	-threads	The number of threads to use. Default is the number of cores available on your machine. 
				An input of '0' will result in default parameters. 
    -help		Show this message.\n";
GetOptions(
    'i=s' => \$TAG_LIST, 
	's=s' => \$SINGLETONS_IN,
    'dist=s' => \$AVE_DIST_IN, 
    'clust=s' => \$CLUSTERS_IN,  
    'size=s' => \$SIZE_IN, 
	'error=s' => \$ERROR_THRESHOLD,
    'o_seq=s' => \$TAG_OUTPUT_SEQ, 
    'o_id=s' => \$TAG_OUTPUT_ID, 
	'quiet=s' => \$QUIET,
	'threads=s' => \$THREADS,
	 help => sub{pod2usage($usage); },
) or pod2usage(2);
unless ($TAG_LIST){ die "\nProvide a text file of tags, barcodes, or umis. One per line in either fasta or list format. $usage"; } 
unless ($SINGLETONS_IN){ die "\nProvide a text file of singletons, produced by matrix building. $usage"; }
unless ($AVE_DIST_IN){ die "\nProvide a text file of cluster distances, produced by matrix building. $usage"; }
unless ($SIZE_IN){ die "\nProvide a text file of the hoa size for lookup, produced by matrix building. $usage"; }
unless ($CLUSTERS_IN){ die "\nProvide a text file of clusters, produced by matrix building. $usage"; }
unless ($ERROR_THRESHOLD){ die "\nProvide the total number of errors allowed between related sequences. $usage"; } 
unless ($TAG_OUTPUT_SEQ){ die "\nProvide the output name for the file containing sequences. $usage"; } 
unless ($TAG_OUTPUT_ID){ die "\nProvide the output name for the file containing header information. $usage"; }

my $ThreeGold_MacOS_Matrix_Clustering = Matrix_Clustering_SLD->new($ERROR_THRESHOLD, $TAG_OUTPUT_SEQ, $TAG_OUTPUT_ID, $SINGLETONS_IN, $AVE_DIST_IN, $CLUSTERS_IN, $SIZE_IN, $QUIET, $THREADS);

$ThreeGold_MacOS_Matrix_Clustering->loadDataset($TAG_LIST);

$ThreeGold_MacOS_Matrix_Clustering->mergeSequenceClusters();

$ThreeGold_MacOS_Matrix_Clustering->locateSequenceClusterCentroids();

$ThreeGold_MacOS_Matrix_Clustering->mergeSingletons();

$ThreeGold_MacOS_Matrix_Clustering->generateResultSet();
