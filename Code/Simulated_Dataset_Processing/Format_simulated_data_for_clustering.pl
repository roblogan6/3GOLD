#!/usr/bin/perl 
use strict; 
use warnings; 
use List::Util qw(shuffle); 

# 4th in the series.

my $dir = $ARGV[0]; 
my @files = glob "$dir/*";

my $number_of_clusters = $ARGV[1]; 
my $number_of_reads = $ARGV[2]; 
my @output; 

foreach my $file ( (shuffle(@files)) ){ 
    if( $file =~ /(.*)([TAGC]{20})(.*sequences_to_cluster)(\._)(\d+)/){
        my $tag = $2; 
        my $sequence_count = $5; 
        if ($sequence_count >= $number_of_reads) { 
            push(@output, $file, "\n"); 
        } 
    }
    my $size_of_output = @output; 
    if ($size_of_output == ($number_of_clusters * 2)){ 
        last; 
    }
}

# now go through all the files in output, open them, and format them into: ref and to cluster files
my @ref_output; 
my @to_cluster_output; 
foreach my $filtered_file(@output){ 
    # Get centroid information
    if ($filtered_file =~ /(.*)([TAGC]{20})(.*sequences_to_cluster)(\._)(\d+)/){
        my $centroid = $2;
        push(@ref_output, "\n>", $centroid, "\n", $centroid, "\t"); 
        push(@to_cluster_output, "\n", ">", $centroid, "\n", $centroid);
        my $count = 1; 
        # Populate output files
        open(IN, '<', $filtered_file) or die $!; 
        while (my $line = <IN>){ 
            chomp $line; 
            if ($line =~ /^>/){ 
                next; 
            } elsif ($line =~ /([TAGC]{20})/){ 
                push(@to_cluster_output, "\n", ">", $1, "\n", $1); 
                push(@ref_output, $1, "\t");
                $count++; 
            }
            if ($count == ($number_of_reads)){ 
                last; 
            }
        }
    }
}
open(REF_OUT, '>', "/Users/roblogan/Desktop/3GOLD_Simulated/Sequences_To_Cluster/50x10_ref.txt") or die $!; 
open(TO_CLUSTER_OUT, '>', "/Users/roblogan/Desktop/3GOLD_Simulated/Sequences_To_Cluster/50x10_to_cluster.txt") or die $!; 
print REF_OUT @ref_output; 
print TO_CLUSTER_OUT @to_cluster_output; 