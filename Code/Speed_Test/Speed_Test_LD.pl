#!/usr/bin/perl
use strict; 
use warnings;
use List::Util 'min';
use Getopt::Long;

# Input parameters 
my $DATASET = ''; 
my $usage = "\n\n$0 [options] \n
Options: 
    -i      Dataset input
    -help   Show this message.\n";
GetOptions( 
    'i=s' => \$DATASET,
    help => sub{pod2usage($usage); },
) or pod2usage(2);
unless($DATASET){die "\nProvide a dataset for processing. $usage";}
my $insertion_weight = 1; 
my $deletion_weight = 1; 
my $substitution_weight = 1; 

# Get the 'reference' sequence 
open (IN, '<', $DATASET) or die $!; 
my $reference; 
my $count = 0; 
while (my $line = <IN>){ 
    chomp $line; 
    $count ++; 
    $reference = $line; 
    if ($count > 0){ 
        last; 
    } 
}
close $DATASET; 

# Process edit distance calculations 
open(IN, '<', $DATASET) or die $!; 
while (my $sequence = <IN>){ 
    chomp $sequence; 
    my $dist_1 = _levdist($reference, $sequence); 
} 

sub _levdist {
	my ( $seq1, $seq2 ) = (@_)[0,1];

	# Do we have valid inputs?
	if (!defined($seq1) || $seq1 eq '') {
		print "Invalid ", $seq1, " input to _levdist, Empty or undefined!\n";
		return undef;
	}
	if (!defined($seq2) || $seq2 eq '') {
		print "Invalid ", $seq2, " input to _levdist, Empty or undefined!\n";
		return undef;
	}

	# Do we have matching inputs?
	if ($seq1 eq $seq2){ 
        return (0, 0); 
    }

	# Do our lengths match?
	if ( length($seq1) != length($seq2) )  {
		print "_levdist error: strings seq1 ,", $seq1, " and seq2 ,", $seq2, " are not the same size!\n";
		return undef;
	}
	
    my $l1 = length($seq1);
    my $l2 = length($seq2);
    my @s1 = split '', $seq1;
    my @s2 = split '', $seq2;
    my $distances;     
    for (my $i = 0; $i <= $l1; $i++) {
	$distances->[$i]->[0] = $i;
    }
    for (my $j = 0; $j <= $l2; $j++) {
	$distances->[0]->[$j] = $j;
    }
    for (my $i = 1; $i <= $l1; $i++) {
	    for (my $j = 1; $j <= $l2; $j++) {
	        my $cost; 
            if ( $s1[$i-1] eq $s2[$j-1] ) {
		        $cost = 0;
	        } else {
		        $cost = 1;
	        }
            $distances->[$i]->[$j] = min($distances->[$i-1]->[$j-1] + $cost, 
                                        $distances->[$i]->[$j-1]+1, 
                                        $distances->[$i-1]->[$j]+1 );
        }
    }
    my $classic_distance = $distances->[$l1]->[$l2];
    return $classic_distance;
}
