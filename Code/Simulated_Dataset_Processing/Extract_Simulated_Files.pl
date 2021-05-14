#!/usr/bin/perl 
use warnings; 
use strict; 

# 1st script in the series. Since PaSS makes many reads of the input, sometimes many simulations of the same input sequence, 
# You need to fish out a simulated sequence to serve as a representative. That is what this script does. It reduces all simulated 
# reads to one per input sequence. 

# put the whole simulated data into memory 
my $simulated_file = '/Users/roblogan/PaSS/50_Centroids_Out.fasta'; 
open my $fh, '<', $simulated_file or die; 
$/ = undef; 
my $simulated_data = <$fh>; 
close $fh; 


my $dir = $ARGV[0]; 
my @ArrayofFiles = glob "$dir/*"; 

#opendir(DH, $ARGV[0]); 
#my @files = readdir(DH); 
#closedir(DH); 

foreach my $file (@ArrayofFiles){ 
    my $centroid;
    if ($file =~ m/(.*)([TAGC]{20})(\.txt)$/){  
        $centroid = $2; 
    } 
    #print $centroid, "\n"; 
    if ($simulated_data =~ /(>m.*)(\n)([TAGC]*)($centroid)([TAGC]*)/){ 
        #print "fasta: \n"; 
        #print $1.$2.$3.$4.$5, "\n"; 
        my @output; 
        push(@output, $1.$centroid.$2.$3.$4.$5); 
        open(OUT, '>', "$ARGV[1]/$centroid.txt") or die $!; 
        print OUT @output; 
    }
}