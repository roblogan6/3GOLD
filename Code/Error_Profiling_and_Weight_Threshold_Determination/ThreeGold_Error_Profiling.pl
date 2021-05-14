#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes qw(gettimeofday);
use Time::Date;
use Math::Round; 
use Scalar::Util qw(looks_like_number); 
use List::Util qw( min );
use List::Util qw( max ); 
use POSIX; 
 
my $REF_SEQ = ''; 
my $MINUS = ''; 
my $SEQUENCE_LIST = ''; 
my $START = ''; 
my $OUTPUT = ''; 
my $QUIET = '';
my $TOLERATED_FRAME_SHIFT = ''; 
my $usage = "\n\n$0 [options] \n
Required Options:
	-ref		A text file of a reference sequence in linearized fasta format 
    -minus      Type 'Y' if the aligned BLAST fasta output is in the minus direction. Type 'N' if it is in the plus direction. 
	-i          A text file with sequences to compare to the reference in linearized multifasta format
                The input sequences need to all be the same length as the reference sequence
    -start      The depth of the starting position, indicated by the aligned 'query' start position. 
	-o          The address to print out the error profile report to 
Optional Options:
	-quiet		Silences progress reporting by assigning 'on' to the argument
    -fs     	The number of prefix or suffix bases allowed to not match between strings. Default is 15% of string length. 
    -help		Show this message.\n";
GetOptions(
	'ref=s' => \$REF_SEQ, 
    'minus=s' => \$MINUS, 
	'i=s' => \$SEQUENCE_LIST, 
    'start=s' => \$START, 
    'o=s' => \$OUTPUT,
	'quiet=s' => \$QUIET,
    'fs=s' => \$TOLERATED_FRAME_SHIFT, 
	 help => sub{pod2usage($usage); },
) or pod2usage(2);
unless ($REF_SEQ){ die "\nProvide a text file of a reference sequence in linearized fasta format. $usage"; } 
unless ($MINUS){ die "\nProvide the direction of the reference sequence. $usage"; }
unless ($SEQUENCE_LIST){ die "\nProvide a text file with sequences to compare to the reference in linearized multifasta format. The input sequences need to all be the same length as the reference sequence$usage"; } 
unless ($START){ die "\nProvide the aligned query's start position. $usage"; }
unless ($OUTPUT){ die "\nProvide an address to print out the error profile report to. $usage"; } 

if ($MINUS ne "Y"){ 
    if ($MINUS ne "N"){ 
        die "\nPlease assign either exactly 'Y' or 'N' for the 'minus' value. $usage"; 
    }
}

print "Does a single substitution error, compared to both an insertion and deletion error, occur less frequently (higher weight), more frequently (lower weight) or the same? Answer exactly 'less', 'more' or 'same': ";
my $SUB = <STDIN>; 
chomp $SUB; 
unless ($SUB eq 'less' || $SUB eq 'more' || $SUB eq 'same'){ die "\n Either answer with an exact 'less', 'more' or 'same';"}
print "Does an insertion error, compared to a deletion error, occur less frequently (higher weight), more frequently (lower weight) or the same? Answer exactly 'less' or 'more' or 'same': "; 
my $INS = <STDIN>; 
chomp $INS; 
unless ($INS eq 'less' || $INS eq 'more' || $INS eq 'same'){ die "\n Either answer with an exact 'less', 'more', or 'same';"}

sub _writeOut {
	my ($message) = @_;
	return if ( $QUIET eq "on");
    my $s = Time::Date->now; 
    print "Gregorian Time: ", $s, "\t", "Unix Time: "; 
    my ($seconds, $microseconds) = gettimeofday;
    printf("[%.3f]",$seconds + ($microseconds / 1000000));
	print "\tStatus: ".$message."\n";
}

loadDataset($REF_SEQ, $SEQUENCE_LIST); 
sub loadDataset {
	my ($REF_SEQ, $SEQUENCE_LIST) = @_;
    _writeOut("Loading sequences");
    my $ref;
    my $trimmed_seq; 
    my @sequences; 

	# Open the dataset
	open (SEQLIST, '<', $SEQUENCE_LIST) or die $!; 
    # Parse the dataset
	my $current_header;
	while (my $seq = <SEQLIST>){ 
        chomp $seq; 
        if ($seq =~ /(^>)(.*)/){
		    next;    
	    } elsif ($seq !~ /^>/){ 
		    # Remove leading and trailing blank spaces
		    $seq =~ s/^\s+|\s+$//g;
            # Skip blank tags and # prefixed
		    next if ( $seq eq '' || $seq =~ qr/^#/ );
            # Trim the sequence 
            my $n = ($START - 1); 
            $trimmed_seq = substr($seq, $n); 
            # Create an array
            #push(@sequences, $seq_2);
	    }  
    }
	close(SEQLIST);

    # Open the ref seq
    open(REF_SEQ, '<', $REF_SEQ) or die $!;  
    # Parse the dataset
	while (my $line = <REF_SEQ>){ 
        chomp $line; 
        if ($line =~ /(^>)(.*)/){
		    next;   
	    } elsif ($line !~ /^>/){ 
		    # Remove leading and trailing blank spaces
		    $line =~ s/^\s+|\s+$//g;
            # Skip blank tags and # prefixed
		    next if ( $line eq '' || $line =~ qr/^#/ );
            $ref = $line; 
            if ($MINUS eq "Y"){ 
                $ref = reverse $ref; 
                $ref =~ tr/ATGCatgc/TACGtacg/;
            }
	    }  
    }
	close(REF_SEQ);

    my $ref_length_test = length($ref);
    my $seq_length_test = length($trimmed_seq); 
    my $largest = max($ref_length_test, $seq_length_test);
    my $smallest = min($ref_length_test, $seq_length_test); 
    my $difference = ($largest - $smallest); 
    my $new_sequence; 


    if ($ref_length_test == $largest){ 
        my @output; 
        $new_sequence = substr $ref, 0, -$difference; 
        my $error_threshold; 
        my $var = looks_like_number($TOLERATED_FRAME_SHIFT);
        my @bases = split '', $new_sequence; 
        my $new_length = @bases;
        if ($var == 0){
            $TOLERATED_FRAME_SHIFT = round($new_length * 0.15); 
            $error_threshold = ($new_length - $TOLERATED_FRAME_SHIFT);
        } elsif ($var == 1) {
            $error_threshold = ($new_length - $TOLERATED_FRAME_SHIFT);
        }
        my $frameshift = $TOLERATED_FRAME_SHIFT; 
        _writeOut("Determining Error Profile");
        my @error_profile = _levdist($new_sequence, $trimmed_seq, $error_threshold, $frameshift);
        if (@error_profile){ 
            my $error_rate_percent = (100 * ($error_profile[6]/$new_length)); 
            my $rounded_error_rate = sprintf("%.2f", $error_rate_percent);
            push(@output, "Length (bases analyzed):\t", $new_length, "\t", "Total_errors:\t\t", $error_profile[6], "\t", "Error rate (percentage):\t", $rounded_error_rate, "\t", "Frameshift:\t", $frameshift, "\n"); 
            push(@output, "Number of insertions:\t\t", $error_profile[0], "\t", "Number of deletions:\t", $error_profile[2], "\t", "Number of substitutions:\t", $error_profile[4], "\n"); 
            my $percent_ins = (100* ($error_profile[0]/$error_profile[6])); 
            my $percent_del = (100 * ($error_profile[2]/$error_profile[6])); 
            my $percent_sub = (100 * ($error_profile[4]/$error_profile[6])); 
            my $highest_percentage = max($percent_ins, $percent_del, $percent_sub); 
            my $ins_weight; 
            my $rounded_ins_weight; 
            my $del_weight; 
            my $rounded_del_weight; 
            my $sub_weight; 
            my $rounded_sub_weight; 
            if ($highest_percentage eq $percent_ins){ 
                $ins_weight = 1;
                if ($percent_del != 0){ 
                    $del_weight = ($percent_ins/$percent_del); 
                    $rounded_del_weight = sprintf("%.2f", $del_weight); 
                }
                if ($percent_sub != 0){ 
                    $sub_weight = ($percent_ins/$percent_sub);
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight);
                }
                if ($percent_del == 0){ 
                    $del_weight = $ins_weight + $sub_weight + 1; 
                    $rounded_del_weight = sprintf("%.2f", $del_weight);
                }
                if ($percent_sub == 0 ){ 
                    $sub_weight = $ins_weight + $del_weight + 1; 
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight);
                }
                push(@output, "Insertion weight:\t\t", $ins_weight, "\tDeletion weight:\t", $rounded_del_weight, "\tSubstitution weight:\t\t", $rounded_sub_weight, "\n"); 
            }elsif ($highest_percentage eq $percent_del){ 
                $del_weight = 1; 
                if ($percent_ins != 0){ 
                    $ins_weight = ($percent_del/$percent_ins);
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight);  
                }
                if ($percent_sub != 0){ 
                    $sub_weight = ($percent_del/$percent_sub);
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight);
                }
                if ($percent_ins == 0){ 
                    $ins_weight = $del_weight + $sub_weight + 1; 
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight);
                }
                if ($percent_sub == 0 ){ 
                    $sub_weight = $del_weight + $ins_weight + 1;
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight); 
                }
                push(@output, "Insertion weight:\t\t", $rounded_ins_weight, "\tDeletion weight:\t", $del_weight, "\tSubstitution weight:\t\t", $rounded_sub_weight, "\n"); 
            }elsif ($highest_percentage eq $percent_sub){ 
                $sub_weight = 1; 
                if ($percent_ins != 0){ 
                    $ins_weight = ($percent_sub/$percent_ins);
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight); 
                }
                if ($percent_del != 0){ 
                    $del_weight = ($percent_sub/$percent_del);
                    $rounded_del_weight = sprintf("%.2f", $del_weight);
                }
                if ($percent_ins == 0){ 
                    $ins_weight = $del_weight + $sub_weight + 1; 
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight);
                }
                if ($percent_del == 0 ){ 
                    $del_weight = $sub_weight + $ins_weight + 1; 
                    $rounded_del_weight = sprintf("%.2f", $del_weight);
                }
                push(@output, "Insertion weight:\t\t", $rounded_ins_weight, "\tDeletion weight:\t", $rounded_del_weight, "\tSubstitution weight:\t\t", $sub_weight, "\n"); 
            }

        }
        push(@output, "Experimental sequence:\n", $trimmed_seq, "\n"); 
        push(@output, "Reference sequence:\n", $new_sequence, "\n"); 

        _writeOut("Printing Output");
        open(OUT, '>', $OUTPUT) or die $!; 
        print OUT @output; 
        _writeOut("Finished");
    }
    if ($seq_length_test == $largest){ 
        my @output; 
        $new_sequence = substr $trimmed_seq, 0, -$difference; 
        my $error_threshold; 
        my $var = looks_like_number($TOLERATED_FRAME_SHIFT);
        my @bases = split '', $new_sequence; 
        my $new_length = @bases;
        if ($var == 0){
            $TOLERATED_FRAME_SHIFT = round($new_length * 0.15); 
            $error_threshold = ($new_length - $TOLERATED_FRAME_SHIFT);
        } elsif ($var == 1) {
            $error_threshold = ($new_length - $TOLERATED_FRAME_SHIFT);
        }
        my $frameshift = $TOLERATED_FRAME_SHIFT; 
        _writeOut("Determining Error Profile");
        my @error_profile = _levdist($ref, $new_sequence, $error_threshold, $frameshift);
        if (@error_profile){ 
            my $error_rate_percent = (100 * ($error_profile[6]/$new_length)); 
            my $rounded_error_rate = sprintf("%.2f", $error_rate_percent);
            push(@output, "Length (bases analyzed):\t", $new_length, "\t", "Total_errors:\t\t", $error_profile[6], "\t", "Error rate (percentage):\t", $rounded_error_rate, "\t", "Frameshift:\t", $frameshift, "\n"); 
            push(@output, "Number of insertions:\t\t", $error_profile[0], "\t", "Number of deletions:\t", $error_profile[2], "\t", "Number of substitutions:\t", $error_profile[4], "\n"); 
            my $percent_ins = (100* ($error_profile[0]/$error_profile[6])); 
            my $percent_del = (100 * ($error_profile[2]/$error_profile[6])); 
            my $percent_sub = (100 * ($error_profile[4]/$error_profile[6])); 
            my $highest_percentage = max($percent_ins, $percent_del, $percent_sub); 
            my $ins_weight; 
            my $rounded_ins_weight; 
            my $del_weight; 
            my $rounded_del_weight; 
            my $sub_weight; 
            my $rounded_sub_weight; 
            if ($highest_percentage eq $percent_ins){ 
                $ins_weight = 1;
                if ($percent_del != 0){ 
                    $del_weight = ($percent_ins/$percent_del); 
                    $rounded_del_weight = sprintf("%.2f", $del_weight); 
                }
                if ($percent_sub != 0){ 
                    $sub_weight = ($percent_ins/$percent_sub);
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight);
                }
                if ($percent_del == 0){ 
                    $del_weight = $ins_weight + $sub_weight + 1; 
                    $rounded_del_weight = sprintf("%.2f", $del_weight);
                }
                if ($percent_sub == 0 ){ 
                    $sub_weight = $ins_weight + $del_weight + 1; 
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight);
                }
                push(@output, "Insertion weight:\t\t", $ins_weight, "\tDeletion weight:\t", $rounded_del_weight, "\tSubstitution weight:\t\t", $rounded_sub_weight, "\n"); 
            }elsif ($highest_percentage eq $percent_del){ 
                $del_weight = 1; 
                if ($percent_ins != 0){ 
                    $ins_weight = ($percent_del/$percent_ins);
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight);  
                }
                if ($percent_sub != 0){ 
                    $sub_weight = ($percent_del/$percent_sub);
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight);
                }
                if ($percent_ins == 0){ 
                    $ins_weight = $del_weight + $sub_weight + 1; 
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight);
                }
                if ($percent_sub == 0 ){ 
                    $sub_weight = $del_weight + $ins_weight + 1;
                    $rounded_sub_weight = sprintf("%.2f", $sub_weight); 
                }
                push(@output, "Insertion weight:\t\t", $rounded_ins_weight, "\tDeletion weight:\t", $del_weight, "\tSubstitution weight:\t\t", $rounded_sub_weight, "\n"); 
            }elsif ($highest_percentage eq $percent_sub){ 
                $sub_weight = 1; 
                if ($percent_ins != 0){ 
                    $ins_weight = ($percent_sub/$percent_ins);
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight); 
                }
                if ($percent_del != 0){ 
                    $del_weight = ($percent_sub/$percent_del);
                    $rounded_del_weight = sprintf("%.2f", $del_weight);
                }
                if ($percent_ins == 0){ 
                    $ins_weight = $del_weight + $sub_weight + 1; 
                    $rounded_ins_weight = sprintf("%.2f", $ins_weight);
                }
                if ($percent_del == 0 ){ 
                    $del_weight = $sub_weight + $ins_weight + 1; 
                    $rounded_del_weight = sprintf("%.2f", $del_weight);
                }
                push(@output, "Insertion weight:\t\t", $rounded_ins_weight, "\tDeletion weight:\t", $rounded_del_weight, "\tSubstitution weight:\t\t", $sub_weight, "\n"); 
            }

        }
        push(@output, "Experimental sequence:\n", $new_sequence, "\n"); 
        push(@output, "Reference sequence:\n", $ref, "\n"); 

        _writeOut("Printing Output");
        open(OUT, '>', $OUTPUT) or die $!; 
        print OUT @output; 
        _writeOut("Finished");
    }
}

sub _levdist {
	my ($seq1, $seq2, $expected_number_of_errors, $tolerated_frame_shift) = @_;

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
	
	my @s1_old = split '', $seq1;
    my @s2_old = split '', $seq2;
    # Remove prefixes from sequences
    foreach (my $i = 0; $i < @s1_old; $i++){
        if ($s1_old[$i] eq $s2_old[$i]){
            undef $s1_old[$i];
            undef $s2_old[$i];
        }elsif ($s1_old[$i] ne $s2_old[$i]){
            last;
        }
    }
    my @s1 = grep{ defined } @s1_old;
    my @s2 = grep{ defined } @s2_old;
    my $l1 = @s1;
    my $l2 = @s2;
    my $dist;
	# Initialize matrix borders
	my $number_of_rows_before_diagonal = $expected_number_of_errors; 
	my $distance_threshold = $expected_number_of_errors + $tolerated_frame_shift;
    if ($expected_number_of_errors <= $l1){
        for (my $i = 0; $i <= $number_of_rows_before_diagonal; $i++) {
		    $dist->[$i]->[0] = $i;
	    }
        for (my $j = 0; $j <= $number_of_rows_before_diagonal; $j++) {
            $dist->[0]->[$j] = $j;
        }
    }elsif ($expected_number_of_errors > $l1){
        for (my $i = 0; $i <= $l1; $i++){
            $dist->[$i]->[0] = $i;
        }
    	for (my $j = 0; $j <= $l1; $j++){
        	$dist->[0]->[$j] = $j;
        }
    }
	my $terminal_position = $number_of_rows_before_diagonal - 1; 
	my $j_upstream_border = 0;
    my $j_downstream_border = $expected_number_of_errors*2+2;
    my $j_start = 1;
    my $j_end = $expected_number_of_errors*2+1;

	# Fill in matrix table with Levenshtein distance calculations
	for (my $i = 1; $i <= $l1; $i++){
		$terminal_position++; 
		if ($i <= $number_of_rows_before_diagonal){
			if ($terminal_position < $l2){  
				for (my $j = 1; $j <= $terminal_position; $j++){ 
					my $cost;
					if ( $s1[$i-1] eq $s2[$j-1] ) {
						$cost = 0;
					} else {
						$cost = 1;
					}
					$dist->[$i]->[$j] = min($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1 );
					$dist->[$i]->[$terminal_position + 1] = $number_of_rows_before_diagonal; 
					if ($i == $j){ 
						if ($dist->[$i]->[$j] > $distance_threshold){ 
							return; 
						} 
					} 
				}
			} elsif ($terminal_position >= $l2){ 
				my $j_end_and_string_length_difference = $terminal_position - $l2;
                my $new_j_end = $terminal_position - $j_end_and_string_length_difference;		
				for (my $j = 1; $j <= $new_j_end; $j++){
                	my $cost;
                    if ( $s1[$i-1] eq $s2[$j-1] ) {
                        $cost = 0;
                    } else {
                    	$cost = 1;
                    }
                    $dist->[$i]->[$j] = min($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1);
					if ($i == $j){
                    	if ($dist->[$i]->[$j] > $distance_threshold){
                        	return; 
                        }
                    }
                }
			}	
        }elsif ($i > $number_of_rows_before_diagonal){ 
			$j_upstream_border++;
            $j_downstream_border++;
            $j_start++;
            $j_end++;
			if ($j_end < $l2){ 
				for (my $j = $j_start; $j <= $j_end; $j++){ 
					$dist->[$i]->[$j_upstream_border] = $number_of_rows_before_diagonal;
					my $cost;
                	if ( $s1[$i-1] eq $s2[$j-1] ) {
                    	$cost = 0;
                    } else {
                        $cost = 1;
                    }
                    $dist->[$i]->[$j] = min($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1);
                	$dist->[$i]->[$j_downstream_border] = $number_of_rows_before_diagonal;
					if ($i == $j){
                    	if ($dist->[$i]->[$j] > $distance_threshold){ 
							return; 
                        }
                    }
				}
			} elsif ($j_end >= $l2){ 
				my $j_end_and_string_length_difference = $terminal_position - $l2;
            	my $new_j_end = $terminal_position - $j_end_and_string_length_difference;
				for (my $j = $j_start; $j <= $new_j_end; $j++){
					$dist->[$i]->[$j_upstream_border] = $number_of_rows_before_diagonal;
                    my $cost;
                    if ( $s1[$i-1] eq $s2[$j-1] ) {
                    	$cost = 0;
                    } else {
                    	$cost = 1;
                	}
                    $dist->[$i]->[$j] = min($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1);
					if ($i == $j){
                    	if ($dist->[$i]->[$j] > $distance_threshold){ 
							return; 
                        }
                    }
				}	
			}
		}   
    }    
	# Calculate error types 
	my $classic_distance = $dist->[$l1]->[$l2];
	return $classic_distance if ( $classic_distance == 0 ); # There is no difference between strings

	my $vert_min_dist;
	if ($expected_number_of_errors > $l1){ 
		if ($l1 >= 1){ 
			$vert_min_dist = $dist->[0]->[$l2];
			my $last_i = 0; 
			my $start_i = $l1-1; 
			for my $i ($last_i..$start_i){
				$vert_min_dist = min($vert_min_dist, $dist->[$i]->[$l2]); 
			}
		} elsif ($l1 <= 1) { 
			$vert_min_dist = $dist->[0]->[$l2];	
			my $last_i = 0;
            my $start_i = $l1;
            for my $i ($last_i..$start_i){
            	$vert_min_dist = min($vert_min_dist, $dist->[$i]->[$l2]);
            }
		}
	} else { 
		if ($l1 >= 1){ 
			$vert_min_dist = $dist->[$l1-$expected_number_of_errors]->[$l2];
			my $last_i = $l1-$expected_number_of_errors; 
			my $start_i = $l1-1; 
			for my $i ($last_i..$start_i){ 
				$vert_min_dist = min($vert_min_dist, $dist->[$i]->[$l2]); 
			}
		} elsif ($l1 <= 1) { 
			$vert_min_dist = $dist->[$l1-$expected_number_of_errors]->[$l2];
			my $last_i = $l1-$expected_number_of_errors;
            my $start_i = $l1;
            for my $i ($last_i..$start_i){
            	$vert_min_dist = min($vert_min_dist, $dist->[$i]->[$l2]);
            }
		}
	} 
	my $hor_min_dist; 
	if ($expected_number_of_errors > $l2){
    	if ($l2 >= 1){
			$hor_min_dist = $dist->[$l1]->[0];
            my $last_j = 0;
            my $start_j = $l2-1;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = min($hor_min_dist, $dist->[$l1]->[$j]);
            }
        } elsif ($l2 <= 1) {
			$hor_min_dist = $dist->[$l1]->[0];
            my $last_j = 0;
            my $start_j = $l2;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = min($hor_min_dist, $dist->[$l1]->[$j]);
            }
        }
    } elsif ($expected_number_of_errors <= $l2){
    	if ($l2 >= 1){
			$hor_min_dist = $dist->[$l1]->[$l2-$expected_number_of_errors];	
            my $last_j = $l1-$expected_number_of_errors;
            my $start_j = $l2-1;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = min($hor_min_dist, $dist->[$l1]->[$j]);
            }
    	} elsif ($l2 <= 1) {
			$hor_min_dist = $dist->[$l1]->[$l2-$expected_number_of_errors];
            my $last_j = $l1-$expected_number_of_errors;
            my $start_j = $l2;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = min($hor_min_dist, $dist->[$l1]->[$j]);
            }
    	}
    }

	my $corner_min_dist = $classic_distance;
	my $modified_distance = min($vert_min_dist, $hor_min_dist, $corner_min_dist);

	my $insertion_sum; 
    my $deletion_sum; 
    my $substitution_sum;  

	# If the minimum distance is found on the vertical border...
	if (($vert_min_dist < $corner_min_dist && $vert_min_dist < $hor_min_dist) ||($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $INS eq 'more') || ($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $INS eq 'same')) { 		
		# Find the location of the modified minimum distance from the corner [i][j]
		for ( my $i = $l1-1; $i >= ($l1-($expected_number_of_errors+1)); $i-- ) {
			if ($dist->[$i]->[$l2] == $modified_distance) {
				# If the modified distance is not consecutive
                if ($dist->[$i-1]->[$l2] != $dist->[$i]->[$l2]){
                    my $size_of_s2 = @s2; 
					my $distance_from_corner = ($size_of_s2 - $i); 
					# If only insertions are present
					if ($distance_from_corner == $modified_distance){ 
						$insertion_sum += (1*$distance_from_corner); 
					# If insertions are mixed with other error types
					}elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_insertions = $distance_from_corner; 
						$insertion_sum += (1*$base_number_of_insertions);
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
                        # If there is an even number of errors, there will either be subs or
                        # indel pairs based on probability 
						if ($remaining_steps % 2 == 0){
                            if ($SUB eq 'less'){ 
                                $insertion_sum += (1 * ($remaining_steps/2));
                                $deletion_sum += (1 * ($remaining_steps/2)); 
                            }elsif ($SUB eq 'more'){ 
                                $substitution_sum += (1 * $remaining_steps); 
                            }elsif ($SUB eq 'same'){ 
                                if ($remaining_steps > 2){ 
                                    my $half_remaining_steps = ($remaining_steps/2); 
                                    $substitution_sum += $half_remaining_steps; 
                                    
                                    if ($half_remaining_steps % 2 == 0){
                                        $insertion_sum += (1 * ($half_remaining_steps/2)); 
                                        $deletion_sum += (1 * ($half_remaining_steps/2)); 
                                    }else{ 
                                        if ($INS eq 'more' || $INS eq 'same'){ 
                                            $insertion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += floor(1 * ($half_remaining_steps/2)); 
                                        }else{ 
                                            $insertion_sum += floor(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                        }
                                    }
                                }elsif ($remaining_steps == 2){ 
                                    if ($INS eq 'more' || $INS eq 'same'){ 
                                        $insertion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }else{ 
                                        $deletion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }
                                }
                            }
                        # If there is an odd number of remaining steps, greater than one (i.e. 3, 5, 7, etc.)
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){ 
							# Every remainder is a substitution. every dividend is an insertion-deletion pair. 
							use integer; 
							my $number_of_indel_pairs = $remaining_steps/2; # will have a decimal
                            if ($SUB eq 'less'){
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += ceil(1 * ($number_of_indel_pairs/2));
                                }
                            }elsif ($SUB eq 'more'){
                                $substitution_sum += $remaining_steps; 
                            }elsif ($SUB eq 'same'){ 
                                # need to split an odd number into sub, ind, del "evenly" 
                                $substitution_sum += (1 * $number_of_indel_pairs); 
                                my $to_be_split_between_ins_and_del = ($remaining_steps - $number_of_indel_pairs); 
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2));
                                }
                            }
                        # If there is only one remaining step
						} elsif ($remaining_steps == 1){ 
                            # CAN ONLY BE A SUBSTITUTION
							$substitution_sum += 1;  
						}
					}	
					last; 
				# If the modified distance is consecutive
                } elsif ($dist->[$i-1]->[$l2] == $dist->[$i]->[$l2]) {
                    my $size_of_s2 = @s2; 
                	my $distance_from_corner = ($size_of_s2 - $i);
                    # If only insertions are present
					if ($distance_from_corner == $modified_distance){
                        $insertion_sum += (1 * $distance_from_corner);
                    	# If insertions are mixed with other error types
					}elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_insertions = $distance_from_corner; 
						$insertion_sum += (1*$base_number_of_insertions);
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
                        # If there is an even number of errors, there will either be subs or
                        # indel pairs based on probability 
						if ($remaining_steps % 2 == 0){
                            if ($SUB eq 'less'){ 
                                $insertion_sum += (1 * ($remaining_steps/2));
                                $deletion_sum += (1 * ($remaining_steps/2)); 
                            }elsif ($SUB eq 'more'){ 
                                $substitution_sum += (1 * $remaining_steps); 
                            }elsif ($SUB eq 'same'){ 
                                if ($remaining_steps > 2){ 
                                    my $half_remaining_steps = ($remaining_steps/2); 
                                    $substitution_sum += $half_remaining_steps; 
                                    
                                    if ($half_remaining_steps % 2 == 0){ 
                                        $insertion_sum += (1 * ($half_remaining_steps/2)); 
                                        $deletion_sum += (1 * ($half_remaining_steps/2)); 
                                    }else{ 
                                        if ($INS eq 'more' || $INS eq 'same'){ 
                                            $insertion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += floor(1 * ($half_remaining_steps/2)); 
                                        }else{ 
                                            $insertion_sum += floor(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                        }
                                    }
                                }elsif ($remaining_steps == 2){ 
                                    if ($INS eq 'more' || $INS eq 'same'){ 
                                        $insertion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }else{ 
                                        $deletion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }
                                }
                            }
                        # If there is an odd number of remaining steps, greater than one (i.e. 3, 5, 7, etc.)
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){ 
							# Every remainder is a substitution. every dividend is an insertion-deletion pair. 
							use integer; 
							my $number_of_indel_pairs = $remaining_steps/2; # will have a decimal
                            if ($SUB eq 'less'){
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += ceil(1 * ($number_of_indel_pairs/2));
                                }
                            }elsif ($SUB eq 'more'){
                                $substitution_sum += $remaining_steps; 
                            }elsif ($SUB eq 'same'){ 
                                # need to split an odd number into sub, ind, del "evenly" 
                                $substitution_sum += (1 * $number_of_indel_pairs); 
                                my $to_be_split_between_ins_and_del = ($remaining_steps - $number_of_indel_pairs); 
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2));
                                }
                            }
                        # If there is only one remaining step
						} elsif ($remaining_steps == 1){ 
                            # CAN ONLY BE A SUBSTITUTION
							$substitution_sum += 1;  
						}
					}
                    last;
                }
	    	}
		} 
	# If the minimum distance is found on the horizontal border...
	} elsif (($hor_min_dist < $corner_min_dist && $hor_min_dist < $vert_min_dist) ||($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $INS eq 'less')){
		# Find the location of the modified minimum distance from corner [i][j]
		for ( my $j = $l2-1; $j >= ($l2-($expected_number_of_errors+1)); $j-- ){ 
			if ( $dist->[$l1]->[$j] == $modified_distance) {
				# If the modified distance is not consecutive
                if ($dist->[$l1]->[$j] != $dist->[$l1]->[$j-1]){
                    my $size_of_s2 = @s2; 
					my $distance_from_corner = ($size_of_s2 - $j); 
					# If only deletions are present 
					if ($distance_from_corner == $modified_distance){ 
                       $deletion_sum += (1 * $distance_from_corner);
					# If deletions are mixed with other error types 
					} elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_deletions = $distance_from_corner; 
                        $deletion_sum += (1 * $base_number_of_deletions);
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
                        # If there is an even number of errors, there will either be subs or
                        # indel pairs based on probability 
						if ($remaining_steps % 2 == 0){
                            if ($SUB eq 'less'){ 
                                $insertion_sum += (1 * ($remaining_steps/2));
                                $deletion_sum += (1 * ($remaining_steps/2)); 
                            }elsif ($SUB eq 'more'){ 
                                $substitution_sum += (1 * $remaining_steps); 
                            }elsif ($SUB eq 'same'){ 
                                if ($remaining_steps > 2){ 
                                    my $half_remaining_steps = ($remaining_steps/2); 
                                    $substitution_sum += $half_remaining_steps; 
                                    
                                    if ($half_remaining_steps % 2 == 0){ 
                                        $insertion_sum += (1 * ($half_remaining_steps/2)); 
                                        $deletion_sum += (1 * ($half_remaining_steps/2)); 
                                    }else{ 
                                        if ($INS eq 'more' || $INS eq 'same'){ 
                                            $insertion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += floor(1 * ($half_remaining_steps/2)); 
                                        }else{ 
                                            $insertion_sum += floor(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                        }
                                    }
                                }elsif ($remaining_steps == 2){ 
                                    if ($INS eq 'more' || $INS eq 'same'){ 
                                        $insertion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }else{ 
                                        $deletion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }
                                }
                            }
                        # If there is an odd number of remaining steps, greater than one (i.e. 3, 5, 7, etc.)
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){ 
							# Every remainder is a substitution. every dividend is an insertion-deletion pair. 
							use integer; 
							my $number_of_indel_pairs = $remaining_steps/2; # will have a decimal
                            if ($SUB eq 'less'){
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += ceil(1 * ($number_of_indel_pairs/2));
                                }
                            }elsif ($SUB eq 'more'){
                                $substitution_sum += $remaining_steps; 
                            }elsif ($SUB eq 'same'){ 
                                # need to split an odd number into sub, ind, del "evenly" 
                                $substitution_sum += (1 * $number_of_indel_pairs); 
                                my $to_be_split_between_ins_and_del = ($remaining_steps - $number_of_indel_pairs); 
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2));
                                }
                            }
                        # If there is only one remaining step
						} elsif ($remaining_steps == 1){ 
                            # CAN ONLY BE A SUBSTITUTION
							$substitution_sum += 1;  
						}
					}
                    last; 
				# If the modified distance is consecutive
                } elsif ( $dist->[$l1]->[$j] == $dist->[$l1]->[$j-1]) {
                    my $size_of_s2 = @s2; 
                    my $distance_from_corner = ($size_of_s2 - $j);
                    # If only deletions are present
                    if ($distance_from_corner == $modified_distance){
                        $deletion_sum += (1 * $distance_from_corner); 
                    # If deletions are mixed with other error types
                    }elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_deletions = $distance_from_corner; 
						$deletion_sum += (1*$base_number_of_deletions);
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
                        # If there is an even number of errors, there will either be subs or
                        # indel pairs based on probability 
						if ($remaining_steps % 2 == 0){
                            if ($SUB eq 'less'){ 
                                $insertion_sum += (1 * ($remaining_steps/2));
                                $deletion_sum += (1 * ($remaining_steps/2)); 
                            }elsif ($SUB eq 'more'){ 
                                $substitution_sum += (1 * $remaining_steps); 
                            }elsif ($SUB eq 'same'){ 
                                if ($remaining_steps > 2){ 
                                    my $half_remaining_steps = ($remaining_steps/2); 
                                    $substitution_sum += $half_remaining_steps; 
                                    
                                    if ($half_remaining_steps % 2 == 0){ 
                                        $insertion_sum += (1 * ($half_remaining_steps/2)); 
                                        $deletion_sum += (1 * ($half_remaining_steps/2)); 
                                    }else{ 
                                        if ($INS eq 'more' || $INS eq 'same'){ 
                                            $insertion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += floor(1 * ($half_remaining_steps/2)); 
                                        }else{ 
                                            $insertion_sum += floor(1 * ($half_remaining_steps/2)); 
                                            $deletion_sum += ceil(1 * ($half_remaining_steps/2)); 
                                        }
                                    }
                                }elsif ($remaining_steps == 2){ 
                                    if ($INS eq 'more' || $INS eq 'same'){ 
                                        $insertion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }else{ 
                                        $deletion_sum += 1; 
                                        $substitution_sum += 1; 
                                    }
                                }
                            }
                        # If there is an odd number of remaining steps, greater than one (i.e. 3, 5, 7, etc.)
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){ 
							# Every remainder is a substitution. every dividend is an insertion-deletion pair. 
							use integer; 
							my $number_of_indel_pairs = $remaining_steps/2; # will have a decimal
                            if ($SUB eq 'less'){
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($number_of_indel_pairs/2)); 
                                    $deletion_sum += ceil(1 * ($number_of_indel_pairs/2));
                                }
                            }elsif ($SUB eq 'more'){
                                $substitution_sum += $remaining_steps; 
                            }elsif ($SUB eq 'same'){ 
                                # need to split an odd number into sub, ind, del "evenly" 
                                $substitution_sum += (1 * $number_of_indel_pairs); 
                                my $to_be_split_between_ins_and_del = ($remaining_steps - $number_of_indel_pairs); 
                                if ($INS eq 'more' || $INS eq 'same'){ 
                                    $insertion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                }else{ 
                                    $insertion_sum += floor(1 * ($to_be_split_between_ins_and_del/2)); 
                                    $deletion_sum += ceil(1 * ($to_be_split_between_ins_and_del/2));
                                }
                            }
                        # If there is only one remaining step
						} elsif ($remaining_steps == 1){ 
                            # CAN ONLY BE A SUBSTITUTION
							$substitution_sum += 1;  
						}
					}
                    last;
                } 
			} 
		} 

	# If the minimum distance is found in the corner [i][j]... 
	} elsif (($corner_min_dist <= $vert_min_dist) || ($corner_min_dist <= $hor_min_dist)){ 
		if ($corner_min_dist > 1){ 
			# If the number at [i][j] is even, calculate the lowest weights. Will either be all subs or indels or a mix of the two.
			if ($corner_min_dist % 2 == 0){ 
                if ($SUB eq 'less'){ 
                    $insertion_sum += (1 * ($corner_min_dist/2)); 
                    $deletion_sum += (1 * ($corner_min_dist/2));
                }elsif ($SUB eq 'more'){ 
                    $substitution_sum += (1 * $corner_min_dist); 
                }elsif ($SUB eq 'same'){ 
                    # If the corner value is 2,4,6,8,10,etc...
                    if ($corner_min_dist > 2){
                        my $half_corner_value = ($corner_min_dist/2); 
                        $substitution_sum += $half_corner_value; 

                        if ($half_corner_value % 2 == 0){ 
                            $insertion_sum += (1 * ($half_corner_value/2)); 
                            $deletion_sum += (1 * ($half_corner_value/2)); 
                        }else{
                            if ($INS eq 'more' || $INS eq 'same'){ 
                                $insertion_sum += ceil(1 * ($half_corner_value/2)); 
                                $deletion_sum += floor(1 * ($half_corner_value/2)); 
                            }else{ 
                                $insertion_sum += floor(1 * ($half_corner_value/2)); 
                                $deletion_sum += ceil(1 * ($half_corner_value/2)); 
                            }
                        }
                    } elsif ($corner_min_dist == 2 ){ 
                        if ($INS eq 'more' || $INS eq 'same'){ 
                            $insertion_sum += 1; 
                            $substitution_sum += 1; 
                        } else { 
                            $deletion_sum += 1; 
                            $substitution_sum += 1; 
                        }
                    }
                }
        	}
        	# If the number at [i][j] is odd, calculate the number of indels and subs needed.
        	if ($corner_min_dist % 2 == 1){
                if ( $SUB eq 'less' ){ 
                    $insertion_sum += (1 * (($corner_min_dist/2) - 0.5)); 
                    $deletion_sum += (1 * (($corner_min_dist/2) - 0.5)); 
                    $substitution_sum += 1;
                }elsif ($SUB eq 'same') { 
                    $insertion_sum += (1 * (($corner_min_dist/2) - 0.5)); 
                    $deletion_sum += (1 * (($corner_min_dist/2) - 0.5)); 
                    $substitution_sum += 1;
                } elsif ($SUB eq 'more'){
                    $substitution_sum += (1 * $corner_min_dist);
                }
        	}
		} elsif ($corner_min_dist == 1){
            # If there is only one edit distance at [i][j] it is a substitution
        	$substitution_sum += 1;  
		}
	}
	
	# Return
    
    #$insertion_sum // $insertion_sum = 0; 
    #$deletion_sum // $deletion_sum = 0; 
    #$substitution_sum // $substitution_sum = 0; 
    
    if (!defined $insertion_sum){ 
        $insertion_sum = 0;
    } 
    if (!defined $deletion_sum){ 
        $deletion_sum = 0;
    } 
    if (!defined $substitution_sum){ 
        $substitution_sum = 0;
    } 
    my $sum = $insertion_sum + $deletion_sum + $substitution_sum; 
	# If we have nothing, just return.
	return if ($sum == 0);

    my $error_sum = ($insertion_sum + $deletion_sum + $substitution_sum); 
    my @results = ($insertion_sum, "\t", $deletion_sum, "\t", $substitution_sum, "\t", $error_sum);
	return @results;
}