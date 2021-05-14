#!/usr/bin/perl 
use warnings; 
use strict; 

# Collects weights for weight threshold determination. 
# To be used on biological datasets
# usage inputs: ./Weight_Threshold_Test_Biological_Dataset.pl < Multifasta of all simulated sequences for searching > < output file name >

# Will need to manually put in the error profiles as well as the name of the output file. 

my $substitution_weight = 3; 
my $insertion_weight = 1; 
my $deletion_weight = 1; 
my $reference_sequence = "CACAAAGACACCGACAACTTTCTTCAGCACCTAAGAAAGTTGTCGGTGTCTTTGTG"; 

my $file = $ARGV[0]; 
open(IN, '<', $file) or die $!; 
my @output;
while (my $line = <IN>){ 
    chomp $line; 
    if ($line =~ /(^>)(.*)/){   
        push(@output, $line, "\n");
    } else { 
        my $windowstart; 
        my $window_size = 56; 
        my $step_size = 1; 
        my $length = length($line);  
        for ($windowstart = 0; $windowstart <= $length; $windowstart += $step_size){ 
            my $seq_capture = substr($line, $windowstart, $window_size);
            my $error_threshold = 3; 
            my $frame_shift = 11; # allows a 20% frameshift 
            my $weight_threshold = 500; # Chose an improbabily high weight theshold, which effectively makes for "no threshold"
            my $rev_search_probe = reverse $reference_sequence;
            my $rev_seq_capture = reverse $seq_capture; 
            my ($dist_1, $dist_2) = levdist($seq_capture, $reference_sequence, $error_threshold, $frame_shift);
			my ($rev_dist_1, $rev_dist_2) = levdist($rev_search_probe, $rev_seq_capture, $error_threshold, $frame_shift);
            if ((defined $dist_1) && (defined $dist_2) && (defined $rev_dist_1) && (defined $rev_dist_2)){ 
                if (($dist_1 < ($weight_threshold+1)) || ($dist_2 < ($weight_threshold+1)) || ($rev_dist_1 < ($weight_threshold+1)) || ($rev_dist_2 < ($weight_threshold+1))){
                    my @distances;
					push(@distances, $dist_1, $dist_2, $rev_dist_1, $rev_dist_2);
					my @sorted = sort { $a <=> $b } @distances;
                    if ($sorted[0] <= $weight_threshold){ 
                        push(@output, $seq_capture, "\t", $sorted[0], "\n"); 
                    }
                }
            }
        }
    }
}
close IN;
open(OUT, '>', "$ARGV[1]") or die $!; 
print OUT @output; 

sub levdist {
	my @err; 
	my @bidir_err; 
	my ( $seq1, $seq2, $expected_number_of_errors, $tolerated_frame_shift) = (@_)[0,1,2,3];
	if ($seq1 =~ $seq2){ 
        return (0, 0); 
    }
	my $l1_old = length($seq1);
    my $l2_old = length($seq2);
    if ($l1_old != $l2_old){
        #print "Error: strings ", $seq1, " and ", $seq2, " are not the same size\n";
        return;
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
	my $number_of_rows_before_diagonal = $expected_number_of_errors + 1; 
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
					$dist->[$i]->[$j] = minimum($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1 );
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
                    $dist->[$i]->[$j] = minimum($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1);
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
                    $dist->[$i]->[$j] = minimum($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1);
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
                    $dist->[$i]->[$j] = minimum($dist->[$i-1]->[$j-1] + $cost, $dist->[$i]->[$j-1]+1, $dist->[$i-1]->[$j]+1);
					if ($i == $j){
                    	if ($dist->[$i]->[$j] > $distance_threshold){ 
							return; 
                        }
                    }
				}	
			}
		} 
	} 
	# Calculate weighted errors 
	my $classic_distance = $dist->[$l1]->[$l2];
	if ( $classic_distance == 0 ) {
		return $classic_distance;
		# There is no difference between strings
	}
	my $vert_min_dist;
	if ($expected_number_of_errors > $l1){ 
		if ($l1 >= 1){ 
			$vert_min_dist = $dist->[0]->[$l2];
			my $last_i = 0; 
			my $start_i = $l1-1; 
			for my $i ($last_i..$start_i){
				$vert_min_dist = minimum($vert_min_dist, $dist->[$i]->[$l2]); 
			}
		} elsif ($l1 <= 1){ 
			$vert_min_dist = $dist->[0]->[$l2];	
			my $last_i = 0;
            my $start_i = $l1;
            for my $i ($last_i..$start_i){
            	$vert_min_dist = minimum($vert_min_dist, $dist->[$i]->[$l2]);
            }
		}
	} elsif ($expected_number_of_errors <= $l1){ 
		if ($l1 >= 1){ 
			$vert_min_dist = $dist->[$l1-$expected_number_of_errors]->[$l2];
			my $last_i = $l1-$expected_number_of_errors; 
			my $start_i = $l1-1; 
			for my $i ($last_i..$start_i){ 
				$vert_min_dist = minimum($vert_min_dist, $dist->[$i]->[$l2]); 
			}
		}elsif ($l1 <= 1){ 
			$vert_min_dist = $dist->[$l1-$expected_number_of_errors]->[$l2];
			my $last_i = $l1-$expected_number_of_errors;
            my $start_i = $l1;
            for my $i ($last_i..$start_i){
            	$vert_min_dist = minimum($vert_min_dist, $dist->[$i]->[$l2]);
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
            	$hor_min_dist = minimum($hor_min_dist, $dist->[$l1]->[$j]);
            }
        } elsif ($l2 <= 1){
			$hor_min_dist = $dist->[$l1]->[0];
            my $last_j = 0;
            my $start_j = $l2;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = minimum($hor_min_dist, $dist->[$l1]->[$j]);
            }
        }
    } elsif ($expected_number_of_errors <= $l2){
    	if ($l2 >= 1){
			$hor_min_dist = $dist->[$l1]->[$l2-$expected_number_of_errors];	
            my $last_j = $l1-$expected_number_of_errors;
            my $start_j = $l2-1;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = minimum($hor_min_dist, $dist->[$l1]->[$j]);
            }
    	}elsif ($l2 <= 1){
			$hor_min_dist = $dist->[$l1]->[$l2-$expected_number_of_errors];
            my $last_j = $l1-$expected_number_of_errors;
            my $start_j = $l2;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = minimum($hor_min_dist, $dist->[$l1]->[$j]);
            }
    	}
    }
	my $corner_min_dist = $classic_distance;
	my $modified_distance = minimum($vert_min_dist, $hor_min_dist, $corner_min_dist);
	# If the minimum distance is found on the vertical border...
	if (( $vert_min_dist < $corner_min_dist && $vert_min_dist < $hor_min_dist)  || ($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $deletion_weight >= $insertion_weight)) { 
		1;
		# Find the location of the modified minimum distance from the corner [i][j]
		for ( my $i = $l1-1; $i >= ($l1-($expected_number_of_errors+1)); $i-- ) {
			if ($dist->[$i]->[$l2] == $modified_distance) {
				# If the modified distance is not consecutive
                if ($dist->[$i-1]->[$l2] != $dist->[$i]->[$l2]){
					my $size_of_s2 = @s2; 
					my $distance_from_corner = ($size_of_s2 - ($i)); 
					# If only insertions are present
					if ($distance_from_corner == $modified_distance){ 
						push (@err, $insertion_weight*$distance_from_corner); 
						push (@bidir_err, $deletion_weight*$distance_from_corner); 
					# If insertions are mixed with other error types
					}elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_insertions = $distance_from_corner; 
						push(@err, $insertion_weight*$base_number_of_insertions); 
						push(@bidir_err, $deletion_weight*$base_number_of_insertions); 
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
						if ($remaining_steps % 2 == 0){
							if ($substitution_weight > $insertion_weight + $deletion_weight){ 
								my $number_of_indel_pairs = $remaining_steps/2; 
								push(@err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
								push(@bidir_err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
							}elsif ($substitution_weight < $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}elsif ($substitution_weight = $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){ 
							# every remainder is a substitution. every dividend is an insertion-deletion pair. 
							use integer; 
							my $number_of_indel_pairs = $remaining_steps/2; 
							if ($number_of_indel_pairs == 0){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							} else{ 
								push(@err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight)); 
								push(@bidir_err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight));
							}
						} elsif ($remaining_steps == 1){ 
							push (@err, $substitution_weight); 
							push (@bidir_err, $substitution_weight); 
						}
					}	
					last; 
				# If the modified distance is consecutive
                } elsif ($dist->[$i-1]->[$l2] == $dist->[$i]->[$l2]){
                	my $size_of_s2 = @s2;
                	my $distance_from_corner = ($size_of_s2 - ($i));
                    # If only insertions are present
					if ($distance_from_corner == $modified_distance){
                    	push (@err, $insertion_weight*$distance_from_corner);
						push (@bidir_err, $deletion_weight*$distance_from_corner); 
                    # If insertions are mixed with other error types
                    }elsif ($distance_from_corner < $modified_distance){
                    	my $base_number_of_insertions = $distance_from_corner;
                    	push(@err, $insertion_weight*$base_number_of_insertions);
						push(@bidir_err, $deletion_weight*$base_number_of_insertions); 
                        my $remaining_steps = ($modified_distance - $distance_from_corner);
						if ($remaining_steps % 2 == 0){
                        	if ($substitution_weight > $insertion_weight + $deletion_weight){ 
								my $number_of_indel_pairs = $remaining_steps/2; 
								push(@err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
								push(@bidir_err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
							}elsif ($substitution_weight < $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}elsif ($substitution_weight = $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}
                        } elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){
                        	# every remainder is a substitution. every dividend is an insertion-deletion pair.
                            use integer;
                            my $number_of_indel_pairs = $remaining_steps/2;
					    	push(@err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight));
							push(@bidir_err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight));
                        } elsif ($remaining_steps == 1){
                            push (@err, $substitution_weight);
							push(@bidir_err, $substitution_weight); 
                        }
                    }
                    last;
                }

	    	}
		} 
	# If the minimum distance is found on the horizontal border...
	} elsif (($hor_min_dist < $corner_min_dist && $hor_min_dist < $vert_min_dist) || ($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $insertion_weight >= $deletion_weight)){
		1;
		# Find the location of the modified minimum distance from corner [i][j]
		for ( my $j = $l2-1; $j >= ($l2-($expected_number_of_errors+1)); $j-- ){ 
			if ( $dist->[$l1]->[$j] == $modified_distance) {
				# If the modified distance is not consecutive
                if ($dist->[$l1]->[$j] != $dist->[$l1]->[$j-1]){
					my $size_of_s2 = @s2; 
					my $distance_from_corner = ($size_of_s2 - ($j)); 
					# If only deletions are present 
					if ($distance_from_corner == $modified_distance){ 
						push(@err, ($deletion_weight*$distance_from_corner)); 
						push(@bidir_err, ($insertion_weight*$distance_from_corner)); 
					# If deletions are mixed with other error types 
					} elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_deletions = $distance_from_corner; 
						push(@err, ($deletion_weight*$base_number_of_deletions)); 
						push(@bidir_err, ($insertion_weight*$base_number_of_deletions)); 
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
						if ($remaining_steps %2 == 0){ 
							if ($substitution_weight > $insertion_weight + $deletion_weight){ 
								my $number_of_indel_pairs = $remaining_steps/2; 
								push(@err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
								push(@bidir_err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
							}elsif ($substitution_weight < $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}elsif ($substitution_weight = $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}
						} elsif ($remaining_steps == 1){ 
							push (@err, $substitution_weight);
							push (@bidir_err, $substitution_weight);
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){
                            # every remainder is a substitution. every dividend is an insertion-deletion pair.
                            use integer;
                            my $number_of_indel_pairs = $remaining_steps/2;
                            push(@err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight)); 
							push (@bidir_err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight)); 
						} 
					} 
					last; 
				# If the modified distance is consecutive
                } elsif ( $dist->[$l1]->[$j] == $dist->[$l1]->[$j-1]){
                	my $size_of_s2 = @s2;
                    my $distance_from_corner = ($size_of_s2 - ($j));
                    # If only deletions are present
                    if ($distance_from_corner == $modified_distance){
                    	push(@err, ($deletion_weight*$distance_from_corner));
						push(@bidir_err, ($insertion_weight*$distance_from_corner)); 
                    # If deletions are mixed with other error types
                    } elsif ($distance_from_corner < $modified_distance){
                        my $base_number_of_deletions = $distance_from_corner;
                        push(@err, ($deletion_weight*$base_number_of_deletions));
						push(@bidir_err, ($insertion_weight*$base_number_of_deletions)); 
                        my $remaining_steps = ($modified_distance - $distance_from_corner);
                        if ($remaining_steps %2 == 0){
                            if ($substitution_weight > $insertion_weight + $deletion_weight){ 
								my $number_of_indel_pairs = $remaining_steps/2; 
								push(@err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
								push(@bidir_err, (($insertion_weight + $deletion_weight) * $number_of_indel_pairs)); 
							}elsif ($substitution_weight < $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}elsif ($substitution_weight = $insertion_weight + $deletion_weight){ 
								push(@err, $substitution_weight); 
								push(@bidir_err, $substitution_weight); 
							}
                        } elsif ($remaining_steps == 1){
                            push (@err, $substitution_weight);
							push(@bidir_err, $substitution_weight); 
                        } elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){
                        # every remainder is a substitution. every dividend is an insertion-deletion pair.
                            use integer;
                            my $number_of_indel_pairs = $remaining_steps/2;
                            push(@err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight));
							push(@bidir_err, ((($insertion_weight + $deletion_weight) * $number_of_indel_pairs) + $substitution_weight));
                        }
                    }
                    last;
                }
			} 
		} 
	# If the minimum distance is found in the corner [i][j]... 
	} elsif (($corner_min_dist <= $vert_min_dist) || ($corner_min_dist <= $hor_min_dist)){ 
		if ($corner_min_dist > 1){ 
			# If the number at [i][j] is even, calculate the lowest weights. Will either be all subs or indels. 
			if ($corner_min_dist % 2 == 0){ 
				if ($substitution_weight > ($deletion_weight + $insertion_weight)){
					push(@err, (($insertion_weight + $deletion_weight) * ($corner_min_dist/2))); 
					push(@bidir_err, (($insertion_weight + $deletion_weight) * ($corner_min_dist/2))); 
				} elsif ($substitution_weight < ($deletion_weight + $insertion_weight)){
                    push(@err, ($substitution_weight*$corner_min_dist)); 
					push(@bidir_err, ($substitution_weight*$corner_min_dist)); 
                } elsif ($substitution_weight == ($deletion_weight + $insertion_weight)){
                    if ($corner_min_dist > 2){ 
                        my $half_corner_value = ($corner_min_dist/2);
                        push(@err, ($substitution_weight*$half_corner_value)); 
					    push(@bidir_err, ($substitution_weight*$half_corner_value)); 
                        if ($half_corner_value % 2 == 0){ 
                            push(@err, ($insertion_weight*($half_corner_value/2)) + ($deletion_weight*($half_corner_value/2))); 
                            push(@bidir_err, ($insertion_weight*($half_corner_value/2)) + ($deletion_weight*($half_corner_value/2))); 
                        } else { 
                            if ($insertion_weight <= $deletion_weight){ 
                                push(@err, ($insertion_weight * ceil($half_corner_value/2)) + ($deletion_weight * floor($half_corner_value/2))); 
                                push(@bidir_err, ($insertion_weight * ceil($half_corner_value/2)) + ($deletion_weight * floor($half_corner_value/2)));
                            }else{ 
                                push(@err, ($insertion_weight * floor($half_corner_value/2)) + ($deletion_weight * ceil($half_corner_value/2))); 
                                push(@bidir_err, ($insertion_weight * floor($half_corner_value/2)) + ($deletion_weight * ceil($half_corner_value/2)));
                            }
                        }
                    }elsif ($corner_min_dist == 2){
                        if ($insertion_weight <= $deletion_weight){ 
                            push(@err, ($insertion_weight + $substitution_weight)); 
                            push(@bidir_err, ($deletion_weight + $substitution_weight)); 
                        }else{ 
                            push(@err, ($deletion_weight + $substitution_weight)); 
                            push(@bidir_err, ($insertion_weight + $substitution_weight)); 
                        }
                    }
                }
        	}
        	# If the number at [i][j] is odd, calculate the number of indels and subs needed. Preference for lowest weights. 
        	if ($corner_min_dist % 2 == 1){
            	if ($substitution_weight > ($deletion_weight + $insertion_weight)){
                	push(@err, (($insertion_weight + $deletion_weight) * (($corner_min_dist/2) - 0.5)) + $substitution_weight);
					push(@bidir_err, (($insertion_weight + $deletion_weight) * (($corner_min_dist/2) - 0.5)) + $substitution_weight);
                } elsif ($substitution_weight < ($deletion_weight + $insertion_weight)){
					push(@err, ($substitution_weight * $corner_min_dist));
					push(@bidir_err, ($substitution_weight * $corner_min_dist));
                } elsif ($substitution_weight == ($deletion_weight + $insertion_weight)){
                    push(@err, (($insertion_weight + $deletion_weight) * (($corner_min_dist/2) - 0.5)) + $substitution_weight);
					push(@bidir_err, (($insertion_weight + $deletion_weight) * (($corner_min_dist/2) - 0.5)) + $substitution_weight);
                }
        	}
		}elsif ($corner_min_dist == 1){
        # If there is only one edit distance at [i][j] it is a substitution
        	push (@err, $substitution_weight);
			push (@bidir_err, $substitution_weight); 
		}
	}
	my $sum;
	my $bidir_sum; 
	if (!@err){
    	return; 
	}else{
        $sum += $_ for @err;
		$bidir_sum += $_ for @bidir_err; 
	}
	return ($sum, $bidir_sum);
}

sub minimum {
    my $min = shift @_;
    foreach ( @_ ) {
		if ( $_ < $min ) {
	    	$min = $_;
		}
    }
    return $min;
}