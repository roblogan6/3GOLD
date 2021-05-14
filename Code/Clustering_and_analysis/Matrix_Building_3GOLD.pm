package Matrix_Building_3GOLD;

use warnings;
use strict;
use local::lib;
use threads;
use Thread::Pool;
use threads::shared;
use List::Util qw( min );
use List::Util qw( max );
use List::Util qw( reduce );
use List::MoreUtils qw(uniq); 
use Array::Utils qw(:all);
use Array::Compare;
use Data::Dumper;
use Time::HiRes qw(gettimeofday);
use Time::Date; 
use Unix::Processors;
use POSIX; 

# Instantiation and Configuration Methods	

sub new {
    my $class = shift;
    my $self = {
    	'env' => {
			'cpu' => {
				'cores' => undef,
				'threads' => undef,
			},
    	},
    	'config' => {
    		'spawn_threads' => 1,
			'error_threshold' => shift, # Total number of errors allowed between related sequences
			'weight_threshold' => shift, # Cumulative weight of errors allowed between related sequences 
			'tolerated_frame_shift' => shift, # The number of frameshift positions tolerated
			'weight_insertion' => shift, # The weight assigned to insertion errors
			'weight_deletion' => shift, # The weight assigned to deletion errors
			'weight_substitution' => shift, # The weight assigned to substitution errors
			'singletons_out' => shift, # The file name to print out discovered singletons to
			'ave_cluster_distances_out' => shift, # The file name to print out the intracluster distances to
            'clusters_out' => shift, # The file name to print out the clusters of sequences per centroid to
			'size_out' => shift, # The file name to print out the cluster size lookup hash to
			'quiet' => shift,
			'threads' => shift,
    	},
    	'exec_state' => {
    		'return_call_file' => undef,
    		'return_call_package' => undef,
    		'return_call_sub' => undef,
    		'return_call_line' => undef,
    		'return_state_code' => undef,
    		'return_state_name' => undef,
    		'return_state_reason' => undef,
    		'current_exec' => undef,
    	},
    	'data' => {
    		'tag_array' => [],
			'sub_tag_array' => [],
            'array_of_arrays' => [],
    		'hoa' => {},
    		'singletons' => {},
    		'find_centroid' => {},
    		'cluster_meta' => {},
    	},
    };
    bless $self, $class;

    # Get CPU Information
    $self->_getCpuInfo();

    # Calculate threads to use ( Lets use all atm )
    $self->{'config'}{'spawn_threads'} = $self->{'env'}{'cpu'}{'threads'};
	#print "Using ", $self->{'env'}{'cpu'}{'threads'}, " threads", "\n"; 

    # Return the object
    return $self;
}

# Internal Call State Subroutines

sub _clearState {
	my ($self) = @_;
    $self->{'exec-state'}->{'return_call_file'} = undef;
    $self->{'exec-state'}->{'return_call_package'} = undef;
    $self->{'exec-state'}->{'return_call_sub'} = undef;
    $self->{'exec-state'}->{'return_call_line'} = undef;
    $self->{'exec-state'}->{'return_state_code'} = undef;
    $self->{'exec-state'}->{'return_state_name'} = undef;
    $self->{'exec-state'}->{'return_state_reason'} = undef;
}
sub _setState {
	my ($self, $return_code, $return_state_reason) = @_;
	my ($package, $filename, $line, $subroutine, $hasargs) = caller(1);
	if ( $subroutine =~ qr/^setState/ ) {
		($package, $filename, $line, $subroutine, $hasargs) = caller(2);
	}
	my $return_state_name = "UNKNOWN";
	if ( $return_code eq '1' ) {
		$return_state_name = "SUCCESS";
	} elsif ( $return_code eq '2' ) {
		$return_state_name = "WARN";
	} elsif ( $return_code eq '3' ) {
		$return_state_name = "ERROR";
	} elsif ( $return_code eq '4' ) {
		$return_state_name = "FATAL";
	}
    $self->{'exec-state'}->{'return_call_file'} = $filename;
    $self->{'exec-state'}->{'vreturn_call_package'} = $package;
    $self->{'exec-state'}->{'return_call_sub'} = $subroutine;
    $self->{'exec-state'}->{'return_call_line'} = $line;
    $self->{'exec-state'}->{'return_state_code'} = $return_code;
    $self->{'exec-state'}->{'return_state_name'} = $return_state_name;
    $self->{'exec-state'}->{'return_state_reason'} = $return_state_reason;
}
sub _setStateSuccess {
	my $self = shift;
	push(@_, 1);
	$self->_setState(@_);
	return 1;
}
sub _setStateWarn {
	my $self = shift;
	push(@_, 2);
	$self->_setState(@_);
	return 1;
}
sub _setStateError {
	my $self = shift;
	push(@_, 3);
	$self->_setState(@_);
	return 0;
}
sub _setStateFatal {
	my $self = shift;
	push(@_, 4);
	$self->_setState(@_);
	return 0;
}

# Internal Helping Subroutines

sub _writeOut {
	my ($self, $message) = @_;

		return if ( $self->{'config'}{'quiet'} eq "on");

		my $s = Time::Date->now; 
    	print "Gregorian Time: ", $s, "\t", "Unix Time: "; 

		my ($seconds, $microseconds) = gettimeofday;
		printf("[%.3f]",$seconds + ($microseconds / 1000000));
		if ( defined($self->{'exec_state'}{'current_exec'}) && $self->{'exec_state'}{'current_exec'} ne '' ) {
			print $self->{'exec_state'}{'current_exec'};
		}
		print "\tStatus: ".$message."\n";
}

sub _getCpuInfoMac {
    my ($self) = @_;

    if ( $self->{'config'}{'threads'} =~  /^(0{0,2}[1-9]|0?[1-9][0-9]|[1-9][0-9][0-9])$/ ) {
        $self->{'env'}{'cpu'}{'threads'} = $self->{'config'}{'threads'};
    } else {
        my $procs = new Unix::Processors;
        $self->{'env'}{'cpu'}{'cores'}   = $procs->max_physical || 0;
        $self->{'env'}{'cpu'}{'threads'} = $procs->max_online   || 0;
    }

    die "Couldn't determine number of CPU threads on Mac"
        unless $self->{'env'}{'cpu'}{'threads'} > 0;
}

sub _getCpuInfoLinux {
	my ($self) = @_;
	my $cpufile = '/proc/cpuinfo';
	open(my $fh, '<:encoding(UTF-8)', $cpufile) or die "Could not open file /proc/cpuinfo? $!";
 
 	my $cpucount = 0;
 	my $threadcount = 0;
 	my $lastcoreid = undef;
 	my $cpuregex = qr/^\s*core id\s*:\s*([0-9]+)\s*$/;
	while (my $line = <$fh>) {
		if ( $line =~ $cpuregex ) {
			$threadcount++;
			if ( !defined($lastcoreid) || $lastcoreid ne $1 ) {
				$cpucount++;
				$lastcoreid = $1;
			}
		}
	}
	$self->{'env'}{'cpu'}{'cores'}   = $cpucount    || 0;
	$self->{'env'}{'cpu'}{'threads'} = $threadcount || 0;

	close($fh);

    die "Couldn't determine number of CPU threads on Linux"
        unless $self->{'env'}{'cpu'}{'threads'} > 0;
}

sub _getCpuInfo {
    my ($self) = @_;

    if ($^O =~ /linux/i) {
        return $self->_getCpuInfoLinux;
    }

    # Consider $^O is `darwin` by default
    return $self->_getCpuInfoMac;
}

sub _recreateClusterSizes {
	my ($self) = @_;
	$self->{'data'}{'hoa_size_lookup'} = {};
	foreach my $seed_centroid (keys %{$self->{'data'}{'hoa'}}){ 
		my $clustersize = scalar(@{$self->{'data'}{'hoa'}{$seed_centroid}});
		# Create our tag hash
		$self->{'data'}{'hoa_size_lookup'}{$seed_centroid} = $clustersize;
	}
}

# Internal Mathematic Subroutines

sub _smallest_value {
	my ($self, $hash) = @_;
	my @sorted = sort {$hash->{$a} <=> $hash->{$b}} keys %$hash; 
	return $sorted[0];
}

# Levenshtein distance calculation
sub _levdist {
	my ($self, $seq1, $seq2, $expected_number_of_errors, $tolerated_frame_shift) = @_;

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
	my $number_of_rows_before_diagonal = $self->{'config'}{'error_threshold'} + 1; 
	my $distance_threshold = $self->{'config'}{'error_threshold'} + $self->{'config'}{'tolerated_frame_shift'};
    if ($self->{'config'}{'error_threshold'} <= $l1){
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
    my $j_downstream_border = $self->{'config'}{'error_threshold'}*2+2;
    my $j_start = 1;
    my $j_end = $self->{'config'}{'error_threshold'}*2+1;

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
	# Calculate weighted errors 
	my $classic_distance = $dist->[$l1]->[$l2];
	return $classic_distance if ( $classic_distance == 0 ); # There is no difference between strings

	my $vert_min_dist;
	if ($self->{'config'}{'error_threshold'} > $l1){ 
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
			$vert_min_dist = $dist->[$l1-$self->{'config'}{'error_threshold'}]->[$l2];
			my $last_i = $l1-$self->{'config'}{'error_threshold'}; 
			my $start_i = $l1-1; 
			for my $i ($last_i..$start_i){ 
				$vert_min_dist = min($vert_min_dist, $dist->[$i]->[$l2]); 
			}
		} elsif ($l1 <= 1) { 
			$vert_min_dist = $dist->[$l1-$self->{'config'}{'error_threshold'}]->[$l2];
			my $last_i = $l1-$self->{'config'}{'error_threshold'};
            my $start_i = $l1;
            for my $i ($last_i..$start_i){
            	$vert_min_dist = min($vert_min_dist, $dist->[$i]->[$l2]);
            }
		}
	} 
	my $hor_min_dist; 
	if ($self->{'config'}{'error_threshold'} > $l2){
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
    } elsif ($self->{'config'}{'error_threshold'} <= $l2){
    	if ($l2 >= 1){
			$hor_min_dist = $dist->[$l1]->[$l2-$self->{'config'}{'error_threshold'}];	
            my $last_j = $l1-$self->{'config'}{'error_threshold'};
            my $start_j = $l2-1;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = min($hor_min_dist, $dist->[$l1]->[$j]);
            }
    	} elsif ($l2 <= 1) {
			$hor_min_dist = $dist->[$l1]->[$l2-$self->{'config'}{'error_threshold'}];
            my $last_j = $l1-$self->{'config'}{'error_threshold'};
            my $start_j = $l2;
            for my $j ($last_j..$start_j){
            	$hor_min_dist = min($hor_min_dist, $dist->[$l1]->[$j]);
            }
    	}
    }

	my $corner_min_dist = $classic_distance;
	my $modified_distance = min($vert_min_dist, $hor_min_dist, $corner_min_dist);

	my $sum; 
	my $bidir_sum; 

	# If the minimum distance is found on the vertical border...
	if (($vert_min_dist < $corner_min_dist && $vert_min_dist < $hor_min_dist)||($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $self->{'config'}{'weight_deletion'} >= $self->{'config'}{'weight_insertion'} )) { 		
		# Find the location of the modified minimum distance from the corner [i][j]
		for ( my $i = $l1-1; $i >= ($l1-($self->{'config'}{'error_threshold'}+1)); $i-- ) {
			if ($dist->[$i]->[$l2] == $modified_distance) {
				# If the modified distance is not consecutive
                if ($dist->[$i-1]->[$l2] != $dist->[$i]->[$l2]){
                    my $size_of_s2 = @s2; 
					my $distance_from_corner = ($size_of_s2 - $i); 
					# If only insertions are present
					if ($distance_from_corner == $modified_distance){ 
						$sum += ($self->{'config'}{'weight_insertion'}*$distance_from_corner); 
						$bidir_sum += ($self->{'config'}{'weight_deletion'}*$distance_from_corner); 
					# If insertions are mixed with other error types
					}elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_insertions = $distance_from_corner; 
						$sum += ($self->{'config'}{'weight_insertion'}*$base_number_of_insertions); 
						$bidir_sum += ($self->{'config'}{'weight_deletion'}*$base_number_of_insertions); 
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
						if ($remaining_steps % 2 == 0){
                            if ($self->{'config'}{'weight_substitution'} > ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                my $number_of_indel_pairs = $remaining_steps/2; 
							    $sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
							    $bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
                            }elsif ($self->{'config'}{'weight_substitution'} < ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            }elsif($self->{'config'}{'weight_substitution'} = ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            } 
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){ 
							# every remainder is a substitution. every dividend is an insertion-deletion pair. 
							use integer; 
							my $number_of_indel_pairs = $remaining_steps/2; 
							if ($number_of_indel_pairs == 0){ 
								$sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'}); 
							} else{ 
								$sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'})); 
								$bidir_sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'}));
							}
						} elsif ($remaining_steps == 1){ 
							$sum += ($self->{'config'}{'weight_substitution'}); 
							$bidir_sum += ($self->{'config'}{'weight_substitution'}); 
						}
					}	
					last; 
				# If the modified distance is consecutive
                } elsif ($dist->[$i-1]->[$l2] == $dist->[$i]->[$l2]) {
                    my $size_of_s2 = @s2; 
                	my $distance_from_corner = ($size_of_s2 - $i);
                    # If only insertions are present
					if ($distance_from_corner == $modified_distance){
                    	$sum += ($self->{'config'}{'weight_insertion'}*$distance_from_corner);
						$bidir_sum += ($self->{'config'}{'weight_deletion'}*$distance_from_corner); 
                    # If insertions are mixed with other error types
                    }elsif ($distance_from_corner < $modified_distance){
                    	my $base_number_of_insertions = $distance_from_corner;
                    	$sum += ($self->{'config'}{'weight_insertion'}*$base_number_of_insertions);
						$bidir_sum += ($self->{'config'}{'weight_deletion'}*$base_number_of_insertions); 
                        my $remaining_steps = ($modified_distance - $distance_from_corner);
						if ($remaining_steps % 2 == 0){
                        	if ($self->{'config'}{'weight_substitution'} > ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                my $number_of_indel_pairs = $remaining_steps/2; 
							    $sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
							    $bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
                            }elsif ($self->{'config'}{'weight_substitution'} < ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            }elsif($self->{'config'}{'weight_substitution'} = ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            } 
                        } elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){
                        	# every remainder is a substitution. every dividend is an insertion-deletion pair.
                            use integer;
                            my $number_of_indel_pairs = $remaining_steps/2;
					    	$sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'}));
							$bidir_sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'}));
                        } elsif ($remaining_steps == 1){
                            $sum += ($self->{'config'}{'weight_substitution'});
							$bidir_sum += ($self->{'config'}{'weight_substitution'}); 
                        }
                    }
                    last;
                }

	    	}
		} 
	# If the minimum distance is found on the horizontal border...
	} elsif (($hor_min_dist < $corner_min_dist && $hor_min_dist < $vert_min_dist)||($hor_min_dist == $vert_min_dist && $hor_min_dist < $corner_min_dist && $self->{'config'}{'weight_insertion'} >= $self->{'config'}{'weight_deletion'})){
		# Find the location of the modified minimum distance from corner [i][j]
		for ( my $j = $l2-1; $j >= ($l2-($self->{'config'}{'error_threshold'}+1)); $j-- ){ 
			if ( $dist->[$l1]->[$j] == $modified_distance) {
				# If the modified distance is not consecutive
                if ($dist->[$l1]->[$j] != $dist->[$l1]->[$j-1]){
                    my $size_of_s2 = @s2; 
					my $distance_from_corner = ($size_of_s2 - $j); 
					# If only deletions are present 
					if ($distance_from_corner == $modified_distance){ 
						$sum += (($self->{'config'}{'weight_deletion'}*$distance_from_corner)); 
						$bidir_sum += (($self->{'config'}{'weight_insertion'}*$distance_from_corner)); 
					# If deletions are mixed with other error types 
					} elsif ($distance_from_corner < $modified_distance){ 
						my $base_number_of_deletions = $distance_from_corner; 
						$sum += (($self->{'config'}{'weight_deletion'}*$base_number_of_deletions)); 
						$bidir_sum += (($self->{'config'}{'weight_insertion'}*$base_number_of_deletions)); 
						my $remaining_steps = ($modified_distance - $distance_from_corner); 
						if ($remaining_steps %2 == 0){ 
							if ($self->{'config'}{'weight_substitution'} > ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                my $number_of_indel_pairs = $remaining_steps/2; 
							    $sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
							    $bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
                            }elsif ($self->{'config'}{'weight_substitution'} < ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            }elsif($self->{'config'}{'weight_substitution'} = ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            } 
						} elsif ($remaining_steps == 1){ 
							$sum += ($self->{'config'}{'weight_substitution'});
							$bidir_sum += ($self->{'config'}{'weight_substitution'});
						} elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){
                            # every remainder is a substitution. every dividend is an insertion-deletion pair.
                            use integer;
                            my $number_of_indel_pairs = $remaining_steps/2;
                            $sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'})); 
							$bidir_sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'})); 
						} 
					} 
                    last; 
				# If the modified distance is consecutive
                } elsif ( $dist->[$l1]->[$j] == $dist->[$l1]->[$j-1]) {
                    my $size_of_s2 = @s2; 
                    my $distance_from_corner = ($size_of_s2 - $j);
                    # If only deletions are present
                    if ($distance_from_corner == $modified_distance){
                    	$sum += (($self->{'config'}{'weight_deletion'}*$distance_from_corner));
						$bidir_sum += (($self->{'config'}{'weight_insertion'}*$distance_from_corner)); 
                    # If deletions are mixed with other error types
                    } elsif ($distance_from_corner < $modified_distance){
                        my $base_number_of_deletions = $distance_from_corner;
                        $sum += (($self->{'config'}{'weight_deletion'}*$base_number_of_deletions));
						$bidir_sum += (($self->{'config'}{'weight_insertion'}*$base_number_of_deletions)); 
                        my $remaining_steps = ($modified_distance - $distance_from_corner);
                        if ($remaining_steps %2 == 0){
                            if ($self->{'config'}{'weight_substitution'} > ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                my $number_of_indel_pairs = $remaining_steps/2; 
							    $sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
							    $bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs)); 
                            }elsif ($self->{'config'}{'weight_substitution'} < ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            }elsif($self->{'config'}{'weight_substitution'} = ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'})){ 
                                $sum += ($self->{'config'}{'weight_substitution'}); 
								$bidir_sum += ($self->{'config'}{'weight_substitution'});
                            } 
                        } elsif ($remaining_steps == 1){
                            $sum += ($self->{'config'}{'weight_substitution'});
							$bidir_sum += ($self->{'config'}{'weight_substitution'}); 
                        } elsif (($remaining_steps > 1) && ($remaining_steps % 2 == 1)){
                        # every remainder is a substitution. every dividend is an insertion-deletion pair.
                            use integer;
                            my $number_of_indel_pairs = $remaining_steps/2;
                            $sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'}));
							$bidir_sum += (((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * $number_of_indel_pairs) + $self->{'config'}{'weight_substitution'}));
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
				if ($self->{'config'}{'weight_substitution'} > ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_insertion'})){
					$sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * ($corner_min_dist/2))); 
					$bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * ($corner_min_dist/2))); 
				} elsif ($self->{'config'}{'weight_substitution'} < ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_insertion'})){
                    $sum += (($self->{'config'}{'weight_substitution'}*$corner_min_dist)); 
					$bidir_sum += (($self->{'config'}{'weight_substitution'}*$corner_min_dist)); 
                } elsif ($self->{'config'}{'weight_substitution'} == ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_insertion'})){
                    if ($corner_min_dist > 2){ 
						my $half_corner_value = ($corner_min_dist/2);
						$sum += (($self->{'config'}{'weight_substitution'}*$half_corner_value)); 
						$bidir_sum += (($self->{'config'}{'weight_substitution'}*$half_corner_value));
						if ($half_corner_value % 2 == 0){ 
							$sum += ( ($self->{'config'}{'weight_insertion'}*($half_corner_value/2)) + ($self->{'config'}{'weight_deletion'}*($half_corner_value/2)) ); 
							$bidir_sum += ( ($self->{'config'}{'weight_insertion'}*($half_corner_value/2)) + ($self->{'config'}{'weight_deletion'}*($half_corner_value/2)) );
						} else { 
							if ($self->{'config'}{'weight_insertion'} <= $self->{'config'}{'weight_deletion'}){ 
								$sum += ( ($self->{'config'}{'weight_insertion'} * ceil($half_corner_value/2)) + ($self->{'config'}{'weight_deletion'} * floor($half_corner_value/2)) ); 
								$bidir_sum += ( ($self->{'config'}{'weight_insertion'} * ceil($half_corner_value/2)) + ($self->{'config'}{'weight_deletion'} * floor($half_corner_value/2)) );
							}else{ 
								$sum += ( ($self->{'config'}{'weight_insertion'} * floor($half_corner_value/2)) + ($self->{'config'}{'weight_deletion'} * ceil($half_corner_value/2)) ); 
								$bidir_sum += ( ($self->{'config'}{'weight_insertion'} * floor($half_corner_value/2)) + ($self->{'config'}{'weight_deletion'} * ceil($half_corner_value/2)) );
							}
						}
					}elsif ($corner_min_dist == 2){
						if ($self->{'config'}{'weight_insertion'} <= $self->{'config'}{'weight_deletion'}){ 
							$sum += ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_substitution'}); 
							$bidir_sum += ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_substitution'}); 
						}else { 
							$sum += ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_substitution'}); 
							$bidir_sum += ($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_substitution'}); 
						}
					}
                }
        	}
        	# If the number at [i][j] is odd, calculate the number of indels and subs needed. Preference for lowest weights. 
        	if ($corner_min_dist % 2 == 1){
            	if ($self->{'config'}{'weight_substitution'} > ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_insertion'})){
                	$sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * (($corner_min_dist/2) - 0.5)) + $self->{'config'}{'weight_substitution'});
					$bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * (($corner_min_dist/2) - 0.5)) + $self->{'config'}{'weight_substitution'});
                } elsif ($self->{'config'}{'weight_substitution'} < ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_insertion'})){
					$sum += (($self->{'config'}{'weight_substitution'} * $corner_min_dist));
					$bidir_sum += (($self->{'config'}{'weight_substitution'} * $corner_min_dist));
                } elsif ($self->{'config'}{'weight_substitution'} == ($self->{'config'}{'weight_deletion'} + $self->{'config'}{'weight_insertion'})){
                    $sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * (($corner_min_dist/2) - 0.5)) + $self->{'config'}{'weight_substitution'});
					$bidir_sum += ((($self->{'config'}{'weight_insertion'} + $self->{'config'}{'weight_deletion'}) * (($corner_min_dist/2) - 0.5)) + $self->{'config'}{'weight_substitution'});
                }
        	}
		} elsif ($corner_min_dist == 1){
       		# If there is only one edit distance at [i][j] it is a substitution
        	$sum += ($self->{'config'}{'weight_substitution'});
			$bidir_sum += ($self->{'config'}{'weight_substitution'}); 
		}
	}
	
	# If we have nothing, just return.
	return if (!defined($sum));

	# Return
	return ($sum, $bidir_sum);
}

# Internal Parallelization Subroutines
sub _createSequenceCluster_Process {
	my ($self, $dataset_tag) = @_;
	my $is_singleton = 0;
	my $average_intracluster_distance;
    my @cluster_distances; 
    my @clustered_sequences; 

    #print $dataset_tag, "\n"; 

	foreach my $other_tag (@{$self->{'data'}{'tag_array'}}) { 
		if ($other_tag eq $dataset_tag){ 
            # If the two tags match, no need for calculations before populating data structures
			push(@clustered_sequences, $dataset_tag);
			push(@cluster_distances, 0);
		} elsif ($other_tag !~ $dataset_tag){
            if ($self->{'config'}{'weight_insertion'} == $self->{'config'}{'weight_deletion'}){
			    my ($dist_1, $dist_2) = $self->_levdist($other_tag, $dataset_tag, $self->{'config'}{'error_threshold'}, $self->{'config'}{'tolerated_frame_shift'}); 
				# Distance returns undefined if the LD matrix table is aborted due to LD > error + frameshift to save time
				if ((defined $dist_1) && (defined $dist_2)){ 	
				    if ( ($dist_1 < ($self->{'config'}{'weight_threshold'}+1)) && ($dist_2 < ($self->{'config'}{'weight_threshold'}+1)) ){
					    push(@clustered_sequences, $other_tag);
						my $min_dist = min $dist_1, $dist_2; 
						push(@cluster_distances, $min_dist); 
					}
				}
			} else { 
                my $rev_dataset_tag = reverse $dataset_tag; 
				my $rev_other_tag = reverse $other_tag; 
				my ($dist_1, $dist_2) = $self->_levdist($dataset_tag, $other_tag, $self->{'config'}{'error_threshold'}, $self->{'config'}{'tolerated_frame_shift'}); 
				my ($rev_dist_1, $rev_dist_2) = $self->_levdist($rev_dataset_tag, $rev_other_tag, $self->{'config'}{'error_threshold'}, $self->{'config'}{'tolerated_frame_shift'});
			    if ((defined $dist_1) && (defined $dist_2) && (defined $rev_dist_1) && (defined $rev_dist_2)){ 
					if (($dist_1 < ($self->{'config'}{'weight_threshold'}+1)) && ($dist_2 < ($self->{'config'}{'weight_threshold'}+1)) && ($rev_dist_1 < ($self->{'config'}{'weight_threshold'}+1)) && ($rev_dist_2 < ($self->{'config'}{'weight_threshold'}+1))){
						push(@clustered_sequences, $other_tag);
						my $min_dist = min $dist_1, $dist_2, $rev_dist_1, $rev_dist_2; 
						push(@cluster_distances, $min_dist); 
					}
				}
            }  
		}
	} 
	# Indicate if the dataset tag is a singleton
	shift @cluster_distances;
	if (scalar(@cluster_distances) == 0) {
		$is_singleton = 1;
	# If it is not a singleton, collect the average edit distance between the seed centroid and all other sequences in the cluster
	} else {
		my $sum = reduce {$a + $b} @cluster_distances; 
		$average_intracluster_distance = $sum/scalar(@cluster_distances);
	}

	# We are cleaning up as CHILD by doing this -- Rest of FOR block should be considered CHILD
	my %result = (
		'is_singleton' => $is_singleton,
		'dataset_tag' => $dataset_tag,
		'clustered_sequences' => \@clustered_sequences,
		'average_intracluster_distance' => $average_intracluster_distance
	);

	return \%result;
}

##
# Internal Methods, listed in the expected call order
##

# Load data for processing
sub loadDataset {
	my ($self,$loadFile, $subloadFile) = @_;
    $self->{'exec_state'}{'current_exec'} = "\tCurrent State: Loading Dataset     ";
    $self->_writeOut("Loading Dataset");
	
	# Have we already loaded a dataset?
	if ( defined($self->{'config'}->{'input_dataset_file'}) && $self->{'config'}->{'input_dataset_file'} eq '' ) {
		return $self->_setStateFatal("Attempted to loadDataset after one has been loaded, We cannot do this. Create a new ThreeGold instance to load a new dataset");
	}
	# Do we have a loadFile
	if ( !defined($loadFile) || $loadFile eq '' ) {
		return $self->_setStateFatal("No loadFile provided to loadDataset, Cannot proceed.");
	}
	# Set our large Input Dataset
	$self->{'config'}->{'input_dataset_file'} = $loadFile;

	# Set our sub input dataset
	$self->{'config'}->{'sub_input_dataset_file'} = $subloadFile;

	# Can we access it?
	if ( ! -r $self->{'config'}->{'input_dataset_file'} ) {
		return $self->_setStateFatal("processDatasetFile Could not find read file \"".$self->{'config'}->{'input_dataset_file'}."\"! Does this file exist?");
	}	
	if ( ! -r $self->{'config'}->{'sub_input_dataset_file'} ) {
		return $self->_setStateFatal("processDatasetFile Could not find read file \"".$self->{'config'}->{'sub_input_dataset_file'}."\"! Does this file exist?");
	}
	# Lets open our sources
	if ( ! open(DATASET, '<', $self->{'config'}->{'input_dataset_file'}) ) {
		return $self->_setStateFatal("processDatasetFile attempted to open file \"".$self->{'config'}->{'input_dataset_file'}."\" but received error ".$!);
	} 
	if ( ! open(SUBDATASET, '<', $self->{'config'}->{'sub_input_dataset_file'}) ) {
		return $self->_setStateFatal("processDatasetFile attempted to open file \"".$self->{'config'}->{'sub_input_dataset_file'}."\" but received error ".$!);
	} 
    # Parse our large data
	my $current_header;
	while (my $tag = <DATASET>){ 
        chomp $tag; 
        if ($tag =~ /(^>)(.*)/){
		    $current_header = $2;   
	    } elsif ($tag !~ /^>/){ 
		    # Remove leading and trailing blank spaces
		    $tag =~ s/^\s+|\s+$//g;
            # Skip blank tags and # prefixed
		    next if ( $tag eq '' || $tag =~ qr/^#/ );
            # Create an array of all tags including duplicates
		    push(@{$self->{'data'}{'tag_array'}}, $tag);
		}  
    }
	close(DATASET);

	# Parse our subsection data
	while (my $tag = <SUBDATASET>){ 
        chomp $tag; 
        if ($tag =~ /(^>)(.*)/){
		    $current_header = $2;   
	    } elsif ($tag !~ /^>/){ 
		    # Remove leading and trailing blank spaces
		    $tag =~ s/^\s+|\s+$//g;
            # Skip blank tags and # prefixed
		    next if ( $tag eq '' || $tag =~ qr/^#/ );
            # Create an array of all tags including duplicates
		    push(@{$self->{'data'}{'sub_tag_array'}}, $tag);
		}  
    }
	close(SUBDATASET);

	return $self->_setStateSuccess();
}

# Create distance matrix clusters and locate singletons
sub createSequenceClusters {
	my ($self) = @_;
	$self->{'exec_state'}{'current_exec'} = "\tCurrent State: Building Clusters";
	$self->_writeOut("Starting Cluster Generation");
	
	# Create thread pool
	my %threadinfo;
	my $threadpool = Thread::Pool->new({
  		optimize => 'cpu', # default: 'memory'
        pre => sub {shift; $self->_writeOut("Starting Worker Thread ID ".(threads->tid()-1)); },
  		do => sub {
  			shift;
  			$threadinfo{threads->tid()}++;
  			if ( ( (threads->tid()-1) + ($threadinfo{threads->tid()}-1) ) % 1  == 0 ) {
                my $percent_completed = ( ( ($threadinfo{threads->tid()}-1) / scalar @{$self->{'data'}{'sub_tag_array'}} ) * 100 );
				my $percent_completed_actual = ($percent_completed * $self->{'env'}{'cpu'}{'threads'});
  				$self->_writeOut("[Thread ID:".(threads->tid()-1)."] Built ".($threadinfo{threads->tid()}-1)." Clusters, ". sprintf("%.2f", $percent_completed_actual). "% Done");
  			} 
            # Every single tag is passed into this subroutine, one at a time. 
  			$self->_createSequenceCluster_Process(@_);
  		},
  		post => sub {shift; $self->_writeOut("Stopping Worker TID ".(threads->tid()-1));},

  		workers => $self->{'config'}{'spawn_threads'},
  		maxjobs => $self->{'config'}{'spawn_threads'}
 	});

	# Send in our data to run
	my @jobids;
    foreach my $dataset_tag (@{$self->{'data'}{'sub_tag_array'}}) {
       # Every single tag is put into this array as well, one at a time. 
		# Put work into our queue
		my @params = ( $self, $dataset_tag );
		push(@jobids,$threadpool->job(@params));
	}
	# Wait for our threads to complete and join
	$threadpool->shutdown;

	# store our data
	$self->_writeOut("Extracting Information About The Clusters");
	foreach my $jobid (@jobids) {
		# $retdata = "returned data from the cluster process subroutine"
		my $retdata = $threadpool->result($jobid);
		if ( $retdata->{'is_singleton'} ) {
			# We have a singleton, Add an entry into our singletons keyset
			$self->{'data'}{'singletons'}{$retdata->{'dataset_tag'}} = 1;
		}
		# If we don't have a singleton, populate hashes
		# A singleton would have an undefined average cluster distance
		if (defined $retdata->{'average_intracluster_distance'}){ 
			$self->{'data'}{'find_centroid'}{$retdata->{'dataset_tag'}} = $retdata->{'average_intracluster_distance'};
			$self->{'data'}{'hoa'}{$retdata->{'dataset_tag'}} = $retdata->{'clustered_sequences'};
		}
	}

	$self->_recreateClusterSizes();

    open (SINGLETON_OUT, '>', $self->{'config'}{'singletons_out'}) or die $!; 
	open (DISTANCE_OUT, '>', $self->{'config'}{'ave_cluster_distances_out'}) or die $!; 
    open (CLUSTER_OUT, '>', $self->{'config'}{'clusters_out'}) or die $!; 
	open (SIZE_OUT, '>', $self->{'config'}{'size_out'}) or die $!; 

	my @cluster_out_array; 
	while(my ( $key, $val) = each %{$self->{'data'}{'hoa'}}){ 
		push(@cluster_out_array, "$key\t=>\t@$val\n"); 
	}

	my @singleton_out_array; 
	while(my ( $key, $val) = each %{$self->{'data'}{'singletons'}}){ 
		push(@singleton_out_array, "$key\t=>\t$val\n"); 
	}

	my @distance_out_array; 
	while(my ( $key, $val) = each %{$self->{'data'}{'find_centroid'}}){ 
		push(@distance_out_array, "$key\t=>\t$val\n"); 
	}

	my @size_lookup_out_array; 
	while(my ( $key, $val) = each %{$self->{'data'}{'hoa_size_lookup'}}){ 
		push(@size_lookup_out_array, "$key\t=>\t$val\n"); 
	}

    print SINGLETON_OUT @singleton_out_array; 
    print DISTANCE_OUT @distance_out_array; 
    print CLUSTER_OUT @cluster_out_array; 
	print SIZE_OUT @size_lookup_out_array; 

	$self->_writeOut("Initial Cluster Creation Completed");
	$self->{'exec_state'}{'current_exec'} = undef;
}

# End positively to load module
1;
