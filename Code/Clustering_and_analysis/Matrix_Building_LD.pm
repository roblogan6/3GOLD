package Matrix_Building_LD;

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

# Levenshtein distance calculation
sub _levdist {
	my ($self, $seq1, $seq2) = @_;

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
    my $classic_distance;

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
	        $distances->[$i]->[$j] = min(($distances->[$i-1]->[$j-1] + $cost), ($distances->[$i]->[$j-1]+1), ($distances->[$i-1]->[$j]+1) );
                $classic_distance = $distances->[$i]->[$j]; 
            }
    }
    return $classic_distance; 
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
            #my $rev_dataset_tag = reverse $dataset_tag; 
			#my $rev_other_tag = reverse $other_tag; 
			my $dist_1 = $self->_levdist($other_tag, $dataset_tag); 
			#my $dist_2 = $self->_levdist($rev_dataset_tag, $rev_other_tag); 
			# Distance returns undefined if the LD matrix table is aborted due to LD > error + frameshift to save time
			if ((defined $dist_1)){ # && (defined $dist_2)){ 	
				if ( ($dist_1 < ($self->{'config'}{'error_threshold'}+1))){
                    #if ( ($dist_2 < ($self->{'config'}{'error_threshold'}+1))){
					    push(@clustered_sequences, $other_tag);
						#my $min_dist = min $dist_1, $dist_2; 
						#push(@cluster_distances, $min_dist); 
						push(@cluster_distances, $dist_1); 
					#}
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
  				$self->_writeOut("[Thread ID:".(threads->tid()-1)."] Built ".($threadinfo{threads->tid()}-1)." Clusters, ". $percent_completed_actual. "% Done");
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
