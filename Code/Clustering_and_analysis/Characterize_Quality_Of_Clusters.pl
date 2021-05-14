#!/usr/bin/perl
use warnings; 
use strict; 
use List::Compare; 
use Math::Round; 
use List::Util qw( max ); 

# Usuage = < experimental cluster in fasta and tab format >, < reference cluster in the same format >, < the size below which an experiemnt cluster won't be compared to ref clusters >, < The size of the reference clusters >, <output>
my @output; 
my $count = 0;

open(EXP, '<', $ARGV[0]) or die $!; 
while (my $exp_line = <EXP>){ 
    chomp $exp_line; 
    if ($exp_line =~ /(^>)/){ 
        next; 
    #}elsif ($exp_line =~ /((^[\-\w]{36})(.*)([\-\w]{36}$))/gm){
    }elsif ($exp_line =~ /((^[\-\w]{20})(.*)([\-\w]{20}$))/gm){
        my $exp_cluster = $1; 
        my $exp_first_element = $2;
        my $exp_last_element = $4; 
        my @array_elements = split /\t/, $exp_line; 
        my $number_of_array_elements = @array_elements;  
        my $middle_array_element = round ($number_of_array_elements/2);
        my $exp_middle_element = $array_elements[$middle_array_element]; 
        my $first_sensitivity; 
        my $first_specificity; 
        my $last_sensitivity; 
        my $last_specificity; 
        my $middle_sensitivity; 
        my $middle_specificity; 
        my @first_element_cache; 
        my @last_element_cache; 

        my @exp_cluster = split /\t/, $exp_cluster; 
        my $size_of_exp_cluster = @exp_cluster; 
        if ($size_of_exp_cluster < $ARGV[2]){ 
            next; 
        }

        open(REF, '<', $ARGV[1]) or die $!; 
        while (my $ref_line = <REF>){ 
            if ($ref_line =~ /^>/){ 
                next; 
            } elsif ($ref_line =~ /(.*)($exp_last_element|$exp_first_element|$exp_middle_element)(.*)/){
                if ($ref_line =~ /(.*)($exp_last_element)(.*)/){ 
                    my $ref_cluster = $1.$2.$3;
                    my @ref_cluster = split /\t/, $ref_cluster;
                    $count++; 
                    my $lc = List::Compare->new(\@exp_cluster, \@ref_cluster);
                    # Find the sequences which are only found in the reference cluster (sensitivity test)
                    # Sensitivity is the ability to correctly identify those who belong to a group 
                    my @ref_only = $lc->get_complement;
                    my $number_of_ref_only = @ref_only;
                    $last_sensitivity = (100-(100*($number_of_ref_only/$ARGV[3]))); 
                    # Find the sequences which are only found in the experimental cluster (specificity test) 
                    # Specificity is the ability to correctly identify those who do not belong to a group 
                    my @exp_only = $lc->get_unique;
                    my $number_of_exp_only = @exp_only;
                    $last_specificity = (100*($number_of_exp_only/$size_of_exp_cluster)); 
                }
                if ($ref_line =~ /(.*)($exp_first_element)(.*)/){ 
                    my $ref_cluster = $1.$2.$3; 
                    my @ref_cluster = split /\t/, $ref_cluster; 
                    $count++; 
                    my $lc = List::Compare->new(\@exp_cluster, \@ref_cluster);
                    # Find the sequences which are only found in the reference cluster (sensitivity test)
                    # Sensitivity is the ability to correctly identify those who belong to a group 
                    my @ref_only = $lc->get_complement; 
                    my $number_of_ref_only = @ref_only; 
                    $first_sensitivity = (100-(100*($number_of_ref_only/$ARGV[3]))); 
                    # Find the sequences which are only found in the experimental cluster (specificity test) 
                    # Specificity is the ability to correctly identify those who do not belong to a group 
                    my @exp_only = $lc->get_unique; 
                    my $number_of_exp_only = @exp_only;
                    $first_specificity = (100*($number_of_exp_only/$size_of_exp_cluster)); 
                }
                if ($ref_line =~ /(.*)($exp_middle_element)(.*)/){ 
                    my $ref_cluster = $1.$2.$3; 
                    my @ref_cluster = split /\t/, $ref_cluster; 
                    $count++; 
                    my $lc = List::Compare->new(\@exp_cluster, \@ref_cluster);
                    # Find the sequences which are only found in the reference cluster (sensitivity test)
                    # Sensitivity is the ability to correctly identify those who belong to a group 
                    my @ref_only = $lc->get_complement; 
                    my $number_of_ref_only = @ref_only; 
                    $middle_sensitivity = (100-(100*($number_of_ref_only/$ARGV[3]))); 
                    # Find the sequences which are only found in the experimental cluster (specificity test) 
                    # Specificity is the ability to correctly identify those who do not belong to a group 
                    my @exp_only = $lc->get_unique; 
                    my $number_of_exp_only = @exp_only;
                    $middle_specificity = (100*($number_of_exp_only/$size_of_exp_cluster)); 
                } 

                $first_sensitivity = 0 unless defined $first_sensitivity; 
                $last_sensitivity = 0 unless defined $last_sensitivity; 
                $middle_sensitivity = 0 unless defined $middle_sensitivity; 

                # Find the highest sensitivity value 
                my @sensitivity_array; 
                push(@sensitivity_array, $first_sensitivity, $last_sensitivity, $middle_sensitivity); 
                my $max_sensitivity = max @sensitivity_array; 
                if ($max_sensitivity > 19){ 
                    if ($max_sensitivity == $first_sensitivity){ 
                        push(@output, "Percent of the correct sequences included: \t", $first_sensitivity, "\t");      
                        push (@output, "Percent of incorrect sequences included: \t", $first_specificity, "\t", "total clusters included in cluster: \t", $size_of_exp_cluster, "\n"); 
                    }elsif ($max_sensitivity == $last_sensitivity){ 
                        push(@output, "Percent of the correct sequences included: \t", $last_sensitivity, "\t");      
                        push (@output, "Percent of incorrect sequences included: \t", $last_specificity, "\t", "total clusters included in cluster: \t", $size_of_exp_cluster, "\n"); 
                    }elsif ($max_sensitivity == $middle_sensitivity){ 
                        push(@output, "Percent of the correct sequences included: \t", $middle_sensitivity, "\t");      
                        push (@output, "Percent of incorrect sequences included: \t", $middle_specificity, "\t", "total clusters included in cluster: \t", $size_of_exp_cluster, "\n"); 
                    }
                } 
                last if ($max_sensitivity > 19);  
            } 
        }
        close REF; 
    } 
   next; 
}
open(OUT, '>', "$ARGV[4]") or die $!; 
print OUT @output; 