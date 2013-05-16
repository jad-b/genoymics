# class_level_percentage.pl
# Created for Oy* team of MCBS 913, Spring '13 at UNH
# Created by Jeremy Dobbins-Bucklad
# Built and tested on Perl v5.14.2 on x64 linux
#
# Outputs statistics on taxonomic assignments.
#
# Input: .txt file containing taxonomic assignments. Format:
# OTU#  k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus;s__ 0.900
#
# Output: % OTUs classified at each level, w/ filename column-header. Tabs
# are used as delimiters.
#   filename
#   Kingdom %%  averaged confidence value
#   Phylum  %%  averaged confidence value
#   etc...
#
# Sample output:
# tutorial_tax_assns.txt
# Kingdom 3.5799522673031 0.991333333333333
# Phylum  0.715990453460621   0.886666666666667
# Class   2.86396181384248    0.93
# Order   4.29594272076372    0.906666666666667
# Family  65.6324582338902    0.955999999999999
# Genus   16.7064439140811    0.947428571428571
# Species 6.20525059665871    0.952307692307692

use strict;
use warnings;

my $debug = 0;

#filenames
my $taxFile = $ARGV[0];
#arrays
my @classifs=(0,0,0,0,0,0,0);   # classification levels
my @confs;                      # confidence levels
my @names = ("Kingdom","Phylum","Class","Order","Family","Genus","Species");
#variables
my $lines = 0;                  # number of OTUs read in
#regexes
my $classPtrn = "(.)__(.+)";    # letter followed by '__' followed by something



open ( IN, $taxFile) or die "Unable to open ".$taxFile;

my $input;
while( $input = <IN> ){
    $lines++;

    chomp $input;
    my @line = split ("\t",$input); #split into OTU#, classifications, confidence#
    if( @line != 3 ){               # Should be three values
        die "format error: line ".$lines." reads '".$input."'\n";
    }

    #find highest level of classification
    my @classes = split (";",$line[1]); # split by level
    if( @classes < 1 ){  #Something's wrong, nothing's there!
        die "format error: line ".$lines." contains no taxonomic data\n";
    }

    #Iterate species->Kingdom
    my $match = "";
    for (my $i=@classes; $i >= 0; $i-- ){
       if( defined $classes[$i] && $classes[$i] =~ /$classPtrn/ ){ # if defined and a match
            $classifs[$i]++;    # Up count for that classification level
            $confs[$i] += $line[2]; # Store conf level for later processing
            if( $debug > 1 ){print STDERR "$lines: $& @ $confs[$i]\n";}
            last;  # Found highest level of classification; quit looking.
        }
    }

}

# Both arrays should be equal in size; if not, we have an error in our logic above.
if( @classifs != @confs ){
    die "Number of classifications and confidence levels do not match\n";
}

if( $debug ){
    print STDERR "lines:\t$lines\nclassifs:\t@classifs\nconfs:\t@confs\n";   # Output both arrays
}


#At this point, @classifs holds a tally for the number of dead-ends at each level,
#and @confs has their summed confidence level. We need to average the confidence levels,
#and turn the tallies into percentages.
# Process values
for( my $i=0; $i < @classifs; $i++ ){
    $confs[$i] = $confs[$i] / $classifs[$i]; # Average confidence values
    $classifs[$i] = ( $classifs[$i] / $lines ) * 100.0;
}


#Output
print "$taxFile\n";
for( my $i=0; $i < @classifs; $i++ ){
    print "$names[$i]\t$classifs[$i]\t$confs[$i]\n";
}

if($debug){
    my $total = 0;
    $total += $_ for @classifs;
    print STDERR "total %:\t$total\n";
}

close $taxFile;
